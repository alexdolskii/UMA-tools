"""Owned temporary storage and supervision for UMA CLI processes.

No imaging libraries are imported here. The supervisor waits for Python
and JVM shutdown before cleanup, so their file handles are closed.
"""

import functools
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

import psutil

SCHEMA = "uma-runtime-1"
RUN_PATTERN = re.compile(r"^\d{8}T\d{12}Z_[0-9a-f]{12}$")
COMMANDS = {
    "uma_alignment": "alignment",
    "uma_thickness": "thickness",
    "area_analysis": "area",
    "uma_collect_results": "collect_results",
    "uma_report": "report",
}
_CLEANUP_BLOCKS = []


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def runtime_home():
    value = Path(
        os.environ.get("UMA_RUNTIME_HOME", Path.home() / ".uma-tools")
    ).expanduser()
    if not value.is_absolute():
        raise ValueError("UMA_RUNTIME_HOME must be an absolute path")
    return value.resolve()


def write_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + ".writing")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def identity(process=None):
    try:
        process = process or psutil.Process()
        return {"pid": process.pid, "created": process.create_time()}
    except psutil.Error as error:
        return {
            "pid": process.pid if process else os.getpid(),
            "created": None,
            "error": str(error),
        }


def observation_supported():
    """
    A remapped/restricted proc filesystem cannot establish process
    ownership.
    """
    try:
        if sys.platform.startswith("linux"):
            with open("/proc/self/stat", encoding="utf-8") as handle:
                if int(handle.read().split(" ", 1)[0]) != os.getpid():
                    return False
        return identity().get("created") is not None
    except (OSError, ValueError):
        return False


def process_state(record, host=None):
    """
    Do not confuse a reused PID, inaccessible process, or another host
    with exit.
    """
    if not isinstance(record, dict) or not observation_supported():
        return "UNKNOWN"
    if host and host != platform.node():
        return "UNKNOWN"
    if record.get("created") is None:
        return "UNKNOWN"
    try:
        process = psutil.Process(int(record["pid"]))
        if abs(process.create_time() - float(record["created"])) > 0.01:
            return "EXITED"
        if process.status() == psutil.STATUS_ZOMBIE:
            return "EXITED"
        return "ACTIVE"
    except psutil.NoSuchProcess:
        return "EXITED"
    except (psutil.AccessDenied, OSError, ValueError, KeyError, TypeError):
        return "UNKNOWN"


def current_run():
    run_id = os.environ.get("UMA_RUN_ID", "")
    directory = os.environ.get("UMA_RUN_DIR")
    if not directory or not RUN_PATTERN.fullmatch(run_id):
        return None
    return run_id, Path(directory)


def runtime_context():
    run = current_run()
    if run:
        run_id, directory = run
        return (
            f"runtime_run_id={run_id} | registry={directory / 'run.json'} | "
            f"temporary_directory={os.environ.get('TMPDIR', '')}"
        )
    return ""


def record_resource(directory, path):
    # Append once per owned directory, not for every image/checkpoint.
    payload = (
        json.dumps({"path": str(path), "kind": "temporary_directory"}) + "\n"
    )
    with (directory / "resources.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(payload)


def create_owned_directory(path, run_id, directory):
    path = Path(path)
    try:
        path.mkdir(mode=0o700)
    except FileExistsError:
        if path.is_symlink() or not path.is_dir():
            raise OSError(f"Unsafe temporary directory: {path}")
        owner = json.loads(
            (path / ".uma-owner.json").read_text(encoding="utf-8")
        )
        if owner != {"schema": SCHEMA, "run_id": run_id}:
            raise OSError(
                f"Temporary directory belongs to a different run: {path}"
            )
    else:
        # Register first: an interrupted creation remains visible to
        # diagnostics.
        record_resource(directory, path)
        write_json(
            path / ".uma-owner.json", {"schema": SCHEMA, "run_id": run_id}
        )
    return path


def temporary_path(destination, legacy=None):
    """
    Unique staging file on the destination filesystem; never move final
    outputs.
    """
    destination = Path(destination)
    run = current_run()
    if run is None:
        return (
            Path(legacy)
            if legacy is not None
            else destination.with_name("." + destination.name + ".tmp")
        )
    run_id, directory = run
    area = create_owned_directory(
        destination.parent.resolve() / (".uma_tmp_" + run_id),
        run_id,
        directory,
    )
    return area / (uuid4().hex + "_" + destination.name)


def resource_paths(directory, limit=100000):
    path = Path(directory) / "resources.jsonl"
    if not path.exists():
        return [], []
    paths, errors = [], []
    with path.open(encoding="utf-8") as handle:
        for index, line in enumerate(handle, 1):
            if index > limit:
                errors.append(f"Resource registry entry limit reached: {path}")
                break
            try:
                value = json.loads(line)
                if value["kind"] != "temporary_directory":
                    raise ValueError("unknown resource kind")
                item = Path(value["path"])
                if not item.is_absolute():
                    raise ValueError("relative resource path")
                if item not in paths:
                    paths.append(item)
            except (ValueError, KeyError, TypeError) as error:
                errors.append(f"{path}:{index}: {error}")
    return paths, errors


def owned_path(path, record):
    """
    Validate both the namespace and owner marker before inspecting or
    deleting.
    """
    run_id = record.get("run_id", "")
    if not RUN_PATTERN.fullmatch(run_id):
        return False
    if path.name not in (run_id, ".uma_tmp_" + run_id):
        return False
    if path.is_symlink() or path.resolve() != path or not path.is_dir():
        return False
    try:
        marker = path / ".uma-owner.json"
        return not marker.is_symlink() and json.loads(
            marker.read_text(encoding="utf-8")
        ) == {"schema": SCHEMA, "run_id": run_id}
    except (OSError, ValueError):
        return False


def cleanup_resources(directory, record):
    paths, errors = resource_paths(directory)
    removed = []
    for path in paths:
        if not path.exists() and not path.is_symlink():
            continue
        if not owned_path(path, record):
            errors.append(f"Ownership could not be verified: {path}")
            continue
        try:
            shutil.rmtree(path)
            removed.append(str(path))
        except OSError as error:
            errors.append(f"{path}: {error}")
    return {
        "status": "INCOMPLETE" if errors else "CLEANED",
        "removed": removed,
        "errors": errors,
    }


def child_environment(record, directory, temp):
    env = os.environ.copy()
    env.update(UMA_RUN_ID=record["run_id"], UMA_RUN_DIR=str(directory))
    env.update({key: str(temp) for key in ("TMPDIR", "TEMP", "TMP")})
    # imagej.initialize_imagej sets java.io.tmpdir before JVM startup.
    # Preserve heap options, numerical settings and persistent caches.
    return env


def event(directory, message):
    with (directory / "runtime.log").open("a", encoding="utf-8") as handle:
        handle.write(f"{utc_now()} | {message}\n")


def block_cleanup(reason):
    """Keep evidence if shutdown or result publication failed."""
    if current_run() is not None:
        _CLEANUP_BLOCKS.append(str(reason))


def record_completion(code, normal, outcomes=None):
    """Record completion after closing journals and ImageJ workers."""
    run = current_run()
    if run is None:
        return
    run_id, directory = run
    write_json(
        directory / "completion.json",
        {
            "schema": SCHEMA,
            "run_id": run_id,
            "exit_code": code,
            "normal_completion": bool(normal and not _CLEANUP_BLOCKS),
            "process": identity(),
            "outcomes": outcomes or {},
            "cleanup_blocks": list(_CLEANUP_BLOCKS),
            "finished_utc": utc_now(),
        },
    )


def read_completion(directory, record, code, worker_pid):
    """Distinguish normal PARTIAL completion from a crashed worker."""
    try:
        path = directory / "completion.json"
        if path.is_symlink() or path.stat().st_size > 2_000_000:
            raise ValueError("Linked or oversized worker completion record")
        value = json.loads(path.read_text(encoding="utf-8"))
        process = value["process"]
        if (
            value["schema"] != SCHEMA
            or value["run_id"] != record["run_id"]
            or value["exit_code"] != code
            or process["pid"] != worker_pid
            or process_state(process, record["host"]) != "EXITED"
        ):
            raise ValueError(
                "Worker completion identity or exit does not match"
            )
        return value
    except (OSError, ValueError, KeyError, TypeError) as error:
        return {"normal_completion": False, "reason": str(error)}


def managed(step):
    """Run the public console command in its own analysis process."""
    command = next(name for name, value in COMMANDS.items() if value == step)

    def decorate(function):
        @functools.wraps(function)
        def invoke(argv=None):
            from .cli import parse_arguments

            arguments = list(sys.argv[1:] if argv is None else argv)
            parse_arguments(step, arguments)
            return run_command(command, arguments)

        return invoke

    return decorate


def worker():
    """Set the temp directory before importing processing libraries."""
    temp = Path(os.environ["TMPDIR"])
    if Path(tempfile.gettempdir()) != temp:
        raise RuntimeError("Python did not accept the UMA temporary directory")
    command = sys.argv[1]
    if command not in COMMANDS:
        raise ValueError(f"Unknown UMA command: {command}")
    run = current_run()
    if run is None:
        raise RuntimeError("Missing UMA runtime registration")
    _, directory = run
    write_json(
        directory / "worker.json",
        {
            "process": identity(),
            "started_utc": utc_now(),
        },
    )
    from . import cli

    sys.argv = [command, *sys.argv[2:]]
    return getattr(cli, COMMANDS[command]).__wrapped__()


def _observe(child, record):
    """
    Small process-tree sample. Summed RSS includes shared pages more
    than once.
    """
    if record.get("observation_unavailable"):
        return
    try:
        processes = [psutil.Process(child.pid)]
        processes += processes[0].children(recursive=True)
    except psutil.NoSuchProcess:
        return
    except psutil.AccessDenied:
        record["process_observation_incomplete"] = True
        return
    observed = {
        (item["pid"], item.get("created")): item
        for item in record["processes"]
    }
    rss = 0
    for process in processes:
        try:
            item = identity(process)
            observed[(item["pid"], item.get("created"))] = item
            rss += process.memory_info().rss
        except psutil.NoSuchProcess:
            continue
        except psutil.AccessDenied:
            record["process_observation_incomplete"] = True
    record["processes"] = list(observed.values())
    record["sampled_peak_tree_rss_bytes"] = max(
        record.get("sampled_peak_tree_rss_bytes", 0), rss
    )
    try:
        record["sampled_peak_system_swap_bytes"] = max(
            record.get("sampled_peak_system_swap_bytes", 0),
            psutil.swap_memory().used,
        )
    except (OSError, psutil.Error):
        record["swap_observation_unavailable"] = True


def run_command(command, argv, worker_path=None):
    """
    Run one analysis command; cleanup only after its whole observed tree
    exits.
    """
    if command not in COMMANDS:
        raise ValueError(f"Unknown UMA command: {command}")
    try:
        home = runtime_home()
    except (OSError, ValueError) as error:
        print(
            f"Cannot prepare UMA temporary storage: {error}", file=sys.stderr
        )
        return 1
    run_id = (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ_")
        + uuid4().hex[:12]
    )
    directory = home / "runs" / run_id
    temp = home / "tmp" / run_id
    try:
        directory.mkdir(parents=True, mode=0o700)
        temp.parent.mkdir(parents=True, exist_ok=True)
        observable = observation_supported()
        record = {
            "schema": SCHEMA,
            "run_id": run_id,
            "command": command,
            "arguments": list(argv),
            "cwd": str(Path.cwd()),
            "host": platform.node(),
            "started_utc": utc_now(),
            "supervisor": identity()
            if observable
            else {"pid": os.getpid(), "created": None},
            "processes": [],
            "process_observation_incomplete": not observable,
            "observation_unavailable": not observable,
            "process_status": "STARTING",
            "cleanup": {"status": "PENDING"},
        }
        create_owned_directory(temp, run_id, directory)
        write_json(directory / "run.json", record)
        event(
            directory,
            f"STARTED | command={command} | temporary_directory={temp}",
        )
    except (OSError, ValueError) as error:
        print(
            f"Cannot prepare UMA temporary storage: {error}", file=sys.stderr
        )
        return 1
    print(f"Run resources: {directory}", flush=True)
    if not observable:
        print(
            "Process inspection is unavailable; "
            "temporary resources will be retained for review.",
            flush=True,
        )
    worker = (
        [str(worker_path)]
        if worker_path is not None
        else [
            "-c",
            "from uma_tools.runtime import worker; raise SystemExit(worker())",
        ]
    )
    child = None
    terminated = False
    journal_failed = False

    def terminate(signum, frame):
        nonlocal terminated
        terminated = True
        if child is not None and child.poll() is None:
            child.send_signal(signum)

    previous_term = signal.signal(signal.SIGTERM, terminate)
    try:
        child = subprocess.Popen(
            [sys.executable, *worker, command, *argv],
            env=child_environment(record, directory, temp),
        )
        record["process_status"] = "RUNNING"
        while True:
            try:
                _observe(child, record)
                try:
                    write_json(directory / "run.json", record)
                except OSError as error:
                    if not journal_failed:
                        print(
                            "Runtime journal unavailable; "
                            f"temporary files will be retained: {error}",
                            file=sys.stderr,
                        )
                    journal_failed = True
                code = child.wait(timeout=1)
                break
            except subprocess.TimeoutExpired:
                continue
            except KeyboardInterrupt:
                # Both foreground processes receive terminal Ctrl+C.
                # Keep supervising the worker's own cancellation.
                terminated = True
                continue
    except OSError as error:
        print(
            f"UMA runtime error: {error}. Resources retained: {directory}",
            file=sys.stderr,
        )
        # A worker may still be alive after a registry write failed.
        # Keep its files.
        record.update(process_status="SUPERVISION_ERROR", error=str(error))
        try:
            write_json(directory / "run.json", record)
        except OSError:
            pass
        return 1
    finally:
        signal.signal(signal.SIGTERM, previous_term)

    states = [
        process_state(item, record["host"]) for item in record["processes"]
    ]
    record.update(
        process_status="EXITED",
        exit_code=code,
        finished_utc=utc_now(),
        interrupted=terminated,
    )
    completion = read_completion(directory, record, code, child.pid)
    record["completion"] = completion
    # A recorded normal PARTIAL exit is safe; an unrecorded exit is not.
    if (
        completion.get("normal_completion")
        and not journal_failed
        and all(state == "EXITED" for state in states)
        and not record["process_observation_incomplete"]
    ):
        try:
            record["cleanup"] = cleanup_resources(directory, record)
        except OSError as error:
            record["cleanup"] = {
                "status": "INCOMPLETE",
                "errors": [str(error)],
            }
    else:
        record["cleanup"] = {
            "status": "RETAINED",
            "reason": (
                "No verified normal completion, live descendant "
                "or uncertain process state"
            ),
        }
    try:
        write_json(directory / "run.json", record)
        event(
            directory,
            f"PROCESS_EXITED | exit_code={code} | "
            f"cleanup={record['cleanup']['status']}",
        )
    except OSError as error:
        print(f"Could not save runtime completion: {error}", file=sys.stderr)
        return code or 1
    print(
        f"Temporary files: {record['cleanup']['status']}. "
        f"Runtime log: {directory / 'runtime.log'}",
        flush=True,
    )
    return (128 - code if code < 0 else code) or (1 if journal_failed else 0)
