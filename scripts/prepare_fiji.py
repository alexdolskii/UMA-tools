"""Restore the pinned JWlz artifact and optionally validate UMA's real JVM.

Run from a checkout after installing UMA-tools. No analysis code is changed.
The two official archives provide the same bytes; a pinned digest is required
before Maven sees the JAR. TLS validation remains enabled.
"""

import argparse
import hashlib
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

ARTIFACT = "JWlz-1.4.0.jar"
SHA256 = "0662cfe9600c24cbef591e17f0bfbbbebbae8b0e8d83098735454a999166410f"
URLS = (
    "https://sites.imagej.net/Fiji/jars/JWlz-1.4.0.jar-20140304121539",
    "https://downloads.openmicroscopy.org/bio-formats/"
    "5.0.0-rc2/artifacts/JWlz-1.4.0.jar",
)
INSTALL_GOAL = (
    "org.apache.maven.plugins:maven-install-plugin:3.1.3:install-file"
)
TRANSIENT_HTTP = {408, 429, 500, 502, 503, 504}


def verify(data):
    """Reject changed or truncated artifacts before installation."""
    actual = hashlib.sha256(data).hexdigest()
    if actual != SHA256:
        raise ValueError(f"{ARTIFACT}: SHA-256 mismatch: {actual}")


def download(url):
    """Retry temporary network failures, never certificate or digest errors."""
    for attempt in range(3):
        try:
            with urlopen(url, timeout=30) as response:
                return response.read()
        except HTTPError as error:
            if error.code not in TRANSIENT_HTTP or attempt == 2:
                raise
        except URLError as error:
            if not isinstance(error.reason, (TimeoutError, socket.timeout)):
                raise
            if attempt == 2:
                raise
        except TimeoutError:
            if attempt == 2:
                raise
        time.sleep(attempt + 1)
    raise RuntimeError("Download retry limit reached")


def artifact_file(cache):
    """Reuse only verified bytes; publish a download atomically."""
    path = cache / ARTIFACT
    if path.exists():
        verify(path.read_bytes())
        print(f"Verified cached artifact: {path}", flush=True)
        return path
    errors = []
    for url in URLS:
        try:
            data = download(url)
        except (OSError, URLError) as error:
            errors.append(f"{url}: {error}")
            continue
        verify(data)
        cache.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=cache) as temporary:
            staging = Path(temporary) / ARTIFACT
            staging.write_bytes(data)
            staging.replace(path)
        print(f"Verified download: {url}\nSHA-256: {SHA256}", flush=True)
        return path
    raise RuntimeError("Cannot download JWlz:\n" + "\n".join(errors))


def prepare(cache):
    """Install into Maven's configured repository, shared with jgo."""
    executable = shutil.which("mvn")
    if executable is None:
        raise RuntimeError("Maven is missing; activate the UMA environment.")
    subprocess.run([executable, "-version"], check=True, timeout=30)
    artifact = artifact_file(cache)
    # SciJava no longer supplies this artifact or its POM. Generate the
    # coordinate-only POM, retaining the original group, artifact and version.
    subprocess.run(
        [
            executable,
            "--batch-mode",
            "--strict-checksums",
            INSTALL_GOAL,
            f"-Dfile={artifact}",
            "-DgroupId=woolz",
            "-DartifactId=JWlz",
            "-Dversion=1.4.0",
            "-Dpackaging=jar",
            "-DgeneratePom=true",
        ],
        check=True,
        timeout=180,
    )


def check_child():
    """Exercise a real filter and worker shutdown in a disposable process."""
    from scyjava import jimport

    from uma_tools.constants import FIJI_ENDPOINT
    from uma_tools.imagej import initialize_imagej, shutdown_imagej_workers

    context = image = None
    try:
        context = initialize_imagej()
        system = jimport("java.lang.System")
        print(f"Endpoint: {FIJI_ENDPOINT}", flush=True)
        print(f"Java: {system.getProperty('java.version')}", flush=True)
        print(f"Java home: {system.getProperty('java.home')}", flush=True)
        if str(system.getProperty("java.specification.version")) != "11":
            raise RuntimeError("Use OpenJDK 11 in the tested UMA environment.")
        ij = jimport("ij.IJ")
        image = ij.createImage("UMA CI check", "8-bit black", 32, 32, 1)
        image.getProcessor().setValue(17)
        image.getProcessor().fill()
        ij.run(image, "Gaussian Blur...", "sigma=1")
        if image.getStatistics().mean != 17:
            raise RuntimeError("The Fiji filter check returned wrong pixels.")
    finally:
        try:
            if image is not None:
                image.close()
            if context is not None:
                context.dispose()
        finally:
            shutdown_imagej_workers()
    print("Java/Fiji check passed.", flush=True)


def check():
    """Retry only a reported transient Maven bootstrap failure."""
    command = [sys.executable, str(Path(__file__).resolve()), "--check-child"]
    environment = os.environ.copy()
    for attempt in range(3):
        print(f"Starting Fiji check (attempt {attempt + 1}/3)...", flush=True)
        try:
            result = subprocess.run(
                command,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=True,
                timeout=900,
            )
        except subprocess.CalledProcessError as error:
            output = error.stdout or ""
            print(output, end="", flush=True)
            transient = (
                "Failed to bootstrap the artifact." in output
                and re.search(
                    r"status code: (408|429|5\d\d)\b|"
                    r"connect timed out|read timed out|connection reset",
                    output,
                    re.IGNORECASE,
                )
            )
            if not transient or attempt == 2:
                raise
            # Maven caches failed artifact requests. -U retries missing
            # releases; it does not replace pinned Fiji with a new release.
            environment["MAVEN_ARGS"] = (
                os.environ.get("MAVEN_ARGS", "") + " -U"
            ).strip()
            print("Temporary Maven failure; retrying with -U.", flush=True)
            time.sleep(2 * (attempt + 1))
        except subprocess.TimeoutExpired as error:
            output = error.stdout or b""
            if isinstance(output, bytes):
                output = output.decode("utf-8", errors="replace")
            print(output, end="", flush=True)
            raise
        else:
            print(result.stdout, end="", flush=True)
            print("Fiji process exited successfully.", flush=True)
            return


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true", help="Start Fiji and verify clean exit"
    )
    parser.add_argument(
        "--artifact-cache",
        type=Path,
        default=Path.home() / ".cache" / "uma-tools" / "bootstrap",
        help="Cache for the checksum-verified historical JAR",
    )
    parser.add_argument(
        "--check-child", action="store_true", help=argparse.SUPPRESS
    )
    args = parser.parse_args(argv)
    if args.check_child:
        check_child()
        return
    prepare(args.artifact_cache.expanduser().resolve())
    if args.check:
        check()


if __name__ == "__main__":
    main()
