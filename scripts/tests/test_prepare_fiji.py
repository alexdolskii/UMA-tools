"""Dependency integrity, bounded retries and failure propagation."""

import hashlib
import importlib.util
import io
import subprocess
from pathlib import Path
from unittest.mock import Mock
from urllib.error import HTTPError

import pytest

SPEC = importlib.util.spec_from_file_location(
    "prepare_fiji", Path(__file__).parents[1] / "prepare_fiji.py"
)
bootstrap = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bootstrap)


@pytest.fixture
def artifact(monkeypatch):
    data = b"synthetic artifact, never passed to real Maven"
    monkeypatch.setattr(bootstrap, "SHA256", hashlib.sha256(data).hexdigest())
    return data


def test_verified_download_and_offline_cache(tmp_path, monkeypatch, artifact):
    fetch = Mock(return_value=artifact)
    monkeypatch.setattr(bootstrap, "download", fetch)
    path = bootstrap.artifact_file(tmp_path)
    assert path.read_bytes() == artifact
    fetch.assert_called_once()
    fetch.side_effect = AssertionError("A valid cache must not need network")
    assert bootstrap.artifact_file(tmp_path) == path


@pytest.mark.parametrize("cached", [False, True])
def test_wrong_digest_never_reaches_maven(
    tmp_path, monkeypatch, artifact, cached
):
    if cached:
        (tmp_path / bootstrap.ARTIFACT).write_bytes(b"corrupt")
    monkeypatch.setattr(bootstrap, "download", Mock(return_value=b"corrupt"))
    monkeypatch.setattr(bootstrap.shutil, "which", lambda _: "/mvn")
    run = Mock()
    monkeypatch.setattr(bootstrap.subprocess, "run", run)
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        bootstrap.prepare(tmp_path)
    assert len(run.call_args_list) == 1  # Version check only.
    assert run.call_args.args[0] == ["/mvn", "-version"]


def test_official_fallback_preserves_digest_check(
    tmp_path, monkeypatch, artifact
):
    fetch = Mock(
        side_effect=[HTTPError("url", 404, "missing", {}, None), artifact]
    )
    monkeypatch.setattr(bootstrap, "download", fetch)
    assert bootstrap.artifact_file(tmp_path).read_bytes() == artifact
    assert [call.args[0] for call in fetch.call_args_list] == list(
        bootstrap.URLS
    )


@pytest.mark.parametrize("status, calls", [(404, 1), (503, 3)])
def test_retry_is_bounded_and_does_not_repeat_404(monkeypatch, status, calls):
    fetch = Mock(side_effect=HTTPError("url", status, "failure", {}, None))
    monkeypatch.setattr(bootstrap, "urlopen", fetch)
    monkeypatch.setattr(bootstrap.time, "sleep", Mock())
    with pytest.raises(HTTPError):
        bootstrap.download("https://example.invalid/artifact")
    assert fetch.call_count == calls


def test_transient_failure_can_recover(monkeypatch):
    fetch = Mock(
        side_effect=[
            HTTPError("url", 503, "temporary", {}, None),
            io.BytesIO(b"downloaded"),
        ]
    )
    monkeypatch.setattr(bootstrap, "urlopen", fetch)
    monkeypatch.setattr(bootstrap.time, "sleep", Mock())
    assert (
        bootstrap.download("https://example.invalid/artifact") == b"downloaded"
    )
    assert fetch.call_count == 2


def test_maven_failure_stops_preparation(tmp_path, monkeypatch, artifact):
    monkeypatch.setattr(bootstrap.shutil, "which", lambda _: "/mvn")
    monkeypatch.setattr(bootstrap, "download", Mock(return_value=artifact))
    run = Mock(side_effect=[None, subprocess.CalledProcessError(1, "mvn")])
    monkeypatch.setattr(bootstrap.subprocess, "run", run)
    with pytest.raises(subprocess.CalledProcessError):
        bootstrap.prepare(tmp_path)
    command = run.call_args.args[0]
    assert bootstrap.INSTALL_GOAL in command
    assert "-DgroupId=woolz" in command
    assert "-DartifactId=JWlz" in command
    assert "-Dversion=1.4.0" in command
    assert "--strict-checksums" in command


@pytest.mark.parametrize(
    "failure",
    [
        subprocess.CalledProcessError(1, "python"),
        subprocess.TimeoutExpired("python", 900),
    ],
)
def test_child_failure_or_hang_is_not_success(tmp_path, monkeypatch, failure):
    monkeypatch.setattr(bootstrap, "prepare", Mock())
    run = Mock(side_effect=failure)
    monkeypatch.setattr(bootstrap.subprocess, "run", run)
    with pytest.raises(type(failure)):
        bootstrap.main(["--check", "--artifact-cache", str(tmp_path)])
    assert run.call_args.kwargs["check"] is True
    assert run.call_args.kwargs["timeout"] == 900
