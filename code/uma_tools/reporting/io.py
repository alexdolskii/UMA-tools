"""Report output policies over shared file and log utilities."""

from __future__ import annotations

from functools import partial

from ..common.files import save_csv as _save_csv
from ..common.files import save_json as _save_json
from ..common.files import sha256_file
from ..common.run import RunLog
from ..common.run import utc_now as _utc_now

save_csv = partial(_save_csv, encoding="utf-8-sig")
save_json = partial(
    _save_json, atomic=False, allow_nan=False, trailing_newline=False
)
utc_now = partial(_utc_now, timespec="seconds")


def save_details(path, rows):
    """
    Write every diagnostic field without dropping affected records.
    """
    columns = list(dict.fromkeys(key for row in rows for key in row))
    save_csv(path, columns or ["Issue"], rows)


__all__ = [
    "RunLog",
    "save_csv",
    "save_json",
    "sha256_file",
    "utc_now",
    "save_details",
]
