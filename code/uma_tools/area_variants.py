"""Stable Area variant identities and collection snapshot paths."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from .config import ValidationError

COLLECTION_SCHEMA = "uma-area-variants-v1"
FLOAT32_MAX = float.fromhex("0x1.fffffep+127")


def digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, allow_nan=False).encode("utf-8")
    ).hexdigest()


def inside(root, relative):
    """Resolve a manifest-relative regular path without following links."""
    root = Path(root)
    path = Path(relative)
    if path.is_absolute() or not path.parts or ".." in path.parts:
        raise ValidationError(f"Unsafe collection path: {relative!r}")
    current = root
    for part in path.parts:
        current /= part
        if current.is_symlink():
            raise ValidationError(f"Symbolic collection path: {current}")
    return current


def profile(rows):
    """Read the recorded effective limits, never infer them from a name."""
    profiles = set()
    for row in rows:
        try:
            lower = float(row["Threshold_Lower"])
            upper = float(
                row.get("Effective_Threshold_Upper")
                or row.get("Threshold_Upper")
                or FLOAT32_MAX
            )
            projection = row["Projection_Method"].strip().upper()
            channel_text = row.get("Channel_Index", "").strip()
            channel = int(channel_text) if channel_text else None
        except (KeyError, ValueError, TypeError) as error:
            raise ValidationError(
                "Area CSV must record threshold limits and projection"
            ) from error
        if (
            not math.isfinite(lower)
            or not math.isfinite(upper)
            or upper < lower
            or projection not in {"SUM", "MEAN", "MAX"}
            or (channel is not None and channel < 1)
        ):
            raise ValidationError("Invalid Area threshold/projection/channel")
        profiles.add((lower, upper, projection, channel))
    if len(profiles) != 1:
        raise ValidationError("One Area run contains mixed threshold settings")
    lower, upper, projection, channel = profiles.pop()
    return {
        "lower": lower,
        "upper": upper,
        "projection": projection,
        "channel": channel,
    }


def variant_id(settings):
    return digest(settings)[:16]


def variant_label(settings):
    lower = format(settings["lower"], ".9g")
    upper = (
        "max"
        if settings["upper"] == FLOAT32_MAX
        else format(settings["upper"], ".9g")
    )
    channel = settings["channel"] or "unknown"
    return f"AreaThr_{lower}_to_{upper}_{settings['projection']}_C{channel}"


def sort_key(settings):
    return (
        settings["lower"],
        settings["upper"],
        settings["projection"],
        settings["channel"] or 0,
    )
