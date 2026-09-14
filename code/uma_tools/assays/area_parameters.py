"""Raw float32 threshold configuration shared by the Area workflow."""

from __future__ import annotations

import math
import struct

FLOAT32_MAX = 3.4028234663852886e38
DEFAULT_THRESHOLD_LOWER = 2000.0


class ValidationError(Exception):
    """
    An input or selection problem that must not be silently bypassed.
    """


def float32_limit(value, label):
    """
    Report the effective precision used by ImageJ's 32-bit thresholding.
    """
    try:
        number = float(value)
        if not math.isfinite(number) or number < 0 or number > FLOAT32_MAX:
            raise ValueError("outside the finite non-negative float32 range")
        effective = struct.unpack("!f", struct.pack("!f", number))[0]
    except (ValueError, TypeError, OverflowError, struct.error) as error:
        raise ValidationError(
            f"{label} must be finite, non-negative, "
            "and representable in 32-bit float."
        ) from error
    return number, effective


def threshold_settings(lower, upper):
    if lower is None:
        if upper is not None:
            raise ValidationError(
                "An upper threshold without a lower threshold is not valid."
            )
        return None
    requested_lower, effective_lower = float32_limit(lower, "Lower threshold")
    unbounded = upper is None or str(upper).strip().lower() in (
        "inf",
        "+inf",
        "infinity",
        "none",
    )
    if unbounded:
        requested_upper, effective_upper = None, FLOAT32_MAX
    else:
        requested_upper, effective_upper = float32_limit(
            upper, "Upper threshold"
        )
        if requested_upper < requested_lower:
            raise ValidationError(
                "Upper threshold must be greater than or equal "
                "to the lower threshold."
            )
    return {
        "requested_lower": requested_lower,
        "requested_upper": requested_upper,
        "lower": effective_lower,
        "upper": effective_upper,
        "upper_unbounded": unbounded,
    }
