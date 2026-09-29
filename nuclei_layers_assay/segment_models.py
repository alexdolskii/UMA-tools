"""Find and validate explicitly selected local StarDist 3D models."""

from __future__ import annotations

import math
from pathlib import Path

from uma_tools.config import load_json

MODEL_ROOT = Path(__file__).resolve().parent / "stardist_models_nuclei_layers"
MODEL_FILES = ("config.json", "thresholds.json", "weights_best.h5")


def model_directory(path: Path) -> Path:
    """Support the one bundled model with an extra 'model' directory."""
    if (
        not (path / "config.json").is_file()
        and (path / "model" / "config.json").is_file()
    ):
        return path / "model"
    return path


def bundled_models() -> dict[str, Path]:
    """List usable assets without importing TensorFlow or StarDist."""
    if not MODEL_ROOT.is_dir():
        return {}
    return {
        path.name: model_directory(path)
        for path in sorted(MODEL_ROOT.iterdir())
        if path.is_dir()
        and not path.name.startswith((".", "_"))
        and all(
            (model_directory(path) / name).is_file() for name in MODEL_FILES
        )
    }


def resolve_model(value: str, config_path: Path) -> Path:
    """Resolve names, moved repository paths and explicit custom paths.

    Existing CWD-relative paths retain the old script's meaning. A second
    different match beside the JSON is rejected rather than guessed.
    """
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Choose a model with model_path or --model")
    path = Path(value).expanduser()
    if path.is_absolute():
        if not path.is_dir():
            raise ValueError(
                f"Model folder does not exist: {path}. Models moved to "
                "nuclei_layers_assay/stardist_models_nuclei_layers; "
                "use --list-models and select its short name."
            )
        return model_directory(path).resolve()
    choices = bundled_models()
    if value in choices:
        return choices[value].resolve()
    candidates = []
    for parent in (Path.cwd(), config_path.resolve().parent):
        candidate = model_directory(parent / path)
        if candidate.is_dir():
            candidates.append(candidate.resolve())
    # Both the old and new repository prefixes identify the moved assets.
    parts = path.parts
    if not candidates and "stardist_models_nuclei_layers" in parts:
        index = parts.index("stardist_models_nuclei_layers")
        suffix = parts[index + 1 :]
        if suffix and ".." not in suffix:
            candidate = model_directory(MODEL_ROOT.joinpath(*suffix))
            if candidate.is_dir():
                candidates.append(candidate.resolve())
    candidates = list(dict.fromkeys(candidates))
    if len(candidates) > 1:
        raise ValueError(
            "Ambiguous model path; use a short bundled name or absolute path: "
            + ", ".join(map(str, candidates))
        )
    if not candidates:
        raise ValueError(
            f"Model not found: {value}. Run uma_nla_segment --list-models."
        )
    return candidates[0]


def validate_model(path: Path) -> dict:
    """Require trained single-channel 3D weights and explicit thresholds."""
    for name in MODEL_FILES:
        asset = path / name
        if not asset.is_file() or asset.stat().st_size == 0:
            raise ValueError(f"Missing nonempty model file: {asset}")
    config = load_json(path / "config.json")
    if (
        config.get("n_dim") != 3
        or config.get("n_channel_in") != 1
        or config.get("axes") not in ("ZYXC", "ZYX")
    ):
        raise ValueError(
            f"Expected a single-channel StarDist 3D model: {path}"
        )
    thresholds = load_json(path / "thresholds.json")
    for key in ("prob", "nms"):
        number = thresholds.get(key)
        if (
            type(number) not in (float, int)
            or not math.isfinite(number)
            or not 0 < number < 1
        ):
            raise ValueError(f"Invalid model threshold {key}: {number}")
    return {key: float(thresholds[key]) for key in ("prob", "nms")}


def choose_model(value: str | None, config_path: Path, interactive: bool):
    """Ask once if unspecified; never substitute a demonstration model."""
    if value is None or value == "":
        choices = bundled_models()
        if not interactive or not choices:
            raise ValueError(
                "Specify model_path in the JSON or --model. "
                "Run uma_nla_segment --list-models to see available names."
            )
        print("Select the trained model appropriate for your images:")
        names = list(choices)
        for index, name in enumerate(names, 1):
            print(f"  {index}: {name}")
        answer = input("Model number or name: ").strip()
        if answer.isdigit() and 1 <= int(answer) <= len(names):
            value = names[int(answer) - 1]
        elif answer in choices:
            value = answer
        else:
            raise ValueError("Invalid model selection")
    path = resolve_model(value, config_path)
    return path, validate_model(path)
