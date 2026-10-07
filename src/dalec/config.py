"""Read config/default.yaml and resolve paths against the repository root.

Every script gets its settings (paths, years, thresholds) from here, so no
script ever types a path or a year directly.
"""

from pathlib import Path
from typing import Any

import yaml

# src/dalec/config.py -> src/dalec -> src -> repository root
REPO_ROOT: Path = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH: Path = REPO_ROOT / "config" / "default.yaml"

def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load the YAML config. Defaults to config/default.yaml."""
    config_path = Path(path) if path is not None else DEFAULT_CONFIG_PATH
    if not config_path.is_file():
        raise FileNotFoundError(f"config file not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError(f"config did not parse to a mapping: {config_path}")
    return config

def resolve_path(value: str | Path) -> Path:
    """Turn a config path into an absolute path (relative = from repo root)."""
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path

def year_block(config: dict[str, Any], key: str) -> tuple[int, int]:
    """Return an inclusive (start, end) year pair from config['years'][key]."""
    block = (config.get("years") or {}).get(key)
    if block is None or len(block) != 2:
        raise ValueError(f"years.{key} must be an inclusive [start, end] pair, got {block!r}")
    start, end = int(block[0]), int(block[1])
    if end < start:
        raise ValueError(f"years.{key} ends before it starts: {block!r}")
    return start, end