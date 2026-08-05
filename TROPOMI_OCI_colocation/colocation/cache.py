"""Output directories and NPZ+JSON stage caches with fingerprints."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np

STAGE_DIRS = ("pairs", "matches", "svd", "rsr", "products")
STAGES = ("pairs", "matches", "svd", "rsr", "simulate")


def ensure_output_dirs(output_dir: Path) -> dict[str, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {name: output_dir / name for name in STAGE_DIRS}
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    return paths


def fingerprint(payload: Any) -> str:
    """Stable SHA256 hex digest of a JSON-serializable payload."""
    blob = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def write_json(path: Path, obj: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str) + "\n")


def read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text())


def save_npz_json(
    npz_path: Path,
    arrays: dict[str, np.ndarray],
    meta_path: Path,
    meta: dict[str, Any],
) -> None:
    npz_path = Path(npz_path)
    meta_path = Path(meta_path)
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(npz_path, **arrays)
    write_json(meta_path, meta)


def load_npz(path: Path) -> dict[str, np.ndarray]:
    with np.load(Path(path), allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def is_fresh(
    meta_path: Path,
    expected_fp: str,
    *,
    required_files: Iterable[Path] | None = None,
    enabled: bool = True,
) -> bool:
    if not enabled:
        return False
    meta_path = Path(meta_path)
    if not meta_path.is_file():
        return False
    if required_files is not None:
        for p in required_files:
            if not Path(p).is_file():
                return False
    try:
        meta = read_json(meta_path)
    except (OSError, json.JSONDecodeError):
        return False
    return meta.get("fingerprint") == expected_fp


def write_run_meta(output_dir: Path, cfg_dict: dict[str, Any], extra: dict[str, Any] | None = None) -> None:
    payload = {
        "written_at": datetime.now(timezone.utc).isoformat(),
        "config": cfg_dict,
        "config_fingerprint": fingerprint(cfg_dict),
    }
    if extra:
        payload.update(extra)
    write_json(Path(output_dir) / "run_meta.json", payload)


def force_set(force: list[str] | None) -> set[str]:
    if not force:
        return set()
    out: set[str] = set()
    for item in force:
        item = item.strip().lower()
        if item == "all":
            return set(STAGES)
        if item not in STAGES:
            raise ValueError(f"Unknown force stage {item!r}; expected one of {STAGES} or 'all'")
        out.add(item)
    return out


def stage_order_index(stage: str) -> int:
    order = {s: i for i, s in enumerate(STAGES)}
    if stage not in order:
        raise ValueError(f"Unknown stage {stage!r}")
    return order[stage]
