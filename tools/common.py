"""Shared infrastructure for every tool family in this server.

Conventions (identical across the cosmo-*-mcp-server family):
- every tool returns ArtifactResult (status, files, message, metadata)
- arrays move between tools as CSV/FITS/NPZ FILE PATHS, never raw arrays
- tools that write files take ``output_dir``
- remote data is downloaded on first use into the data cache
  (``$COSMO_MCP_DATA`` or ``~/.cache/cosmo-mcp/<server>``); small canonical
  tables may be bundled in ``data/`` and are always tried first
- heavy tools set ``fn.weight = "heavy"`` (serialized, heartbeated, HPC
  dispatchable); default weight is "light"
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import os
import shutil
import threading
import urllib.request
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel

REPO_ROOT = Path(__file__).resolve().parents[1]
BUNDLED_DATA = REPO_ROOT / "data"


class ArtifactResult(BaseModel):
    """Uniform result contract returned by every tool."""

    status: Literal["success"]
    files: list[str]
    message: str
    metadata: dict[str, Any]


def ok(message: str, files: list | None = None, **metadata) -> ArtifactResult:
    return ArtifactResult(status="success", files=[str(f) for f in (files or [])],
                          message=message, metadata=metadata)


# ---------------------------------------------------------------- caching
_CACHE: dict[str, Any] = {}
_CACHE_LOCK = threading.Lock()


def get_cached(key: str, factory: Callable[[], Any]) -> Any:
    """Load-once in-process cache for parsed datasets."""
    with _CACHE_LOCK:
        if key not in _CACHE:
            with contextlib.redirect_stdout(io.StringIO()):
                _CACHE[key] = factory()
        return _CACHE[key]


def server_slug() -> str:
    import tomllib
    cfg = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    return cfg["project"]["name"]


def data_cache_dir() -> Path:
    root = os.environ.get("COSMO_MCP_DATA")
    path = Path(root).expanduser() / server_slug() if root else Path.home() / ".cache" / "cosmo-mcp" / server_slug()
    path.mkdir(parents=True, exist_ok=True)
    return path


def fetch_file(url: str, filename: str | None = None, *, subdir: str = "",
               sha256: str | None = None, timeout: float = 120.0) -> Path:
    """Return a local path for ``url``: bundled data/ copy, else cache, else download.

    Downloads are atomic (tmp file + rename) so an interrupted fetch never
    leaves a truncated file in the cache.
    """
    name = filename or url.rstrip("/").split("/")[-1].split("?")[0]
    bundled = BUNDLED_DATA / subdir / name
    if bundled.is_file():
        return bundled
    dest_dir = data_cache_dir() / subdir
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / name
    if dest.is_file() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "cosmo-mcp-server/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as fh:
            shutil.copyfileobj(resp, fh)
    except Exception as exc:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"Download failed for {url}: {exc}. Check network access, "
                           f"or place the file manually at {dest}.") from exc
    if sha256:
        digest = hashlib.sha256(tmp.read_bytes()).hexdigest()
        if digest != sha256:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(f"Checksum mismatch for {url}")
    tmp.replace(dest)
    return dest


# ---------------------------------------------------------------- outputs
def resolve_outdir(output_dir: str) -> Path:
    outdir = Path(output_dir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def param_slug(params: dict) -> str:
    """Short stable hash so scans never overwrite each other's files."""
    blob = ",".join(f"{k}={params[k]}" for k in sorted(params))
    return hashlib.sha1(blob.encode()).hexdigest()[:6]


def mpl():
    """Headless matplotlib pyplot (Agg)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def save_fig(fig, outdir: Path, stem: str) -> Path:
    path = outdir / f"{stem}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return path
