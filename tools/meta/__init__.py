"""Discovery tools every server in the family carries.

Names are prefixed with the server's ``skill_prefix`` from pyproject
([tool.mcp-server]) so several servers can be mounted in one client without
tool-name collisions: list_<p>_datasets, describe_<p>_dataset,
list_<p>_skills, load_<p>_skill.
"""

from __future__ import annotations

import tomllib
from typing import Annotated

from pydantic import Field, validate_call

from ..common import REPO_ROOT, ArtifactResult, ok
from .skills import skill_index, skill_text

_CFG = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["mcp-server"]
PREFIX = _CFG["skill_prefix"]


def _datasets() -> dict:
    from ..registry import DATASETS
    return DATASETS


def _list_datasets() -> ArtifactResult:
    """List the public datasets this server exposes. Call this first.

    Each entry: survey/release, what products are available (summary
    statistics, covariances, chains, catalogs), where the data comes from
    (bundled in the repo vs downloaded on first use, with the source URL),
    the citation, and which tools use it.
    """
    ds = _datasets()
    brief = {k: {kk: v.get(kk) for kk in ("title", "products", "access", "tools")} for k, v in ds.items()}
    return ok(f"{len(ds)} datasets. Use describe_{PREFIX}_dataset for details.", datasets=brief)


def _describe_dataset(
    name: Annotated[str, Field(min_length=1, description="Dataset key from the list tool.")],
) -> ArtifactResult:
    """Full description of one dataset: products, units, conventions, URLs, citation."""
    ds = _datasets()
    if name not in ds:
        raise ValueError(f"Unknown dataset {name!r}. Available: {', '.join(ds)}")
    return ok(f"Dataset {name}.", dataset=ds[name])


def _list_skills() -> ArtifactResult:
    """List this server's skills: named multi-tool recipes.

    A skill is procedural know-how - which tools to chain, in which order,
    with what parameters, and how to sanity-check and report. When a task
    matches one, load it with the load-skill tool and follow it.
    """
    index = skill_index()
    return ok(f"{len(index)} skills available. Load one when a task matches its description.",
              skills=index)


def _load_skill(
    name: Annotated[str, Field(min_length=1, description="Skill name from the list-skills tool.")],
) -> ArtifactResult:
    """Load the full instructions of one of this server's skills and follow them."""
    text = skill_text(name)
    if text is None:
        raise ValueError(f"No skill named {name!r}. Available: {', '.join(sorted(skill_index()))}")
    return ok(f"Loaded skill {name!r}. Follow these instructions.", instructions=text)


_EXPORTS = {
    f"list_{PREFIX}_datasets": _list_datasets,
    f"describe_{PREFIX}_dataset": _describe_dataset,
    f"list_{PREFIX}_skills": _list_skills,
    f"load_{PREFIX}_skill": _load_skill,
}
for _name, _fn in _EXPORTS.items():
    _fn = validate_call(_fn)
    _fn.__name__ = _fn.__qualname__ = _name
    globals()[_name] = _fn

__all__ = list(_EXPORTS)  # noqa: PLE0604 - names depend on skill_prefix
