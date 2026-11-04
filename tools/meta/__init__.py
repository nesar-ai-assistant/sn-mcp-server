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

from ..common import REPO_ROOT, TIERS, ArtifactResult, ok
from .skills import skill_index, skill_text

_CFG = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["mcp-server"]
PREFIX = _CFG["skill_prefix"]


def _datasets() -> dict:
    from ..registry import DATASETS
    return DATASETS


def _list_datasets(
    tier: Annotated[int | None, Field(ge=1, le=3, description="Only datasets with products of this data tier "
                                      "(1 raw/catalog, 2 summary statistics, 3 published inference). null = all.")] = None,
) -> ArtifactResult:
    """List the public datasets this server exposes, organised by DATA TIER. Call this first.

    Tier 1 = raw / catalog-level data, tier 2 = processed summary statistics
    with covariances, tier 3 = published posteriors / parameter constraints.
    Analyses should be built on tier-1 and tier-2 products; tier 3 is for
    comparison only. Each entry gives products per tier, access (bundled vs
    downloaded on first use), and which tools read it.
    """
    ds = _datasets()
    brief = {}
    for k, v in ds.items():
        products = v.get("products", {})
        if tier is not None and not products.get(f"tier{tier}"):
            continue
        brief[k] = {kk: v.get(kk) for kk in ("title", "products", "access", "tools")}
    counts = {t: sum(1 for v in ds.values() if v.get("products", {}).get(f"tier{t}")) for t in (1, 2, 3)}
    return ok(f"{len(brief)} datasets{'' if tier is None else f' with tier-{tier} products'}. "
              f"Datasets with tier-1/2/3 products: {counts[1]}/{counts[2]}/{counts[3]}. "
              "Build analyses on tier 1-2; tier 3 is comparison-only. "
              f"Use describe_{PREFIX}_dataset for details.",
              datasets=brief, tiers={str(k): v for k, v in TIERS.items()}, tier_counts=counts)


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

for _name, _fn in list(_EXPORTS.items()):
    globals()[_name].data_tier = "meta"

__all__ = list(_EXPORTS)  # noqa: PLE0604 - names depend on skill_prefix
