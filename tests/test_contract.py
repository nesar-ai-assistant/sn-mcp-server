"""Family-wide contract tests: every server in the cosmo-*-mcp family passes these."""

import importlib
import json

from mcp_server.dispatch import _PACK_MAX_BYTES, _PACK_MAX_FILES, export_dispatch_pack
from mcp_server.server import configured_tool_module_names, server_config
from tools.common import ArtifactResult
from tools.meta import PREFIX
from tools.meta.skills import skill_index


def test_tool_modules_import_and_define_all():
    for name in configured_tool_module_names():
        mod = importlib.import_module(name)
        assert hasattr(mod, "__all__"), name
        for fn in mod.__all__:
            assert callable(getattr(mod, fn)), (name, fn)


def test_meta_tools_and_skills():
    meta = importlib.import_module("tools.meta")
    r = getattr(meta, f"list_{PREFIX}_datasets")()
    assert isinstance(r, ArtifactResult) and r.metadata["datasets"]
    sk = skill_index()
    assert "hpc-dispatch-handoff" in sk
    assert len(sk) >= 2, "each server ships at least a tour skill + dispatch skill"
    r = getattr(meta, f"load_{PREFIX}_skill")("hpc-dispatch-handoff")
    assert "IRI" in r.metadata["instructions"]


def test_dispatch_tools_registered():
    assert "mcp_server.dispatch" in configured_tool_module_names()
    assert server_config()["skill_prefix"] == PREFIX


def test_pack_shape_and_safety():
    pack = export_dispatch_pack()["dispatch_pack"]
    assert pack["name"] == "tools"
    assert 0 < pack["file_count"] <= _PACK_MAX_FILES
    assert 0 < pack["bytes"] <= _PACK_MAX_BYTES
    json.dumps(pack)
    for rel in pack["files"]:
        assert rel.startswith("tools/") and ".." not in rel and "__pycache__" not in rel
    for key, k in pack["kernels"].items():
        assert {"function", "pip_deps", "duration_hint_s"} <= set(k), key



def _science_kinds():
    out = {}
    for name in configured_tool_module_names():
        if name == "mcp_server.dispatch":
            continue
        mod = importlib.import_module(name)
        for fn in mod.__all__:
            out[f"{name}.{fn}"] = getattr(getattr(mod, fn), "data_kind", None)
    return out


def test_every_tool_declares_its_data_kind():
    for tool, kind in _science_kinds().items():
        assert kind in ("catalog", "summary_statistic", "published_result", "meta"), \
            f"{tool} has no @data_kind(...) (got {kind!r})"


def test_registry_offers_catalogs_and_summary_statistics():
    from tools.registry import DATASETS
    n = {"catalog": 0, "summary_statistic": 0}
    for key, d in DATASETS.items():
        products = d.get("products")
        assert isinstance(products, dict), f"{key}: products must be a dict keyed by data kind"
        assert set(products) <= {"catalog", "summary_statistic", "published_result"}, key
        for k in n:
            n[k] += bool(products.get(k))
    assert n["catalog"] >= 1, "server must expose at least one catalog-level dataset"
    assert n["summary_statistic"] >= 1, "server must expose at least one summary-statistic dataset"


def test_catalog_and_summary_statistic_tools_exist():
    kinds = set(_science_kinds().values())
    assert {"catalog", "summary_statistic"} <= kinds, f"need catalog and summary-statistic tools, have {kinds}"


def test_published_results_are_stamped_comparison_only():
    from tools.common import ArtifactResult, data_kind, ok

    @data_kind("published_result")
    def published():
        return ok("x")

    r = published()
    assert isinstance(r, ArtifactResult)
    assert r.metadata["data_kind"] == "published_result" and "comparison only" in r.metadata["comparison_only"]
