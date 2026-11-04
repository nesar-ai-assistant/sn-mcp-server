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
