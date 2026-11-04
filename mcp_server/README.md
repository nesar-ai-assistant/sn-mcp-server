# mcp_server/ - generic wrapper (do not edit per server)

Identical across the cosmo-*-mcp-server family (derived from
desc_mcp_server). Per-server settings live in pyproject.toml:

```toml
[tool.mcp-server]
name = "..."                 # FastMCP server name
skill_prefix = "..."         # list_<p>_skills / load_<p>_skill / list_<p>_datasets
instructions_file = "mcp_server/INSTRUCTIONS.md"
tool_modules = ["tools.meta", "tools.<family>", "mcp_server.dispatch"]
```

- `server.py`: registers every name in each module's `__all__`; tools run
  off the event loop serialized by `fn.weight` ("light" default / "heavy"),
  with MCP progress heartbeats; hosted mode via `MCP_OUTPUT_ROOT`,
  `MCP_ARTIFACT_URL`, `MCP_PUBLIC=1`; skills/*.md also served as prompts.
- `dispatch.py`: optional HPC dispatch through the hep-genesis IRI engine
  (set_dispatch / get_dispatch / auth_status / export_dispatch_pack).
- `cli.py`: `python -m mcp_server --transport stdio|streamable-http --host --port`.
