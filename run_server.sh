#!/bin/bash
# stdio launcher for agent configs (Hermes, Claude Code ...): runs from the repo root
cd "$(dirname "$0")"
exec .venv/bin/python -m mcp_server "$@"
