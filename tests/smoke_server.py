"""Live-server smoke test: start the stdio server and list tools/prompts.

    python tests/smoke_server.py
"""

import sys

import anyio
from mcp import StdioServerParameters
from mcp.client.session import ClientSession
from mcp.client.stdio import stdio_client


async def main() -> int:
    params = StdioServerParameters(command=sys.executable, args=["-m", "mcp_server", "--transport", "stdio"])
    async with stdio_client(params) as (r, w), ClientSession(r, w) as s:
        await s.initialize()
        tools = (await s.list_tools()).tools
        prompts = (await s.list_prompts()).prompts
        print(f"{len(tools)} tools:", ", ".join(t.name for t in tools))
        print(f"{len(prompts)} prompts:", ", ".join(p.name for p in prompts))
        return 0 if tools else 1


if __name__ == "__main__":
    raise SystemExit(anyio.run(main))
