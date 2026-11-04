# cosmo-*-mcp-server family conventions

All servers in this family are built from `~/github-integration/cosmo-mcp-template`
and must stay structurally identical, so an agent that knows one knows them all.
The wrapper is derived from `desc_mcp_server` (HEP-KE pattern, with IRI dispatch).

## Layout
```
<name>-mcp-server/
  pyproject.toml          [tool.mcp-server] name, skill_prefix, instructions_file, tool_modules
  .mcp.json  run_server.sh  .python-version  .gitignore  LICENSE  README.md
  .github/workflows/ci.yml
  mcp_server/             GENERIC - copied verbatim, never edited per server
    server.py cli.py dispatch.py __init__.py __main__.py README.md
    INSTRUCTIONS.md       per-server instructions string (the only per-server file here)
  tools/
    common.py             ArtifactResult, ok(), fetch_file(), data_cache_dir(), resolve_outdir(), mpl(), save_fig()
    registry.py           DATASETS dict - single source of truth (title, products, access, urls, citation, tools)
    meta/                 list_<p>_datasets, describe_<p>_dataset, list_<p>_skills, load_<p>_skill
    <family>/__init__.py  one sub-package per tool family; __all__ = exposed tools; DISPATCH_KERNELS manifest
    <family>/io.py        readers / downloaders (pure python)
    <family>/kernels.py   dispatchable compute (JSON in -> JSON out, files to cwd); imports ONLY numpy/scipy/astropy/declared pip_deps
  skills/*.md             server-side skills (frontmatter name/description); also served as MCP prompts
                          REQUIRED: <p>-tour.md, hpc-dispatch-handoff.md, + >=1 science recipe
  data/                   small bundled tables only (< ~5 MB total); everything else downloads on first use
  tests/test_contract.py  family contract (copied verbatim)
  tests/test_<family>.py  science tests; network tests marked @pytest.mark.network
  tests/smoke_server.py   stdio server smoke test
  docs/CONVENTIONS.md     this file (copied verbatim; per-server notes go in docs/SERVER_NOTES.md)
```

## Data tiers (family-wide; the point of these servers is tier-1/2 analysis)
| tier | what | examples | use in hep-genesis analyses |
|---|---|---|---|
| 1 | raw / catalog-level | galaxy & cluster catalogs, SN light curves / SALT2 fits, shear catalogs, CMB maps | PRIMARY input |
| 2 | processed summary statistics + covariance | xi+-, C_ell, CMB bandpowers, BAO D/r_d, Hubble diagram mu(z), binned counts, n(z), selection functions | PRIMARY input |
| 3 | published inference products | posterior chains, S8/H0/Om constraints, best-fit theory | COMPARISON ONLY |

Rules:
- Every science tool carries `@tier(1|2|3|"meta")` (from tools.common, above `@validate_call`). The server
  prefixes the tool description with the tier; results carry `metadata["data_tier"]`; tier-3 results also
  carry `metadata["tier3_notice"]`.
- Never mix tiers inside one tool: a tool that fits tier-2 data must not silently add a tier-3 posterior as
  a prior. Gaussian priors are explicit, user-supplied arguments, labelled as such in the output.
- `tools/registry.py` products are split `{"tier1": [...], "tier2": [...], "tier3": [...]}`.
- Every server must expose >= 1 tier-1 and >= 1 tier-2 dataset with tools (contract-tested).
- Tier-1 data that is large is downloaded on demand in SUBSETS (sky patch, tile, z range, nside) and heavy
  tier-1 -> tier-2 reductions (e.g. measuring xi(s), xi+- or C_ell from a catalog/map) are DISPATCH_KERNELS.
- Skills build analyses from tier 1-2 and end with a "compare to published (tier 3)" step.

## Tool rules
- Every tool: `@validate_call`, `Annotated[..., Field(description=...)]` args, docstring that says
  WHEN to use it, returns `ArtifactResult` (via `ok(message, files, **metadata)`).
- Files, not arrays: fetch/compute tools write CSV (or npz/FITS) to `output_dir`; plot tools take
  those paths and write PNGs. Keep metadata small (summary numbers, shapes, column names).
- Naming: `fetch_*` (download + write tidy CSV), `compute_*`/`fit_*` (numbers), `plot_*` (PNG),
  `compare_*` (cross-dataset). Dataset chosen by a `dataset` arg keyed into registry.DATASETS.
- Data access: `tools.common.fetch_file(url, subdir=...)` - bundled data/ first, then cache
  (`$COSMO_MCP_DATA` or `~/.cache/cosmo-mcp/<project-name>`), then download. Never download at import.
- Units/conventions documented in INSTRUCTIONS.md and in each describe entry.
- Observational comparisons cite sources (Author+Year, journal / arXiv).

## Weights + HPC dispatch (IRI)
- `fn.weight = "heavy"` for anything slow (likelihood grids, MCMC-lite, big downloads + reductions);
  heavy tools are serialized and heartbeated by the server.
- Heavy tools with a kernel branch like gaia/desc:
  ```python
  from mcp_server.dispatch import remote_site, run_kernel   # lazy, inside the tool
  site = remote_site()
  if site:
      res = run_kernel("<family>.kernels.<fn>", {...json args...}, pip_deps=[...], duration=900)
      payload = res["result"]; files = res.get("artifact_files", [])
  else:
      payload = kernels.<fn>(**args)
  ```
- Each family declares
  `DISPATCH_KERNELS = {"<key>": {"function": "<family>.kernels.<fn>", "args_shape": "...",
  "pip_deps": [...], "duration_hint_s": N, "returns": "..."}}`.
- Tool names set_dispatch / get_dispatch / auth_status / export_dispatch_pack come from
  mcp_server.dispatch (listed last in tool_modules) - never rename.

## Ports (streamable-http defaults; set `port = N` in [tool.mcp-server], read by cli.py)
sn 8010, galaxy-surveys 8011, cmb-data 8012, weak-lensing 8013, clusters 8014.
