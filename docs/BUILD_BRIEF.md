# Build brief for cosmo-*-mcp-server builders (read fully before starting)

You are building ONE standalone MCP server in a family of five, built in parallel by
independent agents. They must be structurally identical. Do NOT touch other servers'
directories, and do NOT edit any Hermes config (~/.hermes/...) - servers are NOT to be enabled.

## Starting point (mandatory)
1. `cp -R ~/github-integration/cosmo-mcp-template ~/github-integration/<repo-name>`, then
   `rm -rf .git && git init`.
2. Read `docs/CONVENTIONS.md` (the family rules) and `mcp_server/README.md`.
   `mcp_server/` is GENERIC: do NOT edit server.py / cli.py / dispatch.py (except
   INSTRUCTIONS.md). `tests/test_contract.py` and `docs/CONVENTIONS.md` stay verbatim.
3. Reference implementations to read (do not modify):
   - `~/github-integration/desc_mcp_server/`: tools/augur_tools/__init__.py (DISPATCH_KERNELS +
     remote_site/run_kernel branch), tools/meta, skills/*.md (skill style), tests/
   - `~/github-integration/gaia-mcp-server-dispatch-update/tools/cmd_tools.py`: fetch tool with
     dispatch branch (lines ~55-120), plot tool style
4. Edit pyproject.toml: `[project].name`, description, dependencies (add only what you need:
   astropy, pandas, h5py, getdist, sacc, ... - prefer pure-Python wheels; heavy/optional
   deps go to optional-dependencies and are lazily imported inside functions);
   `[tool.mcp-server]` name, skill_prefix, tool_modules (tools.meta first,
   mcp_server.dispatch last). Edit .mcp.json server key, README.md, INSTRUCTIONS.md.
5. Delete `tools/example/` and `skills/template-tour.md`; fill in `tools/registry.py` DATASETS.

## Data policy
- Most data DOWNLOADS ON DEMAND through `tools.common.fetch_file(url, subdir=...)` (cache in
  `~/.cache/cosmo-mcp/<project-name>`, override with `$COSMO_MCP_DATA`). Never download at import.
- Bundle in `data/` ONLY small canonical tables (total < ~5 MB) so the core demo works offline,
  e.g. BAO mean+cov, a binned Hubble diagram. Record provenance (URL + citation) for each
  bundled file in `data/README.md`.
- VERIFY EVERY URL by actually downloading it during development. If a URL is dead, find the
  real one (GitHub raw, Zenodo, LAMBDA, survey portal) and note it. Never invent URLs or numbers.
  Survey notes with candidate URLs are in `~/github-integration/cosmo_mcp_survey_report.md`
  (also `~/cmb_mcp_data_survey.md`, `~/bao_rsd_data_survey.md`,
  `~/weak_lensing_survey_data_products.md`, `~/cosmological_data_survey.md`).

## Tools (the gaia pattern: data access + preliminary plots)
- 8-15 science tools. Per dataset family: `fetch_*` (download + tidy CSV in output_dir),
  `compute_*`/`fit_*` (simple, fast, honest analysis, e.g. a chi2 grid fit of LCDM/wCDM
  using numpy/scipy distances - astropy.cosmology is fine), `plot_*` (publication-quality PNG
  that reproduces a recognisable figure from the release paper), and `compare_*` across surveys.
- Every tool follows the CONVENTIONS rules (validate_call, Annotated Field descriptions,
  docstring says when to use it, returns ArtifactResult via ok(), file paths not arrays).
- At least ONE heavy tool (e.g. a likelihood grid / parameter fit) with `fn.weight = "heavy"`,
  a `<family>/kernels.py` pure function, a DISPATCH_KERNELS manifest entry, and the
  `remote_site()/run_kernel()` branch exactly as in CONVENTIONS.md. Kernels import only
  numpy/scipy/astropy + their declared pip_deps, take JSON args (pass data arrays or a
  download URL, NOT local paths), and return JSON.

## Server-side skills (skills/*.md, frontmatter `name:` + `description:`)
- Required: `<prefix>-tour.md` (one call per family, all-local, < 5 min) and the template's
  `hpc-dispatch-handoff.md` (adapt the kernel names), plus >= 2 science recipes
  (e.g. "reproduce Fig. N of <paper>", "compare surveys and quantify tension",
  "fit w0wa"). Each recipe gives exact tool calls, parameter values, expected numbers
  as sanity checks, and what to report. Style: see desc_mcp_server/skills/*.md.

## Environment + verification (must actually run)
- `cd <repo> && uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -e ".[dev]"`
  (uv venvs have no pip; on a macOS cryptography build failure, `uv pip install "cryptography<44"` first).
- `MPLBACKEND=Agg`. Must pass: `.venv/bin/ruff check tools mcp_server tests`,
  `.venv/bin/python -m pytest -q` (including network tests - they must pass here at least once),
  and `.venv/bin/python tests/smoke_server.py`.
- Then run the tour end to end in plain Python (call the tool functions) into
  `<repo>/output/tour/`, and look at the numbers for sanity (e.g. the fitted Om is near 0.3).
  Keep output/ out of git, but save 2-3 representative PNGs to `docs/figures/` for the README.
- `git add -A && git commit` locally. Do NOT push, do NOT create GitHub repos.

## Final report (your last message)
repo path; tool list (name - one line); datasets with access mode (bundled/download) and
verified URLs; skills list; heavy/dispatch kernels; test results (paste the pytest + smoke
summary lines); the tour's key numbers; figure paths; known gaps/TODOs. Be honest: report
anything that failed or that you could not verify.
