---
name: hpc-dispatch-handoff
description: Run this server's heavy kernel (sn.kernels.fit_sn_cosmology - w0waCDM / fine wCDM grids on the full SN covariance) on ALCF Polaris or NERSC Perlmutter through the hep-genesis IRI engine - client-side pack export (hosted) or server-side set_dispatch (own machine)
---

# HPC dispatch handoff

This server holds NO facility credentials, project, workdir or environment
settings. All facility access goes through IRI (job submission) and Globus /
the IRI filesystem API (staging) via the hep-genesis engine; nothing here
writes batch scripts or opens SSH sessions. Tool names `set_dispatch`,
`get_dispatch`, `auth_status`, `export_dispatch_pack` are a cross-server
convention matched by the hep-genesis harness.

## Kernels in this server
| kernel | when worth dispatching | duration_hint_s |
|---|---|---|
| `sn.kernels.fit_sn_cosmology` | `model="w0wacdm"` with n_grid >= 61, or many wCDM fits (cut variations, all samples) | 1800 |
| `sn.kernels.standardize_sn_salt2` | Large datasets with extended SALT2 fits (tier-1 → tier-2 Tripp standardisation) | 120 |

Args (JSON): `{"model": "lcdm"|"wcdm"|"w0wacdm", "source": {"dataset":
"pantheonplus"|"des5yr"|"union3"|"des_dovekie", "zmin": 0.01, "zmax": 10,
"calibrators": "include"}, "ranges": {"w0": [-2, 0], "wa": [-6, 3]} | null,
"n_grid": 0, "refine": true}`. The kernel downloads the public release files
itself (GitHub raw URLs in `kernels.SOURCES`); no local paths are passed.
Alternatively pass `"data": {"z", "zhel", "mu", "cov"}` arrays. pip_deps:
numpy, scipy, astropy. Returns best_fit, chi2_min, 1D/2D gridded marginals
with HPD levels - identical to the local tool's payload.

Local timings for scale (Apple M-series, Pantheon+ 1590 SNe): lcdm 401 pts
2 s; wcdm 121^2 9 s; w0wacdm 41^3 30-50 s; w0wacdm 81^3 ~5 min (scales as n^3).

## Client-side (hosted servers - recommended)
1. `export_dispatch_pack` -> the client saves `files` locally (do NOT echo
   their contents) and reads `kernels`: each has `function`, `args_shape`,
   `pip_deps`, `duration_hint_s`, `returns`.
2. On the client's facility server call
   `run_pack_kernel(pack=<saved>/tools, function="sn.kernels.fit_sn_cosmology",
   args={...}, pip_deps=["numpy","scipy","astropy"], duration>=1800,
   artifact_dir=<run dir>)`.
3. Long or parallel jobs (e.g. one job per dataset): `wait=False`, then
   `find_dispatched_jobs` / `collect_dispatched_job(record=...)`.

## Server-side (server runs on the user's own machine)
1. `auth_status` -> IRI + transfer tokens must be "available".
2. `set_dispatch("perlmutter" | "polaris")`. `fit_sn_cosmology` (and
   `compare_sn_samples mode="contours"`, which calls it per sample) now run each
   fit as one facility job; the output CSVs/JSON are still written locally and
   `computed_on` reports the facility host. Light tools stay local.
   `set_dispatch("local")` to undo.
3. Never blind-retry a failed dispatch: surface the error. A walltime
   rejection may be resubmitted ONCE with a shorter duration.
