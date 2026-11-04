"""Optional HPC dispatch: run this server's compute kernels on DOE facilities.

Local execution is the default and is untouched by this module. Two modes
(identical to the cosmic-emulator / spectra / gaia servers):

CLIENT-SIDE (hosted deployments - recommended): this server holds NO
credentials, project, workdir or facility environment. export_dispatch_pack
hands the tools/ kernels to the client, and the client's hep-genesis harness
- which owns the facility tokens, Globus endpoint, project and workdir -
stages and submits them itself over IRI + Globus (facility-server
run_pack_kernel tool / sidecar POST /jobs pack=). A VM deployment needs
nothing beyond this server.

SERVER-SIDE (running the server on your own machine): set_dispatch("polaris"
| "perlmutter") routes the heavy tools through the hep-genesis dispatch
engine in THIS process: the tools/ package is staged to the facility over
Globus (or the IRI filesystem API on NERSC), the kernel runs on a compute
node, and results come back. Facility settings (project, workdir, tokens)
come from the USER's hep-genesis .env and sign-in - never from this server.
Do not use server-side mode on shared/hosted deployments.

Kernels are pip-kernels: each entry in a tool family's DISPATCH_KERNELS
manifest names a pure function in tools/ (JSON-safe args -> JSON-safe dict,
files written to the job cwd) plus the pip requirements the node needs. If
tools/env/conda-lock.yml exists it is shipped and built on the node instead
(hep-genesis dispatch contract R8). Kernel modules must import only numpy /
scipy / astropy / the pip_deps they declare - never pydantic, matplotlib, mcp.

No SLURM/PBS scripts and no SSH anywhere: the engine generates the job
wrapper and submits it through the IRI compute API.

Tool names (set_dispatch, get_dispatch, auth_status, export_dispatch_pack)
are a cross-server convention matched by name in the hep-genesis harness -
do not rename them.
"""

import importlib
import os
import tomllib
from contextvars import ContextVar
from pathlib import Path

TOOLS_DIR = Path(__file__).resolve().parents[1] / "tools"
# The environment the pack carries (relative to tools/); the engine ships and
# builds it by default (hep_genesis.iri.dispatch.envsetup.ENV_LOCK_RELPATH).
ENV_LOCK = "env/conda-lock.yml"
ENV_FILES = ("env/conda-lock.yml", "env/environment.yml")

_state = {"site": "local"}

# The output_dir of the tool call being served (set by the server's tool
# wrapper; asyncio.to_thread carries it into the worker). A dispatched job's
# files and its dispatch_<job>.json record go to <output_dir>/<job>/ - the
# CALL's own run folder. The process-wide DISPATCH_ARTIFACT_DIR is shared by
# every conversation on this server: one chat's set_dispatch(artifact_dir=)
# pulled other chats' jobs into its run folder (2026-10-14).
CALL_OUTPUT_DIR: ContextVar[str | None] = ContextVar("call_output_dir", default=None)


def _job_placement(function: str, args: dict) -> dict:
    """job_name + artifact_dir for this call's job ({} = engine defaults)."""
    out = CALL_OUTPUT_DIR.get()
    if not out:
        return {}
    try:
        from hep_genesis.iri.dispatch.staging import default_job_name
    except ImportError:  # older engine: keep its defaults
        return {}
    job_name = default_job_name(function)
    return {"job_name": job_name, "artifact_dir": str(Path(out).expanduser() / job_name)}

__all__ = ["auth_status", "export_dispatch_pack", "get_dispatch", "set_dispatch"]


def _server_name() -> str:
    for directory in Path(__file__).resolve().parents:
        path = directory / "pyproject.toml"
        if path.exists():
            return tomllib.loads(path.read_text(encoding="utf-8"))["project"]["name"]
    return "mcp-server"

# Mirror the hep-genesis engine's pack limits so an exported pack is always
# stageable by the client that receives it.
_PACK_MAX_FILES = 200
_PACK_MAX_BYTES = 10 * 1024 * 1024


def _collect_kernel_manifest() -> dict:
    """Gather DISPATCH_KERNELS from every tool family that declares one."""
    from .server import configured_tool_module_names

    kernels: dict = {}
    for name in configured_tool_module_names():
        if name == __name__:
            continue
        try:
            module = importlib.import_module(name)
        except Exception:  # noqa: BLE001 - manifest stays useful without a family
            continue
        kernels.update(getattr(module, "DISPATCH_KERNELS", {}) or {})
    return kernels


def export_dispatch_pack() -> dict:
    """Hand THIS server's compute kernels to the client for HPC runs.

    Returns the ``tools/`` package as in-band text files plus a manifest of
    kernel entry points (dotted function path, args shape, node-side
    pip_deps, walltime hint). The CLIENT - which holds the facility tokens,
    Globus endpoint, project and workdir - stages and submits the kernels
    with ITS OWN credentials (e.g. the hep-genesis facility servers'
    run_pack_kernel tool). Nothing credential-shaped or project-shaped is
    needed on, or returned by, this server.

    Clients that recognise the ``dispatch_pack`` key should save the files
    locally and NOT echo their contents into the conversation.
    """
    files: dict[str, str] = {}
    total = 0
    env_paths = [TOOLS_DIR / rel for rel in ENV_FILES if (TOOLS_DIR / rel).is_file()]
    for path in sorted(TOOLS_DIR.rglob("*.py")) + env_paths:
        rel = path.relative_to(TOOLS_DIR)
        if any(part.startswith((".", "__pycache__")) for part in rel.parts):
            continue
        text = path.read_text()
        total += len(text.encode())
        if len(files) >= _PACK_MAX_FILES or total > _PACK_MAX_BYTES:
            raise RuntimeError(
                f"tools/ package exceeds pack limits ({_PACK_MAX_FILES} files / "
                f"{_PACK_MAX_BYTES} bytes) - prune before exporting."
            )
        files[f"tools/{rel}"] = text

    has_lock = f"tools/{ENV_LOCK}" in files
    return {
        "dispatch_pack": {
            "name": "tools",
            "server": _server_name(),
            "file_count": len(files),
            "bytes": total,
            "files": files,
            "kernels": _collect_kernel_manifest(),
            "env_lock": ENV_LOCK if has_lock else None,
            "usage": (
                "Save these files on the CLIENT machine, then dispatch with the "
                "hep-genesis facility server: run_pack_kernel(pack=<saved dir>/tools, "
                "function=<kernel function>, args={...}, pip_deps=<kernel pip_deps>, "
                "duration=<walltime s >= duration_hint_s>). Kernels write their output "
                "files into the job cwd; they come back as artifact_files. "
                + ("The pack's env/conda-lock.yml is built on the node - pass no pip_deps."
                   if has_lock else "No conda lock: pass each kernel's pip_deps.")
            ),
        }
    }


def _engine():
    """Import the dispatch engine, with install instructions on failure."""
    try:
        from hep_genesis.iri.alcf.dispatch import run_codes_on_polaris
        from hep_genesis.iri.nersc.dispatch import run_codes_on_perlmutter
    except ImportError as exc:
        raise RuntimeError(
            "HPC dispatch needs the hep-genesis backend in this server's "
            "environment. Install it with: pip install -e "
            "<hep-genesis-agent>/backend[iri]  (import failed: " + str(exc) + ")"
        ) from exc
    return {"polaris": run_codes_on_polaris, "perlmutter": run_codes_on_perlmutter}


def remote_site() -> str | None:
    """The active remote site ('polaris'/'perlmutter'), or None for local.

    Not an MCP tool - compute tools call this to branch at call time.
    """
    site = _state["site"]
    return None if site == "local" else site


def run_kernel(function: str, args: dict, pip_deps: list[str] | None = None,
               duration: int = 600, nodes: int = 1, env_lock: str | bool | None = None,
               tasks_per_node: int | None = None, cpus_per_task: int | None = None,
               gpus_per_task: int | None = None, wait: bool = True) -> dict:
    """Run a kernel from this server's tools/ package on the active site.

    ``wait=False`` (hep-genesis contract R10) returns as soon as the
    scheduler accepts the job: the engine's ``{"status": "submitted",
    "job_id", "record", ...}`` handle. :func:`collect_kernel` fetches the
    result later from that record.

    ``nodes`` > 1 is an MPI job (contract R9): the kernel's module must
    declare ``DISPATCH_MPI`` or the engine refuses before staging;
    ``tasks_per_node`` / ``cpus_per_task`` / ``gpus_per_task`` override that
    declaration. The engine dict then carries an ``mpi`` block (size, nodes)
    and ``result`` is rank 0's return.

    Not an MCP tool. ``function`` is a dotted path within tools/ (e.g.
    "bao.kernels.fit_bao"); ``args`` must
    be JSON-safe. Returns the engine dict: result, host, and artifact_files
    (local paths of files the kernel wrote in its job directory, fetched
    back automatically where the staging path supports it).

    ``env_lock``: None = the pack's own environment if tools/env/conda-lock.yml
    exists (built on the node), else a venv with ``pip_deps``; False = never
    build the lock.

    Raises RuntimeError with remediation text on any dispatch failure -
    callers should surface it, and agents must NOT blind-retry.
    """
    site = _state["site"]
    if site == "local":
        raise RuntimeError("run_kernel called with local dispatch - use the kernel directly.")
    run = _engine()[site]
    try:
        shape = {k: v for k, v in (("tasks_per_node", tasks_per_node), ("cpus_per_task", cpus_per_task),
                                   ("gpus_per_task", gpus_per_task)) if v}
        result = run(function=function, args=args, codes=str(TOOLS_DIR),
                     pip_deps=pip_deps, duration=duration, nodes=nodes, env_lock=env_lock,
                     **_job_placement(function, args), **shape, **({} if wait else {"wait": False}))
    except Exception as exc:
        body = ""
        resp = getattr(exc, "response", None)
        if resp is not None:
            try:
                body = (resp.text or "")[:600]
            except Exception:
                pass
        from hep_genesis.iri.dispatch.hints import dispatch_auth_hint
        hint = dispatch_auth_hint(site, body or str(exc))
        walltime = any(k in (body + str(exc)).lower()
                       for k in ("qosmaxwalldurationperjoblimit", "walltime", "qos cap", "queue cap"))
        advice = (
            "This is a WALLTIME rejection, not auth: resubmit ONCE with a duration the queue admits "
            "(shorter walltime_s / fewer samples, or a longer queue) - never the same arguments."
            if walltime else
            "Do NOT resubmit the same arguments - the failure is in the dispatch layer and will recur. "
            "Surface this error to the user."
        )
        raise RuntimeError(
            f"Remote dispatch to {site} failed: {exc}\n"
            f"{('Response body: ' + body) if body else ''}{hint}\n{advice}"
        ) from exc
    if not wait and result.get("status") == "submitted":
        return result
    return _check_result(result, site)


def _check_result(result: dict, site: str) -> dict:
    """Raise the agent-facing error for a failed engine/kernel result."""
    if result.get("status") != "success":
        raise RuntimeError(
            f"Remote job on {site} failed: {result.get('error', 'unknown error')}\n"
            f"{('Stderr: ' + result['stderr']) if result.get('stderr') else ''}\n"
            "Do NOT retry this call - surface this error to the user."
        )
    inner = result.get("result")
    if isinstance(inner, str) and inner.startswith("Error"):
        raise RuntimeError(f"Kernel on {site} reported: {inner}\nDo NOT retry blindly - "
                           "fix the inputs or the env_setup and call again.")
    return result


def collect_kernel(record: str) -> dict:
    """One non-blocking look at a job submitted with ``run_kernel(wait=False)``
    (or one whose waiting call timed out), from its ``dispatch_<job>.json``
    record path. Returns the engine dict like :func:`run_kernel` once the job
    has finished, or ``{"status": "pending" | "unknown", "state", "message"}``
    while it has not; raises RuntimeError when it ended without results.

    Not an MCP tool. Works whatever the active site is: the record names it.
    """
    try:
        from hep_genesis.iri.dispatch.collect import collect_dispatched, load_record
    except ImportError as exc:
        raise RuntimeError("Collecting a no-wait job needs a hep-genesis backend with "
                           "hep_genesis.iri.dispatch.collect (update the engine install).") from exc
    rec = load_record(record)
    result = collect_dispatched(rec)
    if isinstance(result, dict) and result.get("status") in ("pending", "unknown"):
        return result
    return _check_result(result, rec["site"])


def set_dispatch(site: str, artifact_dir: str | None = None) -> str:
    """Set where THIS server's heavy tools execute: 'local' (this machine),
    'polaris' (ALCF) or 'perlmutter' (NERSC).

    Affects only this server's heavy (dispatchable) tools - never other
    servers' dispatch state. Remote sites run each heavy call as one facility
    job (staging + queue + walltime: minutes, not seconds) and need facility
    sign-in - check with auth_status first. Light tools (listing, plotting,
    small reads) always run locally. Optional artifact_dir sets where
    produced files are fetched back to. Returns the active configuration.
    """
    if artifact_dir:
        os.environ["DISPATCH_ARTIFACT_DIR"] = str(Path(artifact_dir).expanduser())
    site = (site or "").strip().lower()
    if site in ("local", "off", "none", ""):
        _state["site"] = "local"
        return "Dispatch: local execution."
    if site in ("polaris", "alcf"):
        site = "polaris"
    elif site in ("perlmutter", "nersc"):
        site = "perlmutter"
    else:
        return f"Unknown site {site!r}. Use 'local', 'polaris', or 'perlmutter'."
    _engine()  # fail now, with install instructions, if the engine is missing
    os.environ["DISPATCH_TARGET"] = site
    os.environ.setdefault(
        "DISPATCH_ARTIFACT_DIR", str(Path.home() / "hep-genesis" / "artifacts")
    )
    _state["site"] = site
    facility = "ALCF" if site == "polaris" else "NERSC"
    caps = ("Perlmutter walltime caps: debug 30 min, regular 12 h"
            if site == "perlmutter" else
            "Polaris walltime caps: debug 60 min, preemptable 72 h")
    return (
        f"Dispatch: remote on {site} ({facility}) - heavy tools will stage this "
        "server's kernels, submit via IRI, and fetch results back. Each call is "
        f"one facility job (minutes to hours; progress heartbeats keep the call alive). {caps}."
    )


def get_dispatch() -> str:
    """Report where THIS server's heavy tools currently execute."""
    site = _state["site"]
    if site == "local":
        return "Dispatch: local execution."
    return f"Dispatch: remote on {site}."


def auth_status() -> str:
    """Report the facility sign-in state THIS server sees for its own
    HPC dispatch (ALCF and NERSC).

    Call this before set_dispatch to a remote site, or when one of this
    server's remote calls fails with an auth error. It says nothing about
    the facility servers' own tokens - ask those servers directly. Sign-in
    happens outside this server (hep-genesis auth CLIs or the desktop app's
    HPC panel); tokens are re-read on every call.
    """
    try:
        from hep_genesis.iri.alcf.client import (
            _live_iri_token as alcf_iri,
        )
        from hep_genesis.iri.alcf.client import (
            _live_transfer_token as alcf_tx,
        )
        from hep_genesis.iri.nersc.client import (
            _live_iri_token as nersc_iri,
        )
        from hep_genesis.iri.nersc.client import (
            _live_transfer_token as nersc_tx,
        )
    except ImportError:
        return (
            "hep-genesis backend not installed in this server's environment - "
            "server-side remote dispatch unavailable (the client-side "
            "export_dispatch_pack handoff still works). Install: pip install -e "
            "<hep-genesis-agent>/backend[iri]"
        )
    lines = []
    for name, iri, tx in (("ALCF", alcf_iri, alcf_tx), ("NERSC", nersc_iri, nersc_tx)):
        lines.append(
            f"{name}: IRI token {'available' if iri() else 'MISSING'}, "
            f"transfer token {'available' if tx() else 'MISSING'}"
        )
    lines.append(f"{get_dispatch()} Sign in via the hep-genesis auth CLIs or the app's HPC panel.")
    return "\n".join(lines)
