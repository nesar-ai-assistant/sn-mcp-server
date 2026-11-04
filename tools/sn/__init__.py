"""Type Ia supernova cosmology tools: fetch Hubble diagrams, bin them, plot them,
fit flat LCDM / wCDM / w0waCDM with the full covariance, plot and compare contours.

Units: z dimensionless (zHD = CMB-frame + peculiar-velocity corrected), mu in mag,
H0 in km/s/Mpc. Absolute magnitude M and H0 are degenerate in SN-only data; every
fit marginalises the additive offset analytically, so absolute mu levels (and the
H0 passed to compute_distance_modulus) never affect cosmological constraints.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Annotated, Literal

import numpy as np
from pydantic import Field, validate_call

from ..common import ArtifactResult, mpl, ok, param_slug, resolve_outdir, save_fig, tier
from . import io, kernels

__all__ = [
    "fetch_sn_sample",
    "compute_distance_modulus",
    "bin_hubble_diagram",
    "plot_hubble_diagram",
    "compute_sn_chi2",
    "fit_sn_cosmology",
    "plot_sn_contours",
    "compare_sn_samples",
    "plot_sn_covariance",
    # tier-1 tools
    "fetch_sn_salt2_params",
    "standardize_sn_salt2",
    "fetch_sn_lightcurve",
    "plot_sn_lightcurve",
]

DISPATCH_KERNELS = {
    "fit_sn_cosmology": {
        "function": "sn.kernels.fit_sn_cosmology",
        "args_shape": "{'model': 'lcdm'|'wcdm'|'w0wacdm', 'source': {'dataset': 'pantheonplus'|'des5yr'|'union3'|"
                      "'des_dovekie', 'zmin': float, 'zmax': float, 'calibrators': 'include'|'exclude'|'only'}, "
                      "'ranges': {param: [lo, hi]} | None, 'n_grid': int (0 = default 401/121/41), 'refine': bool}  "
                      "(or 'data': {'z','zhel','mu','cov'} inline arrays instead of 'source'; the kernel "
                      "downloads the public release files itself)",
        "pip_deps": ["numpy", "scipy", "astropy"],
        "duration_hint_s": 1800,
        "returns": "best_fit, chi2_min, n_sn, dof, best_offset_mag, marginals_1d {param: x, post, chi2_profile, "
                   "mean, sd, lo68..hi95}, marginals_2d {'A-B': x, y, post, chi2_profile, hpd_levels_68_95}",
    },
    "standardize_sn_salt2": {
        "function": "sn.kernels.standardize_sn_salt2",
        "args_shape": "{'data': {'z','zhel','mB','mBERR','x1','x1ERR','c','cERR','COV_x1_c','HOST_LOGMASS',...}, "
                      "'Om': float (default 0.3)}",
        "pip_deps": ["numpy", "scipy"],
        "duration_hint_s": 120,
        "returns": "alpha, beta, gamma, sigma_int, mu, mu_err, z, chi2_per_dof, rms_vs_released_mu",
    },
}

Dataset = Literal["pantheonplus", "des5yr", "union3", "des_dovekie"]
Calib = Literal["include", "exclude", "only"]
Model = Literal["lcdm", "wcdm", "w0wacdm"]

_LABELS = {"pantheonplus": "Pantheon+", "des5yr": "DES-SN5YR", "union3": "Union3/UNITY1.5",
           "des_dovekie": "DES-Dovekie"}
_COLORS = {"pantheonplus": "tab:blue", "des5yr": "tab:red", "union3": "tab:green", "des_dovekie": "tab:purple"}
_TEX = {"Om": r"$\Omega_m$", "w": r"$w$", "w0": r"$w_0$", "wa": r"$w_a$"}


def _write_csv(path: Path, header: list[str], rows) -> Path:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    return path


def _zmin(dataset: str, zmin: float | None) -> float:
    from ..registry import DATASETS
    return DATASETS[dataset].get("default_zmin", 0.0) if zmin is None else zmin


# ======================================================================= fetch
@tier(2)
@validate_call
def fetch_sn_sample(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the CSV + covariance .npy.")],
    dataset: Annotated[Dataset, Field(description="pantheonplus | des5yr | union3 | des_dovekie (see list_sn_datasets).")] = "pantheonplus",
    zmin: Annotated[float | None, Field(ge=0.0, le=3.0, description="Keep z_HD > zmin. None = release standard "
                                         "(Pantheon+ 0.01, others 0).")] = None,
    zmax: Annotated[float, Field(gt=0.0, le=10.0, description="Keep z_HD <= zmax.")] = 10.0,
    calibrators: Annotated[Calib, Field(description="Pantheon+ Cepheid-host calibrators (IS_CALIBRATOR=1): "
                                        "include | exclude | only. Other samples have none.")] = "include",
) -> ArtifactResult:
    """Download (first use) and tidy one SN Ia Hubble diagram with its FULL stat+sys covariance.

    Use this first: it writes `<dataset>_hd.csv` (cid, survey, z, zhel, mu, mu_err_diag,
    is_calibrator) and `<dataset>_cov.npy` (the matching sub-covariance in mag^2, same row
    order). Every other tool can take the dataset name directly, but the CSV is what you
    hand to plotting or external codes. Union3 is a 22-node binned (compressed) Hubble
    diagram - there is no public per-SN Union3 table.
    """
    zlo = _zmin(dataset, zmin)
    s = io.load_sample(dataset, zlo, zmax, calibrators)
    if len(s["z"]) == 0:
        raise ValueError(f"No {dataset} rows pass zmin={zlo}, zmax={zmax}, calibrators={calibrators}.")
    out = resolve_outdir(output_dir)
    tag = f"{dataset}" if (zmin is None and zmax == 10.0 and calibrators == "include") else \
        f"{dataset}_{param_slug({'zmin': zlo, 'zmax': zmax, 'cal': calibrators})}"
    csv_path = _write_csv(out / f"{tag}_hd.csv",
                          ["cid", "survey", "z", "zhel", "mu", "mu_err_diag", "is_calibrator"],
                          zip(s["cid"], s["survey"], s["z"], s["zhel"], np.round(s["mu"], 6),
                              np.round(s["mu_err"], 6), s["is_calibrator"], strict=True))
    cov_path = out / f"{tag}_cov.npy"
    np.save(cov_path, s["cov"])
    from ..registry import DATASETS
    return ok(f"{len(s['z'])} {_LABELS[dataset]} entries with z in ({zlo:g}, {zmax:g}], calibrators={calibrators}; "
              f"covariance {s['cov'].shape[0]}x{s['cov'].shape[0]}.",
              [csv_path, cov_path], dataset=dataset, n=len(s["z"]), zmin=zlo, zmax=zmax,
              z_range=[float(s["z"].min()), float(s["z"].max())],
              n_calibrators=int(s["is_calibrator"].sum()), mu_kind=s["mu_kind"],
              citation=DATASETS[dataset]["citation"])


# ======================================================================= theory
@tier(2)
@validate_call
def compute_distance_modulus(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the CSV.")],
    Om: Annotated[float, Field(gt=0.0, lt=1.0, description="Matter density Omega_m.")] = 0.3,
    w0: Annotated[float, Field(ge=-3.0, le=1.0, description="Dark-energy w0 (CPL).")] = -1.0,
    wa: Annotated[float, Field(ge=-10.0, le=5.0, description="Dark-energy wa (CPL, w = w0 + wa(1-a)).")] = 0.0,
    H0: Annotated[float, Field(gt=20.0, lt=150.0, description="H0 [km/s/Mpc]; only shifts mu by a constant.")] = 70.0,
    Ok: Annotated[float, Field(ge=-0.5, le=0.5, description="Curvature Omega_k (0 = flat).")] = 0.0,
    z_min: Annotated[float, Field(gt=0.0, le=5.0, description="Lowest redshift.")] = 0.001,
    z_max: Annotated[float, Field(gt=0.0, le=10.0, description="Highest redshift.")] = 2.5,
    n_points: Annotated[int, Field(ge=2, le=20000, description="Number of log-spaced redshifts.")] = 300,
) -> ArtifactResult:
    """Theory distance modulus mu(z) = 5 log10(D_L/10pc) for a given cosmology (w0waCDM, curvature).

    Use for theory curves to overlay on a Hubble diagram, or to quantify how much two
    cosmologies differ in mu (the CSV includes mu minus flat LCDM Om=0.3, offset-free at
    low z). Radiation is neglected (< 1e-4 mag for z < 3). Cross-checked against
    astropy.cosmology in the tests.
    """
    if z_max <= z_min:
        raise ValueError("z_max must exceed z_min")
    z = np.geomspace(z_min, z_max, n_points)
    mu = kernels.distance_modulus(z, om=Om, w0=w0, wa=wa, h0=H0, ok=Ok)
    ref = kernels.distance_modulus(z, om=0.3, h0=H0)
    dl = 10 ** ((mu - 25.0) / 5.0)
    out = resolve_outdir(output_dir)
    pars = {"Om": Om, "w0": w0, "wa": wa, "H0": H0, "Ok": Ok}
    path = _write_csv(out / f"mu_theory_{param_slug(pars)}.csv", ["z", "mu", "dL_Mpc", "mu_minus_lcdm_om0p3"],
                      zip(np.round(z, 6), np.round(mu, 6), np.round(dl, 4), np.round(mu - ref, 6), strict=True))
    return ok(f"mu(z) for {pars} at {n_points} redshifts in [{z_min:g}, {z_max:g}].", [path], params=pars,
              mu_at_z1=float(kernels.distance_modulus([1.0], om=Om, w0=w0, wa=wa, h0=H0, ok=Ok)[0]),
              max_abs_diff_vs_lcdm_om0p3=float(np.max(abs(mu - ref))))


# ================================================================= binning
def _residuals(s: dict, Om: float, w0: float, wa: float):
    """mu - mu_fid with the ML offset (full covariance) removed; returns (resid, offset, icov)."""
    th = kernels.distance_modulus(s["z"], s["zhel"], om=Om, w0=w0, wa=wa)
    icov = np.linalg.inv(s["cov"])
    r = s["mu"] - th
    off = kernels.best_offset(r, icov)
    return r - off, off, icov


def _bin(s: dict, resid: np.ndarray, icov: np.ndarray, edges: np.ndarray):
    """Full-covariance GLS compression of residuals into z bins.

    binned = (A^T C^-1 A)^-1 A^T C^-1 r,  Cov_b = (A^T C^-1 A)^-1, A = bin-membership matrix,
    so correlated systematics (calibration) propagate into the bin covariance.
    """
    idx = np.digitize(s["z"], edges) - 1
    keep = [b for b in range(len(edges) - 1) if np.any(idx == b)]
    A = np.zeros((len(s["z"]), len(keep)))
    for j, b in enumerate(keep):
        A[idx == b, j] = 1.0
    F = A.T @ icov @ A
    covb = np.linalg.inv(F)
    rb = covb @ (A.T @ icov @ resid)
    zb = np.array([np.average(s["z"][idx == b]) for b in keep])
    nb = np.array([int(np.sum(idx == b)) for b in keep])
    return zb, rb, covb, nb, [(edges[b], edges[b + 1]) for b in keep]


@tier(2)
@validate_call
def bin_hubble_diagram(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the binned CSV + covariance.")],
    dataset: Annotated[Dataset, Field(description="SN sample.")] = "pantheonplus",
    n_bins: Annotated[int, Field(ge=2, le=60, description="Number of log-spaced z bins.")] = 20,
    zmin: Annotated[float | None, Field(ge=0.0, le=3.0, description="z cut (None = release standard).")] = None,
    zmax: Annotated[float, Field(gt=0.0, le=10.0, description="Upper z cut.")] = 10.0,
    calibrators: Annotated[Calib, Field(description="include | exclude | only (Pantheon+).")] = "exclude",
    fid_Om: Annotated[float, Field(gt=0.0, lt=1.0, description="Fiducial flat-LCDM Om for residuals.")] = 0.3,
    fid_w0: Annotated[float, Field(ge=-3.0, le=1.0, description="Fiducial w0.")] = -1.0,
    fid_wa: Annotated[float, Field(ge=-10.0, le=5.0, description="Fiducial wa.")] = 0.0,
) -> ArtifactResult:
    """Compress a Hubble diagram into redshift bins of mu residuals vs a fiducial cosmology.

    Use to see trends (e.g. the low-z vs high-z offset that drives Om/w differences between
    samples) and to make compact, correlated data vectors. Residuals are mu - mu_fid with the
    best-fit offset removed; binning is a full-covariance GLS compression, so the bin
    covariance (written as .npy) keeps systematic correlations. Union3 is already binned:
    bins with one node are passed through.
    """
    zlo = _zmin(dataset, zmin)
    s = io.load_sample(dataset, zlo, zmax, calibrators)
    resid, off, icov = _residuals(s, fid_Om, fid_w0, fid_wa)
    edges = np.geomspace(max(s["z"].min(), 1e-4) * 0.999, s["z"].max() * 1.001, n_bins + 1)
    zb, rb, covb, nb, ed = _bin(s, resid, icov, edges)
    eb = np.sqrt(np.diag(covb))
    mu_fid = kernels.distance_modulus(zb, om=fid_Om, w0=fid_w0, wa=fid_wa)
    out = resolve_outdir(output_dir)
    tag = f"{dataset}_binned_{param_slug({'n': n_bins, 'zmin': zlo, 'zmax': zmax, 'cal': calibrators, 'om': fid_Om, 'w0': fid_w0, 'wa': fid_wa})}"
    path = _write_csv(out / f"{tag}.csv", ["z_lo", "z_hi", "z_mean", "n_sn", "dmu", "dmu_err", "mu_fid_plus_dmu"],
                      [(round(a, 5), round(b, 5), round(z, 5), n, round(r, 5), round(e, 5), round(m + r, 5))
                       for (a, b), z, n, r, e, m in zip(ed, zb, nb, rb, eb, mu_fid, strict=True)])
    cpath = out / f"{tag}_cov.npy"
    np.save(cpath, covb)
    chi2_flat = float(rb @ np.linalg.solve(covb, rb))
    return ok(f"{len(zb)} non-empty bins of {_LABELS[dataset]} residuals vs flat w0wa(Om={fid_Om}, w0={fid_w0}, wa={fid_wa}).",
              [path, cpath], dataset=dataset, n_sn=len(s["z"]), n_bins=len(zb), offset_removed_mag=off,
              chi2_binned_vs_fiducial=chi2_flat,
              note="Bins are correlated (systematics); use the covariance, not the error bars, for fits.")


# ================================================================= plotting
@tier(2)
@validate_call
def plot_hubble_diagram(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the PNG.")],
    dataset: Annotated[Dataset, Field(description="SN sample.")] = "pantheonplus",
    zmin: Annotated[float | None, Field(ge=0.0, le=3.0, description="z cut (None = release standard).")] = None,
    zmax: Annotated[float, Field(gt=0.0, le=10.0, description="Upper z cut.")] = 10.0,
    calibrators: Annotated[Calib, Field(description="include | exclude | only (Pantheon+).")] = "include",
    fid_Om: Annotated[float, Field(gt=0.0, lt=1.0, description="Fiducial flat-LCDM Om (curve + residual reference).")] = 0.3,
    fid_w0: Annotated[float, Field(ge=-3.0, le=1.0, description="Fiducial w0.")] = -1.0,
    fid_wa: Annotated[float, Field(ge=-10.0, le=5.0, description="Fiducial wa.")] = 0.0,
    n_bins: Annotated[int, Field(ge=0, le=60, description="Binned residuals overlaid (0 = none).")] = 20,
) -> ArtifactResult:
    """Publication-style Hubble diagram: mu(z) with the fiducial curve and a residual panel.

    Use to show a sample (cf. Brout+22 Fig. 4, DES 2024 Fig. 1): per-SN points with
    diagonal errors (grey), full-covariance binned residuals (colour), the fiducial
    model and Pantheon+ calibrators highlighted. Residuals have the best-fit offset removed
    (M/H0 are degenerate), so only the SHAPE of mu(z) is meaningful.
    """
    zlo = _zmin(dataset, zmin)
    s = io.load_sample(dataset, zlo, zmax, calibrators)
    resid, off, icov = _residuals(s, fid_Om, fid_w0, fid_wa)
    plt = mpl()
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True, gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.04})
    col = _COLORS[dataset]
    cal = s["is_calibrator"].astype(bool)
    zz = np.geomspace(s["z"].min() * 0.9, s["z"].max() * 1.1, 400)
    mu_curve = kernels.distance_modulus(zz, om=fid_Om, w0=fid_w0, wa=fid_wa) + off
    mus = s["mu"]
    alpha = 0.35 if len(s["z"]) > 100 else 0.9
    a1.errorbar(s["z"][~cal], mus[~cal], yerr=s["mu_err"][~cal], fmt=".", ms=2.5, color=col, ecolor="0.75", alpha=alpha,
                lw=0.5, label=f"{_LABELS[dataset]} ({(~cal).sum()})")
    a2.errorbar(s["z"][~cal], resid[~cal], yerr=s["mu_err"][~cal], fmt=".", ms=2, color=col, ecolor="0.8", alpha=alpha, lw=0.5)
    if cal.any():
        a1.errorbar(s["z"][cal], mus[cal], yerr=s["mu_err"][cal], fmt="s", ms=3, color="k", lw=0.6,
                    label=f"Cepheid calibrators ({cal.sum()})")
        a2.errorbar(s["z"][cal], resid[cal], yerr=s["mu_err"][cal], fmt="s", ms=3, color="k", lw=0.6)
    lab = f"flat $\\Lambda$CDM $\\Omega_m$={fid_Om:g}" if (fid_w0 == -1 and fid_wa == 0) else \
        f"$\\Omega_m$={fid_Om:g}, $w_0$={fid_w0:g}, $w_a$={fid_wa:g}"
    a1.plot(zz, mu_curve, "k-", lw=1.2, label=lab)
    files = []
    if n_bins and len(s["z"]) > n_bins:
        edges = np.geomspace(s["z"].min() * 0.999, s["z"].max() * 1.001, n_bins + 1)
        zb, rb, covb, nb, _ = _bin(s, resid, icov, edges)
        a2.errorbar(zb, rb, yerr=np.sqrt(np.diag(covb)), fmt="o", ms=5, color="k", mfc="w", mew=1.2, lw=1.2,
                    label="binned (full cov)", zorder=5)
    # reference: empty (Om=1 EdS-like) is too far; show Om +- 0.1 flat LCDM shape changes (offset-free at z->0)
    for dom, ls in ((-0.1, "--"), (0.1, ":")):
        om_alt = min(max(fid_Om + dom, 0.01), 0.99)
        mu_alt = kernels.distance_modulus(zz, om=om_alt) - kernels.distance_modulus(zz, om=fid_Om, w0=fid_w0, wa=fid_wa)
        a2.plot(zz, mu_alt, "k" + ls, lw=0.8, alpha=0.7, label=f"flat $\\Lambda$CDM $\\Omega_m$={om_alt:.2f}")
    a2.axhline(0, color="k", lw=1)
    a2.legend(loc="lower left", fontsize=7, frameon=False, ncol=3)
    a1.set_xscale("log")
    a1.set_ylabel(r"$\mu$ [mag]")
    a2.set_ylabel(r"$\mu-\mu_{\rm fid}$ [mag]")
    a2.set_xlabel(r"$z_{\rm HD}$")
    a2.set_ylim(-1.0, 1.0)
    a1.legend(loc="upper left", fontsize=8, frameon=False)
    a1.set_title(f"{_LABELS[dataset]} Hubble diagram (N={len(s['z'])})", fontsize=11)
    out = resolve_outdir(output_dir)
    files.append(save_fig(fig, out, f"hubble_{dataset}_{param_slug({'zmin': zlo, 'zmax': zmax, 'cal': calibrators, 'om': fid_Om, 'w0': fid_w0, 'wa': fid_wa, 'nb': n_bins})}"))
    return ok(f"Hubble diagram of {len(s['z'])} {_LABELS[dataset]} entries vs fiducial ({lab}).", files,
              dataset=dataset, n=len(s["z"]), offset_removed_mag=off,
              rms_residual=float(np.std(resid[~cal])) if (~cal).any() else None)


@tier(2)
@validate_call
def plot_sn_covariance(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the PNG.")],
    dataset: Annotated[Dataset, Field(description="SN sample.")] = "pantheonplus",
    zmin: Annotated[float | None, Field(ge=0.0, le=3.0, description="z cut (None = release standard).")] = None,
) -> ArtifactResult:
    """Image of the stat+sys correlation matrix (rows sorted by redshift).

    Use to show how strongly SNe are correlated by calibration/systematics (off-diagonal
    structure is what makes diagonal-error fits wrong) and to sanity-check the matrix.
    """
    s = io.load_sample(dataset, _zmin(dataset, zmin), 10.0, "include")
    o = np.argsort(s["z"])
    c = s["cov"][np.ix_(o, o)]
    d = np.sqrt(np.diag(c))
    corr = c / np.outer(d, d)
    plt = mpl()
    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
    fig.colorbar(im, ax=ax, label="correlation")
    ax.set_title(f"{_LABELS[dataset]} STAT+SYS correlation (sorted by z)", fontsize=10)
    ax.set_xlabel("index")
    ax.set_ylabel("index")
    off = corr[~np.eye(len(corr), dtype=bool)]
    path = save_fig(fig, resolve_outdir(output_dir), f"cov_{dataset}")
    ev = np.linalg.eigvalsh(s["cov"])
    return ok(f"{_LABELS[dataset]} correlation matrix ({len(d)}x{len(d)}).", [path], dataset=dataset,
              n=len(d), mean_offdiag_corr=float(off.mean()), max_offdiag_corr=float(off.max()),
              positive_definite=bool(ev.min() > 0), condition_number=float(ev.max() / ev.min()))


# ================================================================= fitting
def _sample_arrays(dataset: str, zmin, zmax, calibrators) -> dict:
    s = io.load_sample(dataset, _zmin(dataset, zmin), zmax, calibrators)
    if len(s["z"]) < 3:
        raise ValueError("Fewer than 3 SNe pass the cuts.")
    return s


@tier(2)
@validate_call
def compute_sn_chi2(
    dataset: Annotated[Dataset, Field(description="SN sample.")] = "pantheonplus",
    Om: Annotated[float, Field(gt=0.0, lt=1.0, description="Omega_m.")] = 0.3,
    w0: Annotated[float, Field(ge=-3.0, le=1.0, description="w0.")] = -1.0,
    wa: Annotated[float, Field(ge=-10.0, le=5.0, description="wa.")] = 0.0,
    zmin: Annotated[float | None, Field(ge=0.0, le=3.0, description="z cut (None = release standard).")] = None,
    zmax: Annotated[float, Field(gt=0.0, le=10.0, description="Upper z cut.")] = 10.0,
    calibrators: Annotated[Calib, Field(description="include | exclude | only (Pantheon+).")] = "include",
) -> ArtifactResult:
    """Offset-marginalised chi2 of ONE cosmology against a sample (full covariance). Fast.

    Use to test a specific model (e.g. Planck LCDM Om=0.315, or DESI DR2 w0waCDM) against
    a SN sample, or to compute Delta-chi2 between two models without a full fit.
    """
    s = _sample_arrays(dataset, zmin, zmax, calibrators)
    r = s["mu"] - kernels.distance_modulus(s["z"], s["zhel"], om=Om, w0=w0, wa=wa)
    icov = np.linalg.inv(s["cov"])
    c2 = kernels.marg_chi2(r, icov)
    r0 = s["mu"] - kernels.distance_modulus(s["z"], s["zhel"], om=0.3)
    c0 = kernels.marg_chi2(r0, icov)
    return ok(f"chi2 = {c2:.2f} for N = {len(r)} ({_LABELS[dataset]}; Om={Om}, w0={w0}, wa={wa}).", [],
              dataset=dataset, chi2=c2, n_sn=len(r), chi2_per_n=c2 / len(r), delta_chi2_vs_lcdm_om0p3=c2 - c0,
              best_offset_mag=kernels.best_offset(r, icov))


def _write_fit(out: Path, tag: str, res: dict) -> list[Path]:
    files = []
    p = out / f"{tag}_summary.json"
    small = {k: v for k, v in res.items() if k not in ("marginals_1d", "marginals_2d")}
    small["marginals_1d"] = {k: {kk: vv for kk, vv in v.items() if kk not in ("x", "post", "chi2_profile")}
                             for k, v in res["marginals_1d"].items()}
    p.write_text(json.dumps(small, indent=2))
    files.append(p)
    for name, m in res["marginals_1d"].items():
        files.append(_write_csv(out / f"{tag}_1d_{name}.csv", [name, "posterior", "chi2_profile"],
                                zip(np.round(m["x"], 6), np.round(m["post"], 8), np.round(m["chi2_profile"], 4), strict=True)))
    for key, m in res["marginals_2d"].items():
        a, b = key.split("-")
        post = np.asarray(m["post"])
        prof = np.asarray(m["chi2_profile"])
        lv68, lv95 = m["hpd_levels_68_95"]
        rows = []
        for i, x in enumerate(m["x"]):
            for j, y in enumerate(m["y"]):
                rows.append((round(x, 6), round(y, 6), f"{post[i, j]:.6e}", round(float(prof[i, j]), 4),
                             int(post[i, j] >= lv68), int(post[i, j] >= lv95)))
        files.append(_write_csv(out / f"{tag}_2d_{a}_{b}.csv", [a, b, "posterior", "chi2_profile", "in_68", "in_95"], rows))
    return files


@tier(2)
@validate_call
def fit_sn_cosmology(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for summary JSON + marginal/contour CSVs.")],
    dataset: Annotated[Dataset, Field(description="SN sample.")] = "pantheonplus",
    model: Annotated[Model, Field(description="lcdm (flat, Om) | wcdm (flat, Om+w) | w0wacdm (flat CPL, Om+w0+wa).")] = "lcdm",
    zmin: Annotated[float | None, Field(ge=0.0, le=3.0, description="z cut (None = release standard: Pantheon+ 0.01).")] = None,
    zmax: Annotated[float, Field(gt=0.0, le=10.0, description="Upper z cut.")] = 10.0,
    calibrators: Annotated[Calib, Field(description="Pantheon+ calibrators include | exclude | only.")] = "include",
    n_grid: Annotated[int, Field(ge=0, le=801, description="Grid points per parameter (0 = 401/121/41 by model).")] = 0,
    Om_range: Annotated[list[float] | None, Field(description="[lo, hi] flat prior on Om (default [0.05, 0.65]).")] = None,
    w_range: Annotated[list[float] | None, Field(description="[lo, hi] prior on w (wcdm; default [-2, -0.3]).")] = None,
    w0_range: Annotated[list[float] | None, Field(description="[lo, hi] prior on w0 (default [-2, 0]).")] = None,
    wa_range: Annotated[list[float] | None, Field(description="[lo, hi] prior on wa (default [-6, 3]).")] = None,
) -> ArtifactResult:
    """Fit flat LCDM / wCDM / w0waCDM to one SN sample on a grid (+ Nelder-Mead refinement).

    Use for SN-only cosmology: full stat+sys covariance, offset M/H0 marginalised
    analytically, flat priors over the grid. Writes `<tag>_summary.json` (best fit, chi2,
    1D marginal mean/sd/68/95%), `<tag>_1d_<p>.csv` and `<tag>_2d_<a>_<b>.csv` (posterior,
    profile chi2 and in_68/in_95 HPD flags) for plot_sn_contours. Heavy: w0wacdm with
    the default 41^3 grid on Pantheon+ takes a few minutes locally. When set_dispatch
    points at an HPC site it runs as one facility job (kernel sn.kernels.fit_sn_cosmology).
    SN alone barely constrain w0wa; expect posteriors hitting prior edges (flagged).
    Sanity: Pantheon+ lcdm Om = 0.334 +- 0.018 (Brout+22); DES-SN5YR 0.352 +- 0.017.
    """
    ranges = {k: v for k, v in {"Om": Om_range, "w": w_range, "w0": w0_range, "wa": wa_range}.items() if v}
    for k, v in ranges.items():
        if len(v) != 2 or v[0] >= v[1]:
            raise ValueError(f"{k}_range must be [lo, hi] with lo < hi")
    zlo = _zmin(dataset, zmin)
    from mcp_server.dispatch import remote_site, run_kernel  # lazy: server-only

    site = remote_site()
    args = {"model": model, "ranges": ranges or None, "n_grid": n_grid, "refine": True}
    if site:
        res_k = run_kernel("sn.kernels.fit_sn_cosmology",
                           {**args, "source": {"dataset": dataset, "zmin": zlo, "zmax": zmax, "calibrators": calibrators}},
                           pip_deps=DISPATCH_KERNELS["fit_sn_cosmology"]["pip_deps"],
                           duration=DISPATCH_KERNELS["fit_sn_cosmology"]["duration_hint_s"])
        res = res_k["result"]
        computed_on = res_k.get("host", site)
    else:
        s = _sample_arrays(dataset, zmin, zmax, calibrators)
        res = kernels.fit_sn_cosmology(data={"z": s["z"], "zhel": s["zhel"], "mu": s["mu"], "cov": s["cov"]}, **args)
        computed_on = "local"
    res.update(dataset=dataset, zmin=zlo, zmax=zmax, calibrators=calibrators, computed_on=computed_on)
    out = resolve_outdir(output_dir)
    tag = f"fit_{dataset}_{model}_{param_slug({'zmin': zlo, 'zmax': zmax, 'cal': calibrators, 'n': n_grid, 'r': json.dumps(ranges, sort_keys=True)})}"
    files = _write_fit(out, tag, res)
    summ = {p: {k: round(v[k], 4) for k in ("mean", "sd", "lo68", "hi68")} | {"hits_prior_edge": v["hits_prior_edge"]}
            for p, v in res["marginals_1d"].items()}
    txt = ", ".join(f"{p} = {v['mean']:.3f} +- {v['sd']:.3f}" for p, v in summ.items())
    return ok(f"{_LABELS[dataset]} {model}: {txt} (chi2_min = {res['chi2_min']:.1f}, N = {res['n_sn']}).", files,
              dataset=dataset, model=model, best_fit=res["best_fit"], marginals=summ, chi2_min=res["chi2_min"],
              n_sn=res["n_sn"], dof=res["dof"], computed_on=computed_on, summary_json=str(files[0]),
              likelihood=res["likelihood"])


fit_sn_cosmology.weight = "heavy"


def _read_2d(path: str):
    import pandas as pd

    d = pd.read_csv(path)
    a, b = d.columns[:2]
    xs = np.unique(d[a].to_numpy())
    ys = np.unique(d[b].to_numpy())
    post = d["posterior"].to_numpy().reshape(len(xs), len(ys))
    return a, b, xs, ys, post / post.max()


def _levels(post):
    return sorted(kernels._hpd_levels(post))  # [95-level, 68-level]


def _draw_contour(ax, path, color, label, filled=True):
    a, b, xs, ys, post = _read_2d(path)
    lv = _levels(post)
    if filled:
        ax.contourf(xs, ys, post.T, levels=[lv[0], lv[1], 1.01], colors=[color, color], alpha=0.25)
        ax.contourf(xs, ys, post.T, levels=[lv[1], 1.01], colors=[color], alpha=0.35)
    ax.contour(xs, ys, post.T, levels=lv, colors=[color], linewidths=[1.0, 1.6])
    ax.plot([], [], color=color, lw=2, label=label)
    return a, b


@tier(2)
@validate_call
def plot_sn_contours(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the PNG.")],
    contour_csvs: Annotated[list[str], Field(min_length=1, max_length=6, description="2D contour CSVs from fit_sn_cosmology "
                                             "(`*_2d_<a>_<b>.csv`), all for the same parameter pair.")],
    labels: Annotated[list[str] | None, Field(description="Legend labels (default: file stems).")] = None,
    show_lcdm: Annotated[bool, Field(description="Mark w = -1 / (w0, wa) = (-1, 0).")] = True,
) -> ArtifactResult:
    """Plot 68%/95% credible contours (HPD of the gridded marginal posterior) for one parameter pair.

    Use after fit_sn_cosmology (wcdm: Om-w; w0wacdm: w0-wa or Om-w0) to reproduce the
    familiar SN contour figures (Brout+22 Fig. 9, DES 2024 Fig. 4, Rubin+23 Fig. 9);
    pass several CSVs to overlay samples.
    """
    plt = mpl()
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    cols = ["tab:blue", "tab:red", "tab:green", "tab:purple", "tab:orange", "tab:brown"]
    names = None
    for i, pth in enumerate(contour_csvs):
        if not Path(pth).is_file():
            raise ValueError(f"No such file: {pth}")
        stem = Path(pth).stem
        ds = next((k for k in _COLORS if f"fit_{k}_" in stem), None)
        col = _COLORS.get(ds, cols[i % len(cols)]) if len(set(contour_csvs)) == len(contour_csvs) else cols[i]
        lab = labels[i] if labels and i < len(labels) else (_LABELS.get(ds, stem) if ds else stem)
        pair = _draw_contour(ax, pth, col, lab)
        if names and pair != names:
            raise ValueError(f"All CSVs must share a parameter pair; got {names} and {pair}")
        names = pair
    a, b = names
    if show_lcdm:
        if b in ("w", "w0"):
            ax.axhline(-1, color="k", ls="--", lw=0.8)
        if (a, b) == ("w0", "wa"):
            ax.axvline(-1, color="k", ls="--", lw=0.8)
            ax.axhline(0, color="k", ls="--", lw=0.8)
    ax.set_xlabel(_TEX.get(a, a))
    ax.set_ylabel(_TEX.get(b, b))
    ax.legend(fontsize=9, frameon=False, loc="best")
    ax.set_title("SN Ia only: 68% / 95% credible regions", fontsize=10)
    path = save_fig(fig, resolve_outdir(output_dir), f"contours_{a}_{b}_{param_slug({'f': '|'.join(contour_csvs)})}")
    return ok(f"Contours of {a}-{b} for {len(contour_csvs)} fit(s).", [path], params=[a, b], n_fits=len(contour_csvs))


# ================================================================= comparison
def _moments_2d(path: str):
    a, b, xs, ys, post = _read_2d(path)
    w = post / post.sum()
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    m = np.array([(X * w).sum(), (Y * w).sum()])
    dx, dy = X - m[0], Y - m[1]
    c = np.array([[(dx * dx * w).sum(), (dx * dy * w).sum()], [(dx * dy * w).sum(), (dy * dy * w).sum()]])
    return (a, b), m, c


def _gauss_tension(m1, c1, m2, c2) -> dict:
    """Gaussian parameter-difference tension, ignoring sample overlap (see notes)."""
    from scipy import stats

    d = np.atleast_1d(m1 - m2)
    C = np.atleast_2d(c1 + c2)
    chi2 = float(d @ np.linalg.solve(C, d))
    p = float(stats.chi2.sf(chi2, len(d)))
    return {"chi2": chi2, "dof": len(d), "p_value": p, "sigma_equiv": float(stats.norm.isf(p / 2))}


@tier(2)
@validate_call
def compare_sn_samples(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for PNG + CSV/JSON outputs.")],
    datasets: Annotated[list[Dataset], Field(min_length=2, max_length=4, description="Samples to compare.")] =
        ["pantheonplus", "des5yr", "union3"],  # noqa: B006 - pydantic copies defaults
    mode: Annotated[Literal["residuals", "contours"], Field(description="residuals: overlay binned Hubble residuals "
                                                            "vs one fiducial (fast). contours: fit each sample "
                                                            "(model) and overlay 68/95% contours + tension (heavy).")] = "residuals",
    model: Annotated[Model, Field(description="Model for mode=contours.")] = "wcdm",
    fid_Om: Annotated[float, Field(gt=0.0, lt=1.0, description="Fiducial flat-LCDM Om for residuals.")] = 0.3,
    n_bins: Annotated[int, Field(ge=3, le=40, description="Bins per sample for residuals.")] = 15,
    n_grid: Annotated[int, Field(ge=0, le=801, description="Grid per parameter for contours (0 = model default).")] = 0,
) -> ArtifactResult:
    """Compare SN samples: binned residual Hubble diagrams, or cosmology contours + tension.

    Use to see whether Pantheon+, DES-SN5YR and Union3 agree. `residuals` overlays each
    sample's full-covariance binned mu residuals vs the SAME flat LCDM fiducial (each sample's
    offset removed separately - only shapes are comparable). `contours` runs fit_sn_cosmology
    per sample, overlays the 2D contours (Om-w for wcdm, w0-wa for w0wacdm, 1D for lcdm) and
    reports pairwise Gaussian tensions. CAUTION: the samples share SNe (low-z CfA/CSP/Foundation;
    Union3 contains much of Pantheon+), so they are NOT independent. The quoted tension uses
    Cov(diff) = C1 + C2, ignoring the positive cross-correlation, so it UNDER-states the
    significance of a real difference - treat it as indicative, not a formal test.
    """
    out = resolve_outdir(output_dir)
    plt = mpl()
    if mode == "residuals":
        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        rows, files, meta = [], [], {}
        for k, ds in enumerate(datasets):
            cal = "exclude"
            s = io.load_sample(ds, _zmin(ds, None), 10.0, cal)
            resid, off, icov = _residuals(s, fid_Om, -1.0, 0.0)
            if len(s["z"]) > n_bins:
                edges = np.geomspace(s["z"].min() * 0.999, s["z"].max() * 1.001, n_bins + 1)
                zb, rb, covb, nb, _ = _bin(s, resid, icov, edges)
            else:
                zb, rb, covb, nb = s["z"], resid, s["cov"], np.ones(len(s["z"]), int)
            eb = np.sqrt(np.diag(covb))
            ax.errorbar(zb * (1 + 0.012 * (k - 1)), rb, yerr=eb, fmt="os^D"[k], ms=5, color=_COLORS[ds], lw=1.2,
                        capsize=2, label=f"{_LABELS[ds]} (N={len(s['z'])})")
            rows += [(ds, round(z, 5), int(n), round(r, 5), round(e, 5)) for z, n, r, e in zip(zb, nb, rb, eb, strict=True)]
            hi = zb > 0.5
            meta[ds] = {"n": len(s["z"]), "offset_mag": off,
                        "mean_resid_z_gt_0p5": float(np.average(rb[hi], weights=1 / eb[hi] ** 2)) if hi.any() else None}
        zz = np.geomspace(0.01, 2.3, 300)
        for om, ls in ((fid_Om - 0.05, "--"), (fid_Om + 0.05, ":")):
            ax.plot(zz, kernels.distance_modulus(zz, om=om) - kernels.distance_modulus(zz, om=fid_Om), "k" + ls, lw=0.8,
                    label=f"flat $\\Lambda$CDM $\\Omega_m$={om:.2f}")
        ax.axhline(0, color="k", lw=1)
        ax.set_xscale("log")
        ax.set_xlabel(r"$z$")
        ax.set_ylabel(r"$\mu-\mu_{\Lambda{\rm CDM}}(\Omega_m=%g)$ [mag]" % fid_Om)
        ax.set_ylim(-0.25, 0.25)
        ax.legend(fontsize=8, frameon=False, ncol=2)
        ax.set_title("Binned Hubble residuals (offset removed per sample; full covariance)", fontsize=10)
        files.append(save_fig(fig, out, f"compare_residuals_{'_'.join(datasets)}"))
        files.append(_write_csv(out / f"compare_residuals_{'_'.join(datasets)}.csv",
                                ["dataset", "z", "n_sn", "dmu", "dmu_err"], rows))
        return ok(f"Binned residuals of {', '.join(_LABELS[d] for d in datasets)} vs flat LCDM Om={fid_Om}.", files,
                  per_dataset=meta)

    fits, files = {}, []
    for ds in datasets:
        r = fit_sn_cosmology(output_dir=str(out), dataset=ds, model=model, n_grid=n_grid)
        fits[ds] = r
        files += r.files
    names = kernels.MODELS[model]
    tension, summary = {}, {}
    if model == "lcdm":
        fig, ax = plt.subplots(figsize=(6, 4))
        for ds, r in fits.items():
            f1 = next(f for f in r.files if f.endswith("_1d_Om.csv"))
            d = np.loadtxt(f1, delimiter=",", skiprows=1)
            ax.plot(d[:, 0], d[:, 1], color=_COLORS[ds], lw=2,
                    label=f"{_LABELS[ds]}: {r.metadata['marginals']['Om']['mean']:.3f} $\\pm$ {r.metadata['marginals']['Om']['sd']:.3f}")
            summary[ds] = r.metadata["marginals"]
        ax.set_xlabel(r"$\Omega_m$ (flat $\Lambda$CDM)")
        ax.set_ylabel("posterior (peak = 1)")
        ax.set_xlim(0.15, 0.55)
        ax.legend(fontsize=9, frameon=False)
        moments = {ds: (np.array([r.metadata["marginals"]["Om"]["mean"]]),
                        np.array([[r.metadata["marginals"]["Om"]["sd"] ** 2]])) for ds, r in fits.items()}
        pair = "Om"
    else:
        pair = "Om-w" if model == "wcdm" else "w0-wa"
        a, b = pair.split("-")
        fig, ax = plt.subplots(figsize=(5.6, 5.0))
        moments = {}
        for ds, r in fits.items():
            f2 = next(f for f in r.files if f.endswith(f"_2d_{a}_{b}.csv"))
            _draw_contour(ax, f2, _COLORS[ds], _LABELS[ds])
            _, m, c = _moments_2d(f2)
            moments[ds] = (m, c)
            summary[ds] = r.metadata["marginals"]
        ax.axhline(-1 if b in ("w", "w0") else 0, color="k", ls="--", lw=0.8)
        if a == "w0":
            ax.axvline(-1, color="k", ls="--", lw=0.8)
        ax.set_xlabel(_TEX[a])
        ax.set_ylabel(_TEX[b])
        ax.legend(fontsize=9, frameon=False)
    ax.set_title(f"SN-only {model}: sample comparison", fontsize=10)
    files.append(save_fig(fig, out, f"compare_{model}_{'_'.join(datasets)}"))
    keys = list(moments)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            (m1, c1), (m2, c2) = moments[keys[i]], moments[keys[j]]
            tension[f"{keys[i]} vs {keys[j]}"] = _gauss_tension(m1, c1, m2, c2) | {"params": pair}
    tpath = out / f"compare_{model}_{'_'.join(datasets)}_tension.json"
    tpath.write_text(json.dumps({"model": model, "marginals": summary, "tension": tension,
                                 "note": "Gaussian approximation from gridded marginal mean/covariance; samples "
                                         "overlap (shared low-z SNe), so these are NOT independent."}, indent=2))
    files.append(tpath)
    txt = "; ".join(f"{k}: {v['sigma_equiv']:.1f} sigma" for k, v in tension.items())
    return ok(f"{model} contours for {', '.join(_LABELS[d] for d in datasets)}. Pairwise ({pair}) tension: {txt}.", files,
              model=model, parameters=names, marginals=summary, tension=tension)


compare_sn_samples.weight = "heavy"


# ================================================================= tier-1 tools

Salt2Dataset = Literal["pantheonplus", "des5yr"]
LCDataset = Literal["pantheonplus", "des5yr"]


@tier(1)
@validate_call
def fetch_sn_salt2_params(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the CSV.")],
    dataset: Annotated[Salt2Dataset, Field(description="pantheonplus or des5yr.")] = "pantheonplus",
    zmin: Annotated[float, Field(ge=0.0, le=3.0, description="Keep z > zmin.")] = 0.0,
    zmax: Annotated[float, Field(gt=0.0, le=10.0, description="Keep z <= zmax.")] = 10.0,
) -> ArtifactResult:
    """Fetch per-SN SALT2/SALT3 light-curve fit parameters: x1, c, mB and their errors and covariances,
    HOST_LOGMASS and bias corrections.

    Use this to get the TIER-1 (catalog-level) data needed for Tripp standardisation via
    standardize_sn_salt2. Pantheon+ uses SALT2 fits (Brout+22); DES-SN5YR uses SALT3 (Kenworthy+21).
    The output CSV has one row per SN and includes the released distance modulus for comparison.
    """
    d = io.load_salt2_params(dataset)
    z = d["z"]
    m = (z > zmin) & (z <= zmax)
    out = resolve_outdir(output_dir)
    cols = ["cid", "z", "zhel", "mB", "mBERR", "x1", "x1ERR", "c", "cERR",
            "COV_x1_c", "COV_mB_x1", "COV_mB_c", "HOST_LOGMASS", "HOST_LOGMASS_ERR",
            "biasCor_mB", "mu_released", "survey", "is_calibrator"]
    rows = []
    for i in range(len(z)):
        if not m[i]:
            continue
        row = []
        for c in cols:
            v = d[c][i] if c in d else ""
            if isinstance(v, (np.floating, float)):
                row.append(round(float(v), 6))
            elif isinstance(v, (np.integer, int)):
                row.append(int(v))
            else:
                row.append(str(v))
        rows.append(row)
    path = _write_csv(out / f"{dataset}_salt2_params.csv", cols, rows)
    from ..registry import DATASETS
    return ok(f"{len(rows)} {_LABELS[dataset]} SALT2/SALT3 fit parameters (z in ({zmin}, {zmax}]).",
              [path], dataset=dataset, n_sn=len(rows), salt_version=d["salt_version"],
              columns=cols, citation=DATASETS[dataset]["citation"])


@tier(1)
@validate_call
def standardize_sn_salt2(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the Hubble diagram CSV + summary JSON.")],
    dataset: Annotated[Salt2Dataset, Field(description="pantheonplus or des5yr.")] = "pantheonplus",
    zmin: Annotated[float, Field(ge=0.0, le=3.0, description="Keep z > zmin (Pantheon+ cosmology cut is 0.01).")] = 0.01,
    zmax: Annotated[float, Field(gt=0.0, le=10.0, description="Keep z <= zmax.")] = 10.0,
    exclude_calibrators: Annotated[bool, Field(description="Drop Cepheid-host calibrators (IS_CALIBRATOR=1).")] = True,
    Om: Annotated[float, Field(gt=0.0, lt=1.0, description="Fiducial flat-LCDM Om for the distance prior.")] = 0.3,
) -> ArtifactResult:
    """Tripp standardisation: build a Hubble diagram from per-SN SALT2 fits (TIER-1 → TIER-2).

    Fits alpha (stretch), beta (colour), gamma (host-mass step at 10^10 Msun) and sigma_int
    by maximum likelihood: mu = mB - M + alpha*x1 - beta*c ± gamma/2 - biasCor.
    Literature: Pantheon+ alpha~0.15, beta~3.0, gamma~0.05; DES-SN5YR alpha=0.161, beta=3.12, gamma=0.038.
    The tool compares the resulting mu with the released Hubble diagram (rms_vs_released_mu should be
    ~few mmag for Pantheon+ and ~20 mmag for DES, since we use a diagonal approximation).
    Uses ONLY tier-1 inputs (SALT2 fit parameters); never tier-3 posteriors.
    """
    d = io.load_salt2_params(dataset)
    z = d["z"]
    m = (z > zmin) & (z <= zmax)
    if exclude_calibrators and "is_calibrator" in d:
        m &= ~(np.asarray(d["is_calibrator"]).astype(bool))
    # Build dict of arrays for the kernel
    keys = ["z", "zhel", "mB", "mBERR", "x1", "x1ERR", "c", "cERR",
            "COV_x1_c", "COV_mB_x1", "COV_mB_c", "HOST_LOGMASS", "biasCor_mB", "VPECERR"]
    data = {}
    for k in keys:
        if k in d:
            arr = np.asarray(d[k])
            data[k] = arr[m].tolist()
    data["mu_released"] = np.asarray(d["mu_released"])[m].tolist()
    cids = [c for c, keep in zip(d["cid"], m, strict=False) if keep]

    res = kernels.standardize_sn_salt2(data, Om=Om)

    out = resolve_outdir(output_dir)
    tag = f"{dataset}_standardized"
    # Write CSV
    hd_rows = list(zip(cids, res["z"], res["zhel"], res["mu"], res["mu_err"], strict=True))
    csv_path = _write_csv(out / f"{tag}_hd.csv",
                          ["cid", "z", "zhel", "mu", "mu_err"],
                          [(c, round(z, 6), round(zh, 6), round(mu, 5), round(me, 5))
                           for c, z, zh, mu, me in hd_rows])
    # Write summary JSON
    summ = {k: v for k, v in res.items() if k not in ("mu", "mu_err", "z", "zhel")}
    summ_path = out / f"{tag}_summary.json"
    summ_path.write_text(json.dumps(summ, indent=2))

    from ..registry import DATASETS
    txt = (f"alpha={res['alpha']:.3f}, beta={res['beta']:.3f}, gamma={res['gamma']:.3f}, "
           f"sigma_int={res['sigma_int']:.4f}")
    return ok(f"{_LABELS[dataset]} Tripp standardisation ({res['n_sn']} SNe): {txt}.",
              [csv_path, summ_path], dataset=dataset, **summ,
              citation=DATASETS[dataset]["citation"])


standardize_sn_salt2.weight = "heavy"


@tier(1)
@validate_call
def fetch_sn_lightcurve(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the light-curve CSV(s).")],
    cid: Annotated[str, Field(min_length=1, description="SN candidate ID (e.g. '2011fe', '523183', '2254').")],
    dataset: Annotated[LCDataset, Field(description="pantheonplus (SNANA text files, per-SN download) "
                                        "or des5yr (extracted from DES5YR_SALT3_LCFIT.LCPLOT.gz).")] = "pantheonplus",
) -> ArtifactResult:
    """Download and parse raw photometric light-curve data for one SN Ia (tier-1).

    Pantheon+ light curves are individual SNANA .dat files from PantheonPlusSH0ES/DataRelease
    (downloaded on demand). DES-SN5YR light curves come from the SALT3 fit LCPLOT file (data +
    model columns; 7.4 MB download on first use). Returns a CSV with mjd, band, flux, flux_err
    (and mag where available). Use plot_sn_lightcurve to visualise.
    """
    out = resolve_outdir(output_dir)
    files = []
    if dataset == "pantheonplus":
        lcs = io.fetch_pantheonplus_lc(cid)
        for lc in lcs:
            tag = f"lc_{dataset}_{cid}_{lc['folder']}"
            rows = list(zip(lc["mjd"], lc["band"], lc["flux"], lc["flux_err"],
                            lc["mag"], lc["mag_err"], strict=True))
            path = _write_csv(out / f"{tag}.csv", ["mjd", "band", "flux", "flux_err", "mag", "mag_err"],
                              [(round(m, 4), b, round(f, 4), round(fe, 4),
                                round(mg, 4) if np.isfinite(mg) else "",
                                round(me, 4) if np.isfinite(me) else "")
                               for m, b, f, fe, mg, me in rows])
            files.append(path)
        meta = {"snid": lcs[0]["snid"], "survey": lcs[0]["survey"], "z_hel": lcs[0]["z_hel"],
                "z_cmb": lcs[0]["z_cmb"], "filters": lcs[0]["filters"], "peak_mjd": lcs[0]["peak_mjd"],
                "n_obs": sum(len(lc["mjd"]) for lc in lcs), "n_files": len(lcs),
                "folders": [lc["folder"] for lc in lcs]}
    elif dataset == "des5yr":
        lc = io.fetch_des5yr_lc(cid)
        if lc is None:
            raise ValueError(f"CID {cid!r} not found in DES-SN5YR LCPLOT file.")
        data_mask = lc["data_model"] == 0
        tag = f"lc_{dataset}_{cid}"
        rows = list(zip(lc["mjd"][data_mask], lc["band"][data_mask],
                        lc["flux"][data_mask], lc["flux_err"][data_mask], strict=True))
        path = _write_csv(out / f"{tag}.csv", ["mjd", "band", "flux", "flux_err"],
                          [(round(m, 4), b, round(f, 4), round(fe, 4)) for m, b, f, fe in rows])
        files.append(path)
        meta = {"snid": cid, "n_obs": int(data_mask.sum()),
                "bands": list(np.unique(lc["band"][data_mask]))}
    else:
        raise ValueError(f"dataset must be 'pantheonplus' or 'des5yr', got {dataset!r}")

    return ok(f"Light curve for SN {cid} ({_LABELS[dataset]}): {meta.get('n_obs', '?')} observations.",
              files, dataset=dataset, **meta)


@tier(1)
@validate_call
def plot_sn_lightcurve(
    output_dir: Annotated[str, Field(min_length=1, description="Directory for the PNG.")],
    cid: Annotated[str, Field(min_length=1, description="SN candidate ID.")],
    dataset: Annotated[LCDataset, Field(description="pantheonplus or des5yr.")] = "pantheonplus",
    flux_space: Annotated[bool, Field(description="Plot in flux (True) or magnitude (False, inverted axis).")] = True,
) -> ArtifactResult:
    """Plot a SN Ia multi-band light curve (tier-1 photometry).

    Downloads the light curve on demand (same as fetch_sn_lightcurve) and makes a publication-style
    plot with one colour per band. Flux space (default) shows the SALT2 flux calibration; magnitude
    space inverts the y-axis and uses mag where available.
    """
    plt = mpl()
    out = resolve_outdir(output_dir)
    band_colors = {"u": "indigo", "g": "green", "r": "red", "i": "darkorange", "z": "brown", "Y": "darkred",
                   "B": "blue", "V": "green", "R": "red", "I": "darkorange",
                   "U": "violet", "J": "saddlebrown", "H": "olive", "K": "black",
                   "c": "cyan", "o": "orange", "w": "gray", "G": "green"}

    if dataset == "pantheonplus":
        lcs = io.fetch_pantheonplus_lc(cid)
        fig, ax = plt.subplots(figsize=(8, 5))
        peak = lcs[0]["peak_mjd"]
        for lc in lcs:
            bands = np.unique(lc["band"])
            for b in bands:
                m = lc["band"] == b
                color = band_colors.get(b, None)
                phase = lc["mjd"][m] - (peak if peak else lc["mjd"][m].min())
                if flux_space:
                    ax.errorbar(phase, lc["flux"][m], yerr=lc["flux_err"][m],
                                fmt="o", ms=3, color=color, label=f"{b} ({lc['folder']})", alpha=0.8)
                else:
                    good = np.isfinite(lc["mag"][m])
                    if good.any():
                        ax.errorbar(phase[good], lc["mag"][m][good], yerr=lc["mag_err"][m][good],
                                    fmt="o", ms=3, color=color, label=f"{b} ({lc['folder']})", alpha=0.8)
        ylabel = "flux (FLUXCAL)" if flux_space else "magnitude"
        title_extra = f"z={lcs[0]['z_hel']:.4f}" if lcs[0]["z_hel"] else ""
    elif dataset == "des5yr":
        lc = io.fetch_des5yr_lc(cid)
        if lc is None:
            raise ValueError(f"CID {cid!r} not found in DES-SN5YR LCPLOT.")
        fig, ax = plt.subplots(figsize=(8, 5))
        bands = np.unique(lc["band"])
        data_mask = lc["data_model"] == 0
        model_mask = lc["data_model"] == 1
        for b in bands:
            bm = lc["band"] == b
            color = band_colors.get(b, None)
            dm = bm & data_mask
            if dm.any():
                ax.errorbar(lc["phase"][dm], lc["flux"][dm], yerr=lc["flux_err"][dm],
                            fmt="o", ms=3, color=color, label=f"{b} data", alpha=0.8)
            mm = bm & model_mask
            if mm.any():
                order = np.argsort(lc["phase"][mm])
                ax.plot(lc["phase"][mm][order], lc["flux"][mm][order], "-", color=color, alpha=0.5, lw=1)
        ylabel = "flux (FLUXCAL)"
        title_extra = ""
        flux_space = True  # DES LCPLOT is always in flux
    else:
        raise ValueError("dataset must be 'pantheonplus' or 'des5yr'")

    ax.set_xlabel("phase [days from peak]" if dataset == "pantheonplus" else "phase [days]")
    ax.set_ylabel(ylabel)
    if not flux_space:
        ax.invert_yaxis()
    ax.legend(fontsize=7, frameon=False, ncol=3, loc="best")
    ax.set_title(f"SN {cid} ({_LABELS[dataset]}) {title_extra}", fontsize=11)
    if flux_space:
        ax.axhline(0, color="k", lw=0.5, ls=":")
    path = save_fig(fig, out, f"lc_{dataset}_{cid}")
    return ok(f"Light curve plot for SN {cid} ({_LABELS[dataset]}).", [path], dataset=dataset, cid=cid)
