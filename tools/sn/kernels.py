"""Dispatchable SN-cosmology kernels (JSON in -> JSON out).

Imports ONLY numpy + scipy (declared pip_deps), so this file runs unchanged
on a facility node via mcp_server.dispatch.run_kernel. Data arrive either as
inline arrays (z, zhel, mu, cov / icov) or as a ``source`` dict naming a
public download URL - never as local paths.

Likelihood: Gaussian in distance modulus with the full covariance, with the
additive offset (M, or equivalently 5 log10 H0) marginalised ANALYTICALLY
(Goliath et al. 2001, A9-12; the DES-SN5YR / Pantheon+ cosmosis convention):

    chi2_marg = r^T C^-1 r - (1^T C^-1 r)^2 / (1^T C^-1 1),   r = mu_obs - mu_th
"""

from __future__ import annotations

import math

import numpy as np

C_KMS = 299792.458

MODELS = {
    "lcdm": ["Om"],            # flat LCDM
    "wcdm": ["Om", "w"],       # flat wCDM
    "w0wacdm": ["Om", "w0", "wa"],  # flat CPL
}


# ------------------------------------------------------------------ theory
def comoving_distance_c_over_h0(z, om: float, w0: float = -1.0, wa: float = 0.0, ok: float = 0.0):
    """Line-of-sight comoving distance in units of c/H0 (flat unless ok != 0).

    Radiation is neglected (negligible for z < 3 at the 1e-4 level in mu).
    Dark energy follows CPL: w(a) = w0 + wa (1 - a).
    """
    from scipy.integrate import cumulative_trapezoid

    z = np.atleast_1d(np.asarray(z, dtype=float))
    zmax = max(float(z.max()), 1e-3) * 1.001
    zg = np.concatenate([np.linspace(0.0, min(zmax, 0.1), 400), np.linspace(min(zmax, 0.1), zmax, 3600)[1:]])
    a = 1.0 / (1.0 + zg)
    ode = 1.0 - om - ok
    de = ode * (1.0 + zg) ** (3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))
    e2 = om * (1.0 + zg) ** 3 + ok * (1.0 + zg) ** 2 + de
    if np.any(e2 <= 0):
        return np.full_like(z, np.nan)
    dc = cumulative_trapezoid(1.0 / np.sqrt(e2), zg, initial=0.0)
    return np.interp(z, zg, dc)


def transverse_distance(dc, ok: float):
    if abs(ok) < 1e-8:
        return dc
    s = math.sqrt(abs(ok))
    return np.sinh(s * dc) / s if ok > 0 else np.sin(s * dc) / s


def distance_modulus(z, zhel=None, om: float = 0.3, w0: float = -1.0, wa: float = 0.0,
                     h0: float = 70.0, ok: float = 0.0):
    """mu = 5 log10(D_L / 10 pc) with D_L = (1 + z_hel) D_M(z_HD) (Pantheon+/DES convention)."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    zhel = z if zhel is None else np.atleast_1d(np.asarray(zhel, dtype=float))
    dm = transverse_distance(comoving_distance_c_over_h0(z, om, w0, wa, ok), ok)
    dl_mpc = (1.0 + zhel) * dm * C_KMS / h0
    return 5.0 * np.log10(dl_mpc) + 25.0


# ------------------------------------------------------------- likelihood
def marg_chi2(resid, icov, icov_rowsum=None, icov_sum=None) -> float:
    """chi2 with the constant offset analytically marginalised."""
    rs = icov.sum(axis=1) if icov_rowsum is None else icov_rowsum
    cs = float(icov.sum()) if icov_sum is None else icov_sum
    a = float(resid @ icov @ resid)
    b = float(rs @ resid)
    return a - b * b / cs


def best_offset(resid, icov) -> float:
    """Maximum-likelihood additive offset (mu_obs - mu_th) - the value that is marginalised."""
    rs = icov.sum(axis=1)
    return float(rs @ resid / rs.sum())


def _model_mu(z, zhel, model: str, p: dict):
    if model == "lcdm":
        return distance_modulus(z, zhel, om=p["Om"])
    if model == "wcdm":
        return distance_modulus(z, zhel, om=p["Om"], w0=p["w"])
    return distance_modulus(z, zhel, om=p["Om"], w0=p["w0"], wa=p["wa"])


# ------------------------------------------------------------- data access
# Public, verified sources (also listed in tools/registry.py). Kernels may
# download these themselves when dispatched; locally tools.sn.io resolves them
# through tools.common.fetch_file (bundled data/ first, then the cache).
SOURCES = {
    "pantheonplus": {
        "hd": "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
              "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat",
        "cov": "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
               "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov",
        "salt2_params": "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
                        "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat",
        "lc_base": "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
                   "Pantheon%2B_Data/1_DATA/photometry/",
    },
    "des5yr": {
        "hd": "https://raw.githubusercontent.com/des-science/DES-SN5YR/1.3/4_DISTANCES_COVMAT/DES-SN5YR_HD.csv",
        "cov": "https://raw.githubusercontent.com/des-science/DES-SN5YR/1.3/4_DISTANCES_COVMAT/STAT%2BSYS.txt.gz",
        "salt2_params": "https://raw.githubusercontent.com/des-science/DES-SN5YR/1.3/"
                        "4_DISTANCES_COVMAT/DES-SN5YR_HD%2BMetaData.csv",
        "lcplot": "https://raw.githubusercontent.com/des-science/DES-SN5YR/1.3/"
                  "0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz",
    },
    "union3": {
        "hd": "https://raw.githubusercontent.com/rubind/union3_release/main/mu_mat_union3_cosmo%3D2_mu.fits",
    },
    "des_dovekie": {
        "hd": "https://raw.githubusercontent.com/des-science/DES-SN5YR/main/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv",
        "cov": "https://raw.githubusercontent.com/des-science/DES-SN5YR/main/4_DISTANCES_COVMAT/STAT%2BSYS.npz",
    },
}

# ------------- light-curve file index (Pantheon+) bundled as CSV
LC_INDEX_COLS = ("cid", "folder", "path")  # data/pantheonplus/lc_index.csv


def parse_snana_lc(text: str) -> dict:
    """Parse an SNANA-format SN light-curve text file into arrays.

    Returns dict(snid, survey, ra, dec, z_hel, z_cmb, mwebv, peak_mjd, filters,
                 mjd, band, flux, flux_err[, mag, mag_err]).
    """
    import re
    meta = {}
    obs = []
    varlist = None
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("END"):
            break
        if line.startswith("VARLIST:"):
            varlist = line.split(":")[1].split()
            continue
        if line.startswith("OBS:"):
            parts = line.split()[1:]
            if varlist:
                row = {varlist[i]: parts[i] for i in range(min(len(varlist), len(parts)))}
            else:
                row = {"raw": parts}
            obs.append(row)
            continue
        m = re.match(r"^(\w[\w\(\)]*)\s*:\s*(.+)", line)
        if m:
            key = m.group(1).strip()
            val = m.group(2).strip()
            if key not in meta:
                meta[key] = val

    # Extract metadata
    def _float(k, default=None):
        v = meta.get(k, "")
        try:
            return float(v.split()[0])
        except (ValueError, IndexError):
            return default

    snid = meta.get("SNID", "").strip()
    survey = meta.get("SURVEY", "").strip()
    ra = _float("RA")
    dec = _float("DECL", _float("DEC"))
    z_hel = _float("REDSHIFT_HELIO")
    z_cmb = _float("REDSHIFT_CMB")
    mwebv = _float("MWEBV", 0.0)
    peak_mjd = _float("SEARCH_PEAKMJD", _float("PEAKMJD"))
    filters = meta.get("FILTERS", "").strip()

    # Parse obs
    mjd, band, flux, flux_err, mag, mag_err = [], [], [], [], [], []
    for row in obs:
        if isinstance(row, dict) and "MJD" in row:
            _mjd = float(row["MJD"])
            if _mjd < 0:
                continue
            mjd.append(_mjd)
            band.append(row.get("FLT", row.get("BAND", "?")))
            flux.append(float(row.get("FLUXCAL", 0)))
            flux_err.append(float(row.get("FLUXCALERR", 0)))
            _mag = row.get("MAG")
            if _mag is not None and _mag != "NULL":
                try:
                    mag.append(float(_mag))
                    mag_err.append(float(row.get("MAGERR", 0)))
                except ValueError:
                    mag.append(np.nan)
                    mag_err.append(np.nan)
            else:
                mag.append(np.nan)
                mag_err.append(np.nan)

    return {
        "snid": snid, "survey": survey, "ra": ra, "dec": dec,
        "z_hel": z_hel, "z_cmb": z_cmb, "mwebv": mwebv, "peak_mjd": peak_mjd,
        "filters": filters,
        "mjd": np.array(mjd), "band": np.array(band),
        "flux": np.array(flux), "flux_err": np.array(flux_err),
        "mag": np.array(mag), "mag_err": np.array(mag_err),
    }


# -------------------------------------------------------- SALT2 standardization kernel
def standardize_sn_salt2(data: dict, Om: float = 0.3) -> dict:
    """Tripp standardisation of per-SN SALT2 light-curve fits to a Hubble diagram (tier-1 → tier-2).

    Fits alpha, beta, gamma (mass step at 10^10 Msun) and sigma_int by maximum likelihood,
    builds mu = mB - M + alpha*x1 - beta*c + gamma*step, and compares with the released mu.

    data = {"z": [...], "zhel": [...], "mB": [...], "mBERR": [...], "x1": [...], "x1ERR": [...],
            "c": [...], "cERR": [...], "COV_x1_c": [...], "COV_mB_x1": [...], "COV_mB_c": [...],
            "HOST_LOGMASS": [...], "biasCor_mB": [...] (optional),
            "VPECERR": [...] (optional, km/s, default 250),
            "mu_released": [...] (optional, for comparison)}

    Returns: alpha, beta, gamma, sigma_int, M (per-redshift-bin offset), mu, mu_err,
             residuals_vs_released (if mu_released given), chi2_per_dof.
    """
    from scipy.optimize import minimize as sp_minimize

    z = np.asarray(data["z"], float)
    zhel = np.asarray(data.get("zhel", data["z"]), float)
    mB = np.asarray(data["mB"], float)
    mBERR = np.asarray(data["mBERR"], float)
    x1 = np.asarray(data["x1"], float)
    x1ERR = np.asarray(data["x1ERR"], float)
    c = np.asarray(data["c"], float)
    cERR = np.asarray(data["cERR"], float)
    cx1c = np.asarray(data["COV_x1_c"], float)
    cmBx1 = np.asarray(data.get("COV_mB_x1", np.zeros(len(z))), float)
    cmBc = np.asarray(data.get("COV_mB_c", np.zeros(len(z))), float)
    mass = np.asarray(data["HOST_LOGMASS"], float)
    bc = np.asarray(data.get("biasCor_mB", np.zeros(len(z))), float)
    vpec_err = np.asarray(data.get("VPECERR", np.full(len(z), 250.0)), float)
    mu_rel = np.asarray(data["mu_released"], float) if "mu_released" in data else None

    n = len(z)
    step = np.where(mass > 10.0, 0.5, -0.5)
    # peculiar velocity error in mag
    sz = 5.0 / np.log(10.0) * vpec_err / C_KMS / np.maximum(z, 1e-4)

    # theory mu at fiducial Om
    mu_th = distance_modulus(z, zhel, om=Om)

    # Use 10 redshift bins for the offset (approximation of M marginalisation)
    n_bins = min(10, max(1, n // 30))
    edges = np.quantile(np.log(np.maximum(z, 1e-6)), np.linspace(0, 1, n_bins + 1))
    edges[0] -= 1
    edges[-1] += 1
    b = np.clip(np.searchsorted(edges, np.log(np.maximum(z, 1e-6))) - 1, 0, n_bins - 1)

    def nll(p):
        alpha, beta, gamma, log_sig = p
        mu_obs = mB + alpha * x1 - beta * c + gamma * step - bc
        var = (mBERR ** 2 + alpha ** 2 * x1ERR ** 2 + beta ** 2 * cERR ** 2
               + 2 * alpha * cmBx1 - 2 * beta * cmBc - 2 * alpha * beta * cx1c
               + np.exp(2 * log_sig) + sz ** 2)
        var = np.maximum(var, 1e-10)  # clamp to avoid log(negative)
        r = mu_obs - mu_th
        w = 1.0 / var
        # per-bin offset
        off = np.bincount(b, w * r, n_bins) / np.maximum(np.bincount(b, w, n_bins), 1e-30)
        r = r - off[b]
        return float(np.sum(r * r / var + np.log(var)))

    res = sp_minimize(nll, [0.15, 3.0, 0.05, np.log(0.1)], method="BFGS")
    alpha, beta, gamma, log_sig = res.x
    sigma_int = float(np.exp(log_sig))

    # Rebuild mu at the best fit
    mu_obs = mB + alpha * x1 - beta * c + gamma * step - bc
    var = (mBERR ** 2 + alpha ** 2 * x1ERR ** 2 + beta ** 2 * cERR ** 2
           + 2 * alpha * cmBx1 - 2 * beta * cmBc - 2 * alpha * beta * cx1c
           + sigma_int ** 2 + sz ** 2)
    var = np.maximum(var, 1e-10)
    w = 1.0 / var
    off = np.bincount(b, w * (mu_obs - mu_th), n_bins) / np.maximum(np.bincount(b, w, n_bins), 1e-30)
    mu_standardized = mu_obs - off[b]
    mu_err = np.sqrt(var)

    chi2_resid = float(np.sum(((mu_standardized - mu_th) / mu_err) ** 2))
    chi2_per_dof = chi2_resid / max(n - 4 - n_bins, 1)

    result = {
        "alpha": float(alpha), "beta": float(beta), "gamma": float(gamma),
        "sigma_int": sigma_int,
        "n_sn": n, "chi2_per_dof": chi2_per_dof,
        "mu": mu_standardized.tolist(), "mu_err": mu_err.tolist(),
        "z": z.tolist(), "zhel": zhel.tolist(),
        "Om_fiducial": Om,
    }

    if mu_rel is not None:
        diff = mu_standardized - mu_rel
        diff -= np.median(diff)
        result["rms_vs_released_mu"] = float(np.std(diff))
        result["max_abs_diff_vs_released_mu"] = float(np.max(np.abs(diff)))

    return result


def _read_matrix_txt(path) -> np.ndarray:
    """SNANA/Pantheon+ covariance text: first line N, then N*N values (optionally gzipped)."""
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        n = int(fh.readline())
        vals = np.loadtxt(fh)
    return vals.reshape(n, n)


def _read_upper_tri_npz(path) -> np.ndarray:
    """DES-Dovekie STAT+SYS.npz: INVERSE covariance stored as upper triangle (see the DES likelihood)."""
    d = np.load(path)
    n = int(np.atleast_1d(d["nsn"])[0])
    inv = np.zeros((n, n))
    inv[np.triu_indices(n)] = d["cov"]
    il = np.tril_indices(n, -1)
    inv[il] = inv.T[il]
    return inv


def _read_table(path, sep=None) -> dict:
    """Minimal header+rows reader (no pandas so the kernel stays numpy/scipy/astropy-only)."""
    with open(path) as fh:
        lines = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
    cols = [c.strip() for c in lines[0].split(sep)]
    rows = [[v.strip() for v in ln.split(sep)] for ln in lines[1:]]
    out = {}
    for j, c in enumerate(cols):
        vals = [r[j] for r in rows]
        try:
            out[c] = np.array(vals, dtype=float)
        except ValueError:
            out[c] = np.array(vals, dtype=object)
    return out


def parse_dataset(dataset: str, hd_path, cov_path=None, cov_matrix=None) -> dict:
    """Parse one release into full-length arrays (no cuts).

    Returns dict(name, z, zhel, mu, mu_err, cov, is_calibrator, survey, cid, mu_kind).
    ``cov`` is the TOTAL (stat+sys) covariance of ``mu``. ``cov_matrix`` optionally passes the
    already-parsed text matrix (Pantheon+/DES-SN5YR) to skip the slow text read.
    """
    if dataset == "pantheonplus":
        d = _read_table(hd_path)
        cov = (_read_matrix_txt(cov_path) if cov_matrix is None else np.asarray(cov_matrix))
        return dict(z=d["zHD"], zhel=d["zHEL"],
                    # MU_SH0ES = m_b_corr + 19.253 (SH0ES M_B); the offset is marginalised anyway
                    mu=d["MU_SH0ES"], mu_err=np.sqrt(np.diag(cov)),
                    cov=cov, is_calibrator=d["IS_CALIBRATOR"].astype(int),
                    survey=d["IDSURVEY"].astype(int), cid=[str(c) for c in d["CID"]],
                    mu_kind="per-SN, MU_SH0ES (= m_b_corr + 19.253)")
    if dataset == "des5yr":
        d = _read_table(hd_path, sep=",")
        csys = (_read_matrix_txt(cov_path) if cov_matrix is None else np.asarray(cov_matrix))
        n = len(d["zHD"])
        if csys.shape[0] != n:
            raise ValueError(f"DES-SN5YR covariance is {csys.shape[0]}x{csys.shape[0]} but HD has {n} rows")
        err = d["MUERR_FINAL"]
        cov = csys + np.diag(err ** 2)  # STAT+SYS.txt.gz is systematics only (DES 5_COSMOLOGY likelihood)
        return dict(z=d["zHD"], zhel=d["zHEL"], mu=d["MU"],
                    mu_err=np.sqrt(np.diag(cov)), cov=cov, is_calibrator=np.zeros(n, int),
                    survey=d["IDSURVEY"].astype(int), cid=[str(c) for c in d["CID"]],
                    mu_kind="per-SN, MU (H0=70 convention)")
    if dataset == "des_dovekie":
        rows = []
        with open(hd_path) as fh:
            for line in fh:
                if line.startswith("SN:"):
                    rows.append(line.split()[1:])
        cid = [r[0] for r in rows]
        arr = np.array([[float(x) for x in r[1:6]] for r in rows])
        inv = _read_upper_tri_npz(cov_path)
        cov = np.linalg.inv(inv)
        return dict(z=arr[:, 1], zhel=arr[:, 2], mu=arr[:, 3], mu_err=np.sqrt(np.diag(cov)), cov=cov,
                    is_calibrator=np.zeros(len(cid), int), survey=arr[:, 0].astype(int), cid=cid,
                    mu_kind="per-SN, MU (H0=70 convention)")
    if dataset == "union3":
        from astropy.io import fits

        with fits.open(hd_path) as h:
            m = np.array(h[0].data, dtype=float)
        z = m[0, 1:]
        icov = m[1:, 1:]
        cov = np.linalg.inv(icov)
        n = len(z)
        return dict(z=z, zhel=z, mu=m[1:, 0], mu_err=np.sqrt(np.diag(cov)), cov=cov,
                    is_calibrator=np.zeros(n, int), survey=np.zeros(n, int),
                    cid=[f"node{i:02d}" for i in range(n)],
                    mu_kind="binned spline nodes (UNITY1.5 compressed), arbitrary offset")
    raise ValueError(f"unknown dataset {dataset!r}; choose from {sorted(SOURCES)}")


def select(sample: dict, zmin: float = 0.0, zmax: float = 10.0, calibrators: str = "include") -> dict:
    """Apply redshift cuts and the calibrator choice; covariance is SUBSET (never re-inverted inverse)."""
    z = sample["z"]
    m = (z > zmin) & (z <= zmax)
    cal = sample["is_calibrator"].astype(bool)
    if calibrators == "exclude":
        m &= ~cal
    elif calibrators == "only":
        m &= cal
    elif calibrators != "include":
        raise ValueError("calibrators must be include | exclude | only")
    out = {k: (np.asarray(v)[m] if k not in ("cov", "mu_kind", "cid") else v) for k, v in sample.items()}
    out["cid"] = [c for c, keep in zip(sample["cid"], m, strict=True) if keep]
    out["cov"] = sample["cov"][np.ix_(m, m)]
    return out


def _download(url: str) -> str:
    import os
    import tempfile
    import urllib.request

    name = url.rstrip("/").split("/")[-1].replace("%2B", "+").replace("%3D", "=")
    path = os.path.join(tempfile.gettempdir(), "sn-kernel-" + name)
    if not os.path.exists(path):
        req = urllib.request.Request(url, headers={"User-Agent": "sn-mcp-kernel/0.1"})
        with urllib.request.urlopen(req, timeout=300) as r, open(path + ".part", "wb") as fh:
            fh.write(r.read())
        os.replace(path + ".part", path)
    return path


def load_from_source(source: dict) -> dict:
    """source = {dataset, zmin, zmax, calibrators[, hd_url, cov_url]} -> selected sample (downloads)."""
    ds = source["dataset"]
    urls = dict(SOURCES[ds])
    hd = _download(source.get("hd_url") or urls["hd"])
    cov_url = source.get("cov_url") or urls.get("cov")
    cov = _download(cov_url) if cov_url else None
    s = parse_dataset(ds, hd, cov)
    return select(s, source.get("zmin", 0.0), source.get("zmax", 10.0), source.get("calibrators", "include"))


# -------------------------------------------------------------- the fit
DEFAULT_RANGES = {"Om": [0.05, 0.65], "w": [-2.0, -0.3], "w0": [-2.0, 0.0], "wa": [-6.0, 3.0]}


def _hpd_levels(post: np.ndarray, masses=(0.6827, 0.9545)) -> list[float]:
    """Density thresholds enclosing the given posterior masses (highest-density region)."""
    p = np.sort(post.ravel())[::-1]
    c = np.cumsum(p) / p.sum()
    return [float(p[min(np.searchsorted(c, m), len(p) - 1)]) for m in masses]


def _summ1d(x: np.ndarray, like: np.ndarray) -> dict:
    """Mean, sd and 68/95% equal-tail intervals of a gridded 1D marginal."""
    w = like / like.sum()
    mean = float((x * w).sum())
    sd = float(np.sqrt(((x - mean) ** 2 * w).sum()))
    cdf = np.cumsum(w)
    q = [float(np.interp(p, cdf, x)) for p in (0.02275, 0.15866, 0.5, 0.84134, 0.97725)]
    edge = bool(w[0] > 0.01 * w.max() or w[-1] > 0.01 * w.max())
    return {"mean": mean, "sd": sd, "median": q[2], "lo68": q[1], "hi68": q[3], "lo95": q[0], "hi95": q[4],
            "hits_prior_edge": edge}


def fit_sn_cosmology(model: str = "lcdm", data: dict | None = None, source: dict | None = None,
                     ranges: dict | None = None, n_grid: int = 0, refine: bool = True) -> dict:
    """Grid (+ optimizer) fit of flat LCDM / wCDM / w0waCDM to a SN Hubble diagram.

    data   = {"z": [...], "zhel": [...], "mu": [...], "cov": [[...]]}  (inline arrays), OR
    source = {"dataset": "pantheonplus"|"des5yr"|"union3"|"des_dovekie", "zmin", "zmax", "calibrators"}
             (the kernel downloads the public files itself).
    ranges = {param: [lo, hi]} flat priors (defaults DEFAULT_RANGES); n_grid points per axis
             (0 -> 401 / 121 / 41 for lcdm / wcdm / w0wacdm).

    The offset M (degenerate with H0) is marginalised analytically. Returns best fit,
    chi2_min, 1D marginal summaries, and gridded 1D/2D marginal posteriors + profile chi2
    (nested lists) with HPD density levels for the 68/95% contours.
    """
    if model not in MODELS:
        raise ValueError(f"model must be one of {list(MODELS)}")
    if data is None:
        if source is None:
            raise ValueError("pass data (arrays) or source (dataset spec)")
        s = load_from_source(source)
        z, zhel, mu, cov = s["z"], s["zhel"], s["mu"], s["cov"]
    else:
        z = np.asarray(data["z"], float)
        zhel = np.asarray(data.get("zhel", data["z"]), float)
        mu = np.asarray(data["mu"], float)
        cov = np.asarray(data["cov"], float)
    icov = np.linalg.inv(cov)
    icov = 0.5 * (icov + icov.T)
    rs = icov.sum(axis=1)
    cs = float(rs.sum())
    names = MODELS[model]
    rng = {**DEFAULT_RANGES, **(ranges or {})}
    n = n_grid or {"lcdm": 401, "wcdm": 121, "w0wacdm": 41}[model]
    axes = [np.linspace(rng[p][0], rng[p][1], n) for p in names]

    def chi2(vals) -> float:
        p = dict(zip(names, vals, strict=True))
        if not (0.0 < p["Om"] < 1.0):
            return np.inf
        th = _model_mu(z, zhel, model, p)
        if not np.all(np.isfinite(th)):
            return np.inf
        return marg_chi2(mu - th, icov, rs, cs)

    # batched evaluation: theory per point (cheap 1D integral), chi2 via BLAS-3 on residual blocks
    shape = tuple(len(a) for a in axes)
    pts = np.stack([m.ravel() for m in np.meshgrid(*axes, indexing="ij")], axis=1)
    flat = np.full(len(pts), np.inf)
    batch = 512
    for start in range(0, len(pts), batch):
        blk = pts[start:start + batch]
        th = np.empty((len(blk), len(z)))
        for j, vals in enumerate(blk):
            th[j] = _model_mu(z, zhel, model, dict(zip(names, vals, strict=True)))
        r = mu[None, :] - th
        a = np.einsum("ij,ij->i", r @ icov, r)
        b = r @ rs
        c2 = a - b * b / cs
        c2[~np.isfinite(c2)] = np.inf
        flat[start:start + batch] = c2
    grid = flat.reshape(shape)
    finite = np.isfinite(grid)
    c2min_grid = float(grid[finite].min())
    ibest = np.unravel_index(np.where(finite, grid, np.inf).argmin(), grid.shape)
    best = [float(axes[k][i]) for k, i in enumerate(ibest)]
    c2min = c2min_grid
    if refine:
        from scipy.optimize import minimize

        bounds = [(rng[p][0], rng[p][1]) for p in names]
        res = minimize(chi2, best, method="Nelder-Mead", bounds=bounds,
                       options={"xatol": 1e-5, "fatol": 1e-4})
        if res.success and res.fun <= c2min_grid + 1e-9:
            best, c2min = [float(v) for v in res.x], float(res.fun)
    like = np.where(finite, np.exp(-0.5 * (grid - c2min_grid)), 0.0)

    marg1d = {}
    for k, p in enumerate(names):
        other = tuple(j for j in range(len(names)) if j != k)
        l1 = like.sum(axis=other) if other else like
        prof = np.where(finite, grid, np.inf).min(axis=other) if other else grid
        marg1d[p] = {"x": axes[k].tolist(), "post": (l1 / l1.max()).tolist(), "chi2_profile": prof.tolist(),
                     **_summ1d(axes[k], l1)}
    marg2d = {}
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            other = tuple(j for j in range(len(names)) if j not in (a, b))
            l2 = like.sum(axis=other) if other else like
            prof = np.where(finite, grid, np.inf).min(axis=other) if other else grid
            l2 = l2 / l2.max()
            marg2d[f"{names[a]}-{names[b]}"] = {
                "x": axes[a].tolist(), "y": axes[b].tolist(), "post": l2.tolist(),
                "chi2_profile": np.where(np.isfinite(prof), prof, 1e30).tolist(),
                "hpd_levels_68_95": _hpd_levels(l2)}
    pbest = dict(zip(names, best, strict=True))
    offset = best_offset(mu - _model_mu(z, zhel, model, pbest), icov)
    return {"model": model, "params": names, "best_fit": pbest, "chi2_min": c2min, "n_sn": int(len(z)),
            "dof": int(len(z) - len(names) - 1), "best_offset_mag": offset,
            "ranges": {p: rng[p] for p in names}, "n_grid": n,
            "marginals_1d": marg1d, "marginals_2d": marg2d,
            "likelihood": "Gaussian, full covariance, offset M/H0 marginalised analytically (flat prior)"}
