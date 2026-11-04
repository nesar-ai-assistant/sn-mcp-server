"""Readers for the SN releases (local side): resolve files via tools.common.fetch_file
(bundled data/ first, then the cache, then download) and parse with the kernel parsers
so local and dispatched runs see byte-identical data."""

from __future__ import annotations

import numpy as np

from ..common import BUNDLED_DATA, data_cache_dir, fetch_file, get_cached
from . import kernels

DATASET_NAMES = ("pantheonplus", "des5yr", "union3", "des_dovekie")

_LOCAL_NAMES = {
    "pantheonplus": ("Pantheon+SH0ES.dat", "Pantheon+SH0ES_STAT+SYS.cov"),
    "des5yr": ("DES-SN5YR_HD.csv", "STAT+SYS.txt.gz"),
    "union3": ("mu_mat_union3_cosmo=2_mu.fits", None),
    "des_dovekie": ("DES-Dovekie_HD.csv", "STAT+SYS.npz"),
}


def dataset_files(dataset: str) -> tuple:
    """Local paths (hd, cov|None), downloading on first use."""
    if dataset not in kernels.SOURCES:
        raise ValueError(f"Unknown dataset {dataset!r}. Available: {', '.join(DATASET_NAMES)}")
    src = kernels.SOURCES[dataset]
    hd_name, cov_name = _LOCAL_NAMES[dataset]
    hd = fetch_file(src["hd"], hd_name, subdir=dataset)
    cov = fetch_file(src["cov"], cov_name, subdir=dataset) if cov_name else None
    return hd, cov


def _load_full(dataset: str) -> dict:
    hd, cov = dataset_files(dataset)
    mat = None
    if cov is not None and str(cov).endswith((".cov", ".gz")):
        # text covariances are slow to parse (33 MB for Pantheon+): cache as .npy once
        npy = data_cache_dir() / dataset / (cov.name + ".npy")
        if npy.is_file():
            mat = np.load(npy)
        else:
            mat = kernels._read_matrix_txt(cov)
            np.save(npy, mat)
    return kernels.parse_dataset(dataset, hd, cov, cov_matrix=mat)


def load_full(dataset: str) -> dict:
    return get_cached(f"sn:{dataset}", lambda: _load_full(dataset))


def load_sample(dataset: str, zmin: float | None = None, zmax: float = 10.0,
                calibrators: str = "include") -> dict:
    """Selected sample; zmin=None uses the release's standard cosmology cut (Pantheon+ z > 0.01)."""
    from ..registry import DATASETS

    if zmin is None:
        zmin = DATASETS[dataset].get("default_zmin", 0.0)
    return kernels.select(load_full(dataset), zmin, zmax, calibrators)


# ================================================================= catalog-level helpers (per-SN fits, light curves)

def load_salt2_params(dataset: str) -> dict:
    """Load per-SN SALT2/SALT3 light-curve fit parameters for Pantheon+ or DES-SN5YR.

    Returns dict with arrays: cid, z, zhel, mB, mBERR, x1, x1ERR, c, cERR,
    COV_x1_c, COV_mB_x1, COV_mB_c, HOST_LOGMASS, HOST_LOGMASS_ERR,
    biasCor_mB, mu_released, survey.
    """
    import csv as csv_mod

    src = kernels.SOURCES[dataset]
    LN10 = np.log(10)

    if dataset == "pantheonplus":
        hd, _ = dataset_files(dataset)  # bundled Pantheon+SH0ES.dat already has all SALT2 params
        d = kernels._read_table(hd)
        # COV_x1_x0 and COV_c_x0 need converting from x0-space to mB-space
        f = -2.5 / (LN10 * d["x0"])
        return {
            "cid": [str(c) for c in d["CID"]], "z": d["zHD"], "zhel": d["zHEL"],
            "mB": d["mB"], "mBERR": d["mBERR"], "x1": d["x1"], "x1ERR": d["x1ERR"],
            "c": d["c"], "cERR": d["cERR"], "COV_x1_c": d["COV_x1_c"],
            "COV_mB_x1": f * d["COV_x1_x0"], "COV_mB_c": f * d["COV_c_x0"],
            "HOST_LOGMASS": d["HOST_LOGMASS"], "HOST_LOGMASS_ERR": d["HOST_LOGMASS_ERR"],
            "biasCor_mB": d["biasCor_m_b"],
            "mu_released": d["MU_SH0ES"],
            "survey": d["IDSURVEY"].astype(int),
            "is_calibrator": d["IS_CALIBRATOR"].astype(int),
            "VPECERR": d["VPECERR"],
            "salt_version": "SALT2",
        }

    if dataset == "des5yr":
        meta_url = src["salt2_params"]
        meta_path = fetch_file(meta_url, "DES-SN5YR_HD+MetaData.csv", subdir="des5yr")
        with open(meta_path) as fh:
            reader = csv_mod.DictReader(fh)
            rows = list(reader)
        def col(name):
            return np.array([float(r[name]) for r in rows])
        x0 = col("x0")
        f = -2.5 / (LN10 * x0)
        return {
            "cid": [r["CID"] for r in rows], "z": col("zHD"), "zhel": col("zHEL"),
            "mB": -2.5 * np.log10(x0), "mBERR": 2.5 / LN10 * col("x0ERR") / x0,
            "x1": col("x1"), "x1ERR": col("x1ERR"),
            "c": col("c"), "cERR": col("cERR"), "COV_x1_c": col("COV_x1_c"),
            "COV_mB_x1": f * col("COV_x1_x0"), "COV_mB_c": f * col("COV_c_x0"),
            "HOST_LOGMASS": col("HOST_LOGMASS"), "HOST_LOGMASS_ERR": col("HOST_LOGMASS_ERR"),
            "biasCor_mB": col("biasCor_mu"),
            "mu_released": col("MU"),
            "survey": col("IDSURVEY").astype(int),
            "is_calibrator": np.zeros(len(rows), int),
            "VPECERR": col("VPECERR"),
            "salt_version": "SALT3",
        }

    raise ValueError(f"SALT2/3 params not available for {dataset!r}; use 'pantheonplus' or 'des5yr'.")


def load_lc_index() -> list[dict]:
    """Load the Pantheon+ CID → light-curve file index (bundled CSV)."""
    import csv as csv_mod
    idx_path = BUNDLED_DATA / "pantheonplus" / "lc_index.csv"
    with open(idx_path) as fh:
        return list(csv_mod.DictReader(fh))


def fetch_pantheonplus_lc(cid: str, folder: str | None = None) -> list[dict]:
    """Download and parse Pantheon+ SNANA light-curve file(s) for one SN.

    Returns a list of parsed LC dicts (one per survey file matching the CID).
    """
    idx = load_lc_index()
    matches = [r for r in idx if r["cid"] == cid]
    if folder:
        matches = [r for r in matches if r["folder"] == folder]
    if not matches:
        raise ValueError(f"No light-curve files found for CID={cid!r} in the Pantheon+ index"
                         + (f" (folder={folder!r})" if folder else "")
                         + ". Use fetch_sn_salt2_params to see available CIDs.")
    base = kernels.SOURCES["pantheonplus"]["lc_base"]
    results = []
    for row in matches:
        import urllib.parse
        url = base + urllib.parse.quote(row["path"])
        lc_path = fetch_file(url, row["path"].replace("/", "_"), subdir="pantheonplus_lc")
        text = lc_path.read_text(errors="replace")
        parsed = kernels.parse_snana_lc(text)
        parsed["source_file"] = row["path"]
        parsed["folder"] = row["folder"]
        results.append(parsed)
    return results


def fetch_des5yr_lc(cid: str) -> dict | None:
    """Extract one SN's light curve from DES5YR_SALT3_LCFIT.LCPLOT.gz (data + SALT3 model).

    Returns dict(snid, mjd, band, flux, flux_err, phase, data_model, chi2) or None.
    """
    import gzip

    lcplot_url = kernels.SOURCES["des5yr"]["lcplot"]
    lcplot_path = fetch_file(lcplot_url, "DES5YR_SALT3_LCFIT.LCPLOT.gz", subdir="des5yr")
    cid_str = str(cid).strip()
    rows = []
    with gzip.open(lcplot_path, "rt") as fh:
        for line in fh:
            if line.startswith(cid_str + " ") or line.startswith(cid_str + "\t"):
                rows.append(line.split())
    if not rows:
        return None
    arr = np.array([[float(x) for x in r[1:6]] for r in rows])
    bands = [r[6] for r in rows]
    chi2 = np.array([float(r[7]) for r in rows]) if len(rows[0]) > 7 else None
    dm = arr[:, 4].astype(int)
    return {
        "snid": cid_str, "mjd": arr[:, 0], "phase": arr[:, 1],
        "flux": arr[:, 2], "flux_err": arr[:, 3], "data_model": dm,
        "band": np.array(bands), "chi2": chi2,
    }
