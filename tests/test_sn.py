"""Science tests for the SN Ia family.

Offline tests use only bundled files (Union3 nodes, Hubble-diagram tables) or
synthetic data. Network tests (covariances, ~6-33 MB, cached after first run)
reproduce the published SN-only flat-LCDM Om values.
"""

import json
import os

import numpy as np
import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

from tools.sn import (  # noqa: E402
    DISPATCH_KERNELS,
    bin_hubble_diagram,
    compare_sn_samples,
    compute_distance_modulus,
    compute_sn_chi2,
    fetch_sn_sample,
    fetch_sn_lightcurve,
    fetch_sn_salt2_params,
    fit_sn_cosmology,
    plot_hubble_diagram,
    plot_sn_contours,
    plot_sn_lightcurve,
    standardize_sn_salt2,
)
from tools.sn import kernels  # noqa: E402


# ------------------------------------------------------------------ theory
def test_distance_modulus_matches_astropy():
    from astropy.cosmology import Flatw0waCDM, LambdaCDM

    z = np.array([0.01, 0.1, 0.5, 1.0, 2.0])
    ref = Flatw0waCDM(H0=70, Om0=0.3, w0=-0.9, wa=-0.5, Tcmb0=0).distmod(z).value
    mu = kernels.distance_modulus(z, om=0.3, w0=-0.9, wa=-0.5, h0=70)
    assert np.max(abs(mu - ref)) < 2e-4
    ref_k = LambdaCDM(H0=70, Om0=0.3, Ode0=0.6, Tcmb0=0).distmod(z).value
    assert np.max(abs(kernels.distance_modulus(z, om=0.3, ok=0.1) - ref_k)) < 2e-4


def test_marginalised_chi2_is_offset_invariant():
    rng = np.random.default_rng(1)
    n = 40
    a = rng.normal(size=(n, n)) * 0.02
    cov = a @ a.T + np.diag(np.full(n, 0.01))
    icov = np.linalg.inv(cov)
    r = rng.normal(size=n) * 0.1
    assert kernels.marg_chi2(r, icov) == pytest.approx(kernels.marg_chi2(r + 3.7, icov), rel=1e-10)
    off = kernels.best_offset(r, icov)
    # marginalised chi2 equals the profiled chi2 for a flat prior on the offset
    assert kernels.marg_chi2(r, icov) == pytest.approx(float((r - off) @ icov @ (r - off)), rel=1e-9)


def test_kernel_recovers_injected_cosmology():
    rng = np.random.default_rng(2)
    z = np.sort(rng.uniform(0.02, 1.5, 400))
    cov = np.diag(np.full(len(z), 0.02 ** 2))
    mu = kernels.distance_modulus(z, om=0.3, w0=-1.0) + 0.17 + rng.normal(0, 0.02, len(z))
    res = kernels.fit_sn_cosmology("lcdm", data={"z": z.tolist(), "mu": mu.tolist(), "cov": cov.tolist()})
    om = res["marginals_1d"]["Om"]
    assert abs(om["mean"] - 0.3) < 3 * om["sd"]
    assert res["best_offset_mag"] == pytest.approx(0.17, abs=0.02)
    json.dumps(res)  # JSON-safe for dispatch


def test_dispatch_manifest():
    k = DISPATCH_KERNELS["fit_sn_cosmology"]
    assert k["function"] == "sn.kernels.fit_sn_cosmology"
    assert set(k["pip_deps"]) <= {"numpy", "scipy", "astropy"}
    src = open(kernels.__file__).read()
    assert "import pandas" not in src and "from ..common" not in src  # kernel stays standalone


def test_compute_distance_modulus_csv(tmp_path):
    r = compute_distance_modulus(str(tmp_path), Om=0.3, z_min=0.01, z_max=2.0, n_points=50)
    assert r.files[0].endswith(".csv")
    d = np.loadtxt(r.files[0], delimiter=",", skiprows=1)
    assert d.shape == (50, 4) and np.allclose(d[:, 3], 0.0, atol=1e-9)
    assert r.metadata["mu_at_z1"] == pytest.approx(44.10, abs=0.02)


# ---------------------------------------------- bundled data, offline
def test_union3_offline_fit(tmp_path):
    """Rubin+23 Union3/UNITY1.5 SN-only flat LCDM: Om = 0.356 +0.028 -0.026."""
    r = fit_sn_cosmology(str(tmp_path), dataset="union3", model="lcdm")
    m = r.metadata["marginals"]["Om"]
    assert m["mean"] == pytest.approx(0.356, abs=0.01)
    assert m["sd"] == pytest.approx(0.027, abs=0.004)
    assert r.metadata["n_sn"] == 22 and r.metadata["computed_on"] == "local"
    s = json.loads(open(r.metadata["summary_json"]).read())
    assert s["best_fit"]["Om"] == pytest.approx(m["mean"], abs=0.01)


def test_union3_plots_and_bins(tmp_path):
    r = plot_hubble_diagram(str(tmp_path), dataset="union3", n_bins=0)
    assert r.files[0].endswith(".png")
    f = fetch_sn_sample(str(tmp_path), dataset="union3")
    assert f.metadata["n"] == 22 and len(f.files) == 2
    assert np.load(f.files[1]).shape == (22, 22)


def test_bad_inputs(tmp_path):
    with pytest.raises(ValueError):
        fetch_sn_sample(str(tmp_path), dataset="union3", zmin=2.9, zmax=3.0)
    with pytest.raises(ValueError):
        fit_sn_cosmology(str(tmp_path), dataset="union3", Om_range=[0.5, 0.1])
    with pytest.raises(ValueError):
        compute_distance_modulus(str(tmp_path), z_min=1.0, z_max=0.5)


# ---------------------------------------------- network (covariances)
@pytest.mark.network
def test_pantheonplus_om(tmp_path):
    """Brout+22: Pantheon+ SN-only flat LCDM Om = 0.334 +- 0.018 (z > 0.01)."""
    r = fit_sn_cosmology(str(tmp_path), dataset="pantheonplus", model="lcdm")
    m = r.metadata["marginals"]["Om"]
    assert r.metadata["n_sn"] == 1590
    assert m["mean"] == pytest.approx(0.334, abs=0.006)
    assert m["sd"] == pytest.approx(0.018, abs=0.002)


@pytest.mark.network
def test_des5yr_om(tmp_path):
    """DES Collaboration 2024: DES-SN5YR SN-only flat LCDM Om = 0.352 +- 0.017."""
    r = fit_sn_cosmology(str(tmp_path), dataset="des5yr", model="lcdm")
    m = r.metadata["marginals"]["Om"]
    assert r.metadata["n_sn"] == 1829
    assert m["mean"] == pytest.approx(0.352, abs=0.006)
    assert m["sd"] == pytest.approx(0.017, abs=0.002)


@pytest.mark.network
def test_pantheonplus_calibrators_and_bins(tmp_path):
    inc = fetch_sn_sample(str(tmp_path), dataset="pantheonplus", zmin=0.0, calibrators="include")
    exc = fetch_sn_sample(str(tmp_path), dataset="pantheonplus", zmin=0.0, calibrators="exclude")
    only = fetch_sn_sample(str(tmp_path), dataset="pantheonplus", zmin=0.0, calibrators="only")
    assert inc.metadata["n"] == 1701 and only.metadata["n"] == 77
    assert exc.metadata["n"] + only.metadata["n"] == 1701
    b = bin_hubble_diagram(str(tmp_path), dataset="pantheonplus", n_bins=12)
    cb = np.load(b.files[1])
    assert cb.shape[0] == b.metadata["n_bins"] and np.all(np.linalg.eigvalsh(cb) > 0)
    assert compute_sn_chi2("pantheonplus", Om=0.334).metadata["chi2_per_n"] < 1.0


@pytest.mark.network
def test_wcdm_contours_and_compare(tmp_path):
    """Brout+22 SN-only flat wCDM: w = -0.90 +- 0.14; DES 2024: w = -0.80 +0.14/-0.16."""
    p = fit_sn_cosmology(str(tmp_path), dataset="pantheonplus", model="wcdm", n_grid=61)
    w = p.metadata["marginals"]["w"]
    assert w["mean"] == pytest.approx(-0.90, abs=0.05) and w["sd"] == pytest.approx(0.15, abs=0.03)
    csvs = [f for f in p.files if "_2d_Om_w" in f]
    assert csvs and plot_sn_contours(str(tmp_path), csvs).files[0].endswith(".png")
    c = compare_sn_samples(str(tmp_path), datasets=["pantheonplus", "des5yr", "union3"], mode="residuals")
    assert set(c.metadata["per_dataset"]) == {"pantheonplus", "des5yr", "union3"}
    t = compare_sn_samples(str(tmp_path), datasets=["pantheonplus", "des5yr"], mode="contours", model="lcdm")
    assert t.metadata["tension"]["pantheonplus vs des5yr"]["sigma_equiv"] < 2.0


@pytest.mark.network
def test_dispatch_branch_with_source(tmp_path, monkeypatch):
    """Heavy-tool remote branch: the kernel receives a JSON `source` (URLs, no local paths)."""
    import mcp_server.dispatch as d

    calls = {}

    def fake_run_kernel(function, args, pip_deps=None, duration=None, **kw):
        calls.update(function=function, args=args)
        json.dumps(args)
        assert function == "sn.kernels.fit_sn_cosmology"
        return {"result": kernels.fit_sn_cosmology(**args), "host": "fake-node"}

    monkeypatch.setattr(d, "remote_site", lambda: "perlmutter")
    monkeypatch.setattr(d, "run_kernel", fake_run_kernel)
    r = fit_sn_cosmology(str(tmp_path), dataset="union3", model="lcdm")
    assert calls["args"]["source"]["dataset"] == "union3"
    assert r.metadata["computed_on"] == "fake-node"
    assert r.metadata["marginals"]["Om"]["mean"] == pytest.approx(0.356, abs=0.01)


# ------------------------------------------------------------------ tier-1: SALT2 params
def test_fetch_salt2_params_pantheonplus(tmp_path):
    """fetch_sn_salt2_params returns a CSV with the expected SALT2 columns (offline, bundled data)."""
    r = fetch_sn_salt2_params(str(tmp_path), dataset="pantheonplus", zmin=0.01, zmax=0.05)
    assert r.status == "success"
    assert r.metadata["salt_version"] == "SALT2"
    assert r.metadata["n_sn"] > 10
    csv_path = r.files[0]
    import csv as csv_mod
    with open(csv_path) as fh:
        reader = csv_mod.DictReader(fh)
        rows = list(reader)
    assert len(rows) == r.metadata["n_sn"]
    # Check required columns
    for col in ("cid", "z", "mB", "x1", "c", "HOST_LOGMASS"):
        assert col in rows[0], f"missing column {col}"
    # Check that z values are in range
    for row in rows:
        z = float(row["z"])
        assert 0.01 < z <= 0.05


def test_fetch_salt2_params_has_tier(tmp_path):
    assert fetch_sn_salt2_params.data_tier == 1


# ------------------------------------------------------------------ tier-1: standardize
def test_standardize_salt2_pantheonplus(tmp_path):
    """standardize_sn_salt2 recovers alpha~0.15, beta~3.0 from Pantheon+ (offline, bundled)."""
    r = standardize_sn_salt2(str(tmp_path), dataset="pantheonplus", zmin=0.01, zmax=2.3)
    assert r.status == "success"
    alpha = r.metadata["alpha"]
    beta = r.metadata["beta"]
    gamma = r.metadata["gamma"]
    assert 0.10 < alpha < 0.20, f"alpha={alpha} outside expected range"
    assert 2.5 < beta < 3.5, f"beta={beta} outside expected range"
    assert -0.02 < gamma < 0.12, f"gamma={gamma} outside expected range"
    assert r.metadata["n_sn"] > 1500
    # rms_vs_released should be modest (diagonal approx without full systematic covariance: < 60 mmag)
    assert r.metadata["rms_vs_released_mu"] < 0.06


def test_standardize_has_tier(tmp_path):
    assert standardize_sn_salt2.data_tier == 1


# ------------------------------------------------------------------ tier-1: light curves
def test_lc_index_loads():
    """The bundled lc_index.csv is readable and has entries."""
    from tools.sn.io import load_lc_index
    idx = load_lc_index()
    assert len(idx) > 1000
    assert all("cid" in r and "folder" in r and "path" in r for r in idx[:5])


@pytest.mark.skipif(os.environ.get("SN_SKIP_NETWORK"), reason="network tests skipped")
def test_fetch_lc_pantheonplus(tmp_path):
    """fetch_sn_lightcurve downloads a known SN from Pantheon+ (network)."""
    r = fetch_sn_lightcurve(str(tmp_path), cid="2011fe", dataset="pantheonplus")
    assert r.status == "success"
    assert r.metadata["n_obs"] > 50  # SN 2011fe is well-observed


@pytest.mark.skipif(os.environ.get("SN_SKIP_NETWORK"), reason="network tests skipped")
def test_plot_lc_pantheonplus(tmp_path):
    """plot_sn_lightcurve produces a PNG for a known SN (network)."""
    r = plot_sn_lightcurve(str(tmp_path), cid="2011fe", dataset="pantheonplus")
    assert r.status == "success"
    assert any(str(f).endswith(".png") for f in r.files)


def test_lc_tools_have_tier():
    assert fetch_sn_lightcurve.data_tier == 1
    assert plot_sn_lightcurve.data_tier == 1


# ------------------------------------------------------------------ kernel: parse SNANA LC
def test_parse_snana_lc():
    """parse_snana_lc parses a minimal SNANA-format light curve."""
    text = """\
SURVEY:  CfA3
SNID: 2005el
FILTERS: UBVRri
RA: 18.8333
DECL: 5.4167
REDSHIFT_HELIO: 0.0149 +- 0.0001
MWEBV: 0.034
SEARCH_PEAKMJD: 53634.0

VARLIST: MJD FLT FLUXCAL FLUXCALERR MAG MAGERR
OBS: 53623.0 B 100.0 10.0 20.0 0.11
OBS: 53624.0 V 200.0 15.0 19.0 0.08
OBS: 53625.0 r 300.0 20.0 18.5 0.07
END:
"""
    lc = kernels.parse_snana_lc(text)
    assert lc["snid"] == "2005el"
    assert lc["survey"] == "CfA3"
    assert lc["z_hel"] == pytest.approx(0.0149, abs=0.001)
    assert len(lc["mjd"]) == 3
    assert lc["flux"][0] == pytest.approx(100.0)
    assert lc["band"][1] == "V"


# ------------------------------------------------------------------ dispatch manifest includes standardize
def test_standardize_in_dispatch_manifest():
    assert "standardize_sn_salt2" in DISPATCH_KERNELS
    d = DISPATCH_KERNELS["standardize_sn_salt2"]
    assert "function" in d
    assert "numpy" in d["pip_deps"]
