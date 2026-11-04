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
    fit_sn_cosmology,
    plot_hubble_diagram,
    plot_sn_contours,
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
