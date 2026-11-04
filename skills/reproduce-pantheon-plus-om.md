---
name: reproduce-pantheon-plus-om
description: Reproduce the Pantheon+ SN-only flat-LCDM Om = 0.334 +- 0.018 (Brout+22 Table 3) and wCDM w = -0.90 +- 0.14, with the checks that matter - z > 0.01 cut, calibrators, full vs diagonal covariance, offset marginalisation
---

# Reproduce Pantheon+ Om (Brout et al. 2022)

Target: Brout+22 (ApJ 938, 110), SN-only, flat LCDM: **Om = 0.334 +- 0.018**;
flat wCDM: **Om = 0.309 +0.063/-0.069, w = -0.90 +- 0.14**. The likelihood is
the Pantheon+ STAT+SYS covariance with M marginalised (no SH0ES Cepheids).
Use one `output_dir`.

## 1. The sample
`fetch_sn_sample(output_dir, dataset="pantheonplus")` -> `n` = 1590 (zHD > 0.01,
the cobaya/Brout+22 cut; 1701 rows total). Check `z_range` ~ [0.0102, 2.26].
`n_calibrators` = 10: Cepheid hosts above z = 0.01 stay in by default, as in
the Hubble-flow sample of Brout+22 (they are ordinary SNe with Cepheid
distances that only matter when M is calibrated).

## 2. The fit
`fit_sn_cosmology(output_dir, dataset="pantheonplus", model="lcdm")`.
Expect `marginals.Om` = 0.332 +- 0.018, `chi2_min` = 1402.9 for 1590 SNe
(chi2/N = 0.88, as in Brout+22: the covariance is conservative). The 0.002
offset from 0.334 is within the quoted precision (Brout+22 used an MCMC with
CAMB distances; this is a 401-point grid with no radiation).

## 3. Robustness (each ~2 s)
- `calibrators="exclude"` -> 1580 SNe; Om = 0.333 +- 0.018 (shift ~0.001).
- `zmin=0.0` (all 1701 rows incl. the z < 0.01 calibrators, whose MU_SH0ES
  values are not Hubble-flow distances) -> Om = 0.363: a 1.6 sigma shift -
  this is why the cut exists. Never fit without it.
- `zmin=0.023` (classic Hubble-flow cut) -> 1371 SNe, Om = 0.329 +- 0.019.
- Diagonal-only covariance is NOT offered on purpose; `plot_sn_covariance`
  shows max off-diagonal correlation ~0.93 - systematics dominate.

## 4. wCDM
`fit_sn_cosmology(output_dir, dataset="pantheonplus", model="wcdm")` ->
w = -0.91 +- 0.15, Om = 0.29 +- 0.07 (mean of a banana-shaped posterior;
Brout+22 quote the mode 0.309). `plot_sn_contours(output_dir,
contour_csvs=[<*_2d_Om_w.csv>])` reproduces the SN-only contour of
Brout+22 Fig. 9 (without the CMB/BAO combination).

## 5. Spot checks
`compute_sn_chi2(dataset="pantheonplus", Om=0.315)` (Planck 2018) ->
delta chi2 vs the best fit = 0.85 (< 1 sigma), consistent with Brout+22 saying
Pantheon+ agrees with Planck LCDM.

## Report
Om mean +- sd vs 0.334 +- 0.018; chi2_min/N; the cut/calibrator variations;
w for wCDM; figure paths. State that offset M is marginalised (no H0 here).
