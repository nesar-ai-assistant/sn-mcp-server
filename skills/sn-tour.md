---
name: sn-tour
description: Five-minute all-local tour of the SN Ia server - list datasets, fetch Pantheon+, plot its Hubble diagram, fit flat LCDM to Pantheon+/DES-SN5YR/Union3, overlay residuals of the three samples, and draw a wCDM contour
---

# Supernova cosmology tour (all local, < 5 min)

Goal: one call per tool family so a new user sees every capability and the
headline numbers. Everything runs locally (first run downloads ~55 MB of
covariances into `~/.cache/cosmo-mcp/sn-mcp-server`; later runs are offline).
Use one `output_dir`, e.g. `output/tour`.

## 1. Datasets
`list_sn_datasets()` then `describe_sn_dataset("pantheonplus")`. Note the
conventions block: Pantheon+ uses z_HD > 0.01 for cosmology; Union3 is a
22-node compressed Hubble diagram, not per-SN.

## 2. Fetch
`fetch_sn_sample(output_dir, dataset="pantheonplus")` -> 1590 SNe (z > 0.01),
CSV + 1590x1590 covariance .npy. `calibrators="only"` gives the 77 Cepheid hosts
(with `zmin=0`).

## 3. Theory curve
`compute_distance_modulus(output_dir, Om=0.3)` -> `mu_at_z1` = 44.10 mag (H0 = 70).

## 4. Hubble diagram
`plot_hubble_diagram(output_dir, dataset="pantheonplus")` (cf. Brout+22 Fig. 4).
`bin_hubble_diagram(output_dir, dataset="pantheonplus", n_bins=20)` -> binned
residuals + correlated bin covariance.

## 5. Fit flat LCDM (three samples, ~3 s each)
`fit_sn_cosmology(output_dir, dataset=<ds>, model="lcdm")`. Expect:

| dataset | N | Om (this server) | published |
|---|---|---|---|
| pantheonplus | 1590 | 0.332 +- 0.018 | 0.334 +- 0.018 (Brout+22) |
| des5yr | 1829 | 0.353 +- 0.017 | 0.352 +- 0.017 (DES 2024) |
| union3 | 22 | 0.358 +- 0.027 | 0.356 +0.028/-0.026 (Rubin+23) |

## 6. wCDM contour
`fit_sn_cosmology(output_dir, dataset="pantheonplus", model="wcdm")` (~10 s) ->
w = -0.91 +- 0.15 (Brout+22: -0.90 +- 0.14). Then
`plot_sn_contours(output_dir, contour_csvs=[<the *_2d_Om_w.csv file>])`.

## 7. Compare samples
`compare_sn_samples(output_dir, mode="residuals")` overlays the three binned
residual Hubble diagrams vs flat LCDM Om = 0.3.

## Report
The table of Om values, the wCDM w, the four PNG paths, and one sentence on
whether the samples agree (they do in LCDM: pairwise Om tensions 0.2-0.8 sigma).
