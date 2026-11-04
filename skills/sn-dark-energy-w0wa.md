---
name: sn-dark-energy-w0wa
description: Fit wCDM and w0waCDM to Pantheon+, DES-SN5YR and Union3 separately, overlay their contours, and quantify how far each prefers evolving dark energy and how much the samples disagree (residual shapes, Om shifts, pairwise tension) - with the caveats on shared SNe
---

# SN dark energy: w and (w0, wa) across samples

Context: DESI DR1/DR2 BAO + CMB prefer w0 > -1, wa < 0, with significance
depending on WHICH SN sample is added (Pantheon+ weakest, DES-SN5YR
strongest, Union3 between). This recipe shows the SN-side origin of that
spread. SN alone cannot pin (w0, wa): expect broad, prior-edge-limited
contours along the w0-wa degeneracy. Use one `output_dir`.

## 1. LCDM baseline (fast)
`compare_sn_samples(output_dir, datasets=["pantheonplus","des5yr","union3"],
mode="contours", model="lcdm")` -> Om 0.332 +- 0.018 / 0.353 +- 0.017 /
0.358 +- 0.027; pairwise tensions < 1 sigma. The Om ordering is the LCDM
fingerprint of what drives w0wa.

## 2. Residual shapes
`compare_sn_samples(output_dir, mode="residuals", fid_Om=0.3)`: DES-SN5YR and
Union3 sit low at z > 0.5 relative to their low-z anchors more than Pantheon+
does (see `mean_resid_z_gt_0p5`); the low-z vs high-z offset is what pushes
the fits to higher Om or w0 > -1. Note low-z SNe (CfA/CSP/Foundation) are
shared and DES-SN5YR's low-z set differs from Pantheon+'s (calibration
offset of ~0.04 mag was debated: Efstathiou 2024, Vincenzi+2025).

## 3. wCDM per sample (~10 s each)
`fit_sn_cosmology(output_dir, dataset=<ds>, model="wcdm")`. Expected (SN only):

| dataset | w | published |
|---|---|---|
| pantheonplus | -0.91 +- 0.15 | -0.90 +- 0.14 (Brout+22) |
| des5yr | -0.81 +- 0.15 | -0.80 +0.14/-0.16 (DES 2024) |
| union3 | -0.77 +- 0.17 (Om at prior edge) | (Rubin+23 quote w with CMB/BAO) |

Overlay: `plot_sn_contours(output_dir, contour_csvs=[3 x *_2d_Om_w.csv])`.

## 4. w0waCDM per sample (heavy: 30-50 s each locally at the default 41^3 grid)
`fit_sn_cosmology(output_dir, dataset=<ds>, model="w0wacdm")` (defaults:
Om [0.05, 0.65], w0 [-2, 0], wa [-6, 3]). Or one call:
`compare_sn_samples(output_dir, mode="contours", model="w0wacdm")`.
For publication grids (n_grid >= 61) dispatch to HPC - see
`hpc-dispatch-handoff`. Read for each sample: the w0, wa means and whether
`hits_prior_edge` is true (it will be for wa - report that the constraint is
prior-limited, never quote it as a measurement). Without a CMB/BAO Om prior,
SN-only w0wa is dominated by the Om-w0-wa degeneracy; add an informative
`Om_range` (e.g. [0.28, 0.34], a crude stand-in for CMB+BAO) to show how the
samples then separate in (w0, wa).

## 5. Quantify
- Distance from LCDM per sample: lcdm `chi2_min` minus w0wacdm `chi2_min`
  (2 extra parameters). Measured with default priors: Pantheon+ 0.5,
  Union3 3.4, DES-SN5YR 7.8 (best fit at the wa = -6 prior edge, Om ~ 0.47;
  cf. DES 2024 SN-only w0wa Om = 0.495, w0 = -0.36, wa = -8.8). The DES
  improvement is a degeneracy-direction fit, not a detection on its own.
  Expected w0wacdm means (default priors): Pantheon+ w0 = -0.92 +- 0.15,
  wa = -1.0 +- 1.7; DES-SN5YR w0 = -0.70 +- 0.15, wa = -4.1 +- 1.6;
  Union3 w0 = -0.65 +- 0.20, wa = -2.9 +- 1.9 (all wa prior-edge flagged).
- Between samples: the `tension` block of compare_sn_samples (Gaussian
  approximation, C1 + C2). It ignores shared SNe, so it UNDER-states the
  significance of a genuine difference; quote as indicative.

## Report
Table of Om (LCDM), w (wCDM), w0/wa means with prior-edge flags per sample;
pairwise tensions; delta chi2 vs LCDM; the contour PNGs; a clear statement
that SN-only data do not detect evolving dark energy and that the sample
differences sit mainly in the low-z/high-z relative calibration.
