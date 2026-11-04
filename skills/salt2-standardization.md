---
name: salt2-standardization
description: Build a Hubble diagram from scratch using SALT2 light-curve fit parameters (tier-1 to tier-2) - Tripp formula with alpha, beta, gamma, sigma_int fitted by maximum likelihood
---

# SALT2/SALT3 Tripp Standardisation (tier 1 → tier 2)

Goal: turn per-SN SALT2 fit parameters (x1, c, mB, errors, covariances,
host mass) into a standardised Hubble diagram mu(z), fitting the nuisance
parameters alpha (stretch), beta (colour), gamma (host-mass step) and
sigma_int (intrinsic scatter) by maximum likelihood.

## Background

The Tripp formula (Tripp 1998):
  mu = mB - M + alpha*x1 - beta*c ± gamma/2 - biasCor

where the ± is +gamma/2 for HOST_LOGMASS > 10, -gamma/2 otherwise.
M is the absolute magnitude (degenerate with H0, marginalised per z-bin).

Literature values (Brout+22 / Vincenzi+24):
- alpha ~ 0.15 (Pantheon+), 0.161 (DES-SN5YR)
- beta ~ 3.0 (Pantheon+), 3.118 (DES-SN5YR)
- gamma ~ 0.05 (Pantheon+), 0.038 (DES-SN5YR)

## Steps

### 1. Get the SALT2 parameters
```
fetch_sn_salt2_params(output_dir, dataset="pantheonplus", zmin=0.01)
```
or `dataset="des5yr"`. The CSV includes mB, x1, c, their errors and
cross-covariances (COV_x1_c, COV_mB_x1, COV_mB_c), HOST_LOGMASS,
biasCor_mB, and the released mu for comparison.

### 2. Standardise
```
standardize_sn_salt2(output_dir, dataset="pantheonplus", zmin=0.01)
```
Returns: alpha, beta, gamma, sigma_int, mu, mu_err, chi2_per_dof,
rms_vs_released_mu.

### 3. Validate
- alpha, beta, gamma should match published values within ~0.01-0.02
- rms_vs_released_mu should be < 50 mmag (limited by diagonal approximation
  vs full BEAMS+systematic covariance in the official analysis)
- chi2_per_dof should be ~1

### 4. Use the result
The standardised mu(z) is a tier-2 Hubble diagram. Feed it to
fit_sn_cosmology or compare it with the released HD from fetch_sn_sample.

## Pitfalls
- Pantheon+ uses mB directly; DES-SN5YR stores x0 and we convert mB = -2.5 log10(x0)
- The covariances COV_mB_x1, COV_mB_c are derived from x0-space (COV_x1_x0, COV_c_x0)
  using the chain rule f = -2.5/(ln10 * x0)
- This tool uses a DIAGONAL per-SN variance, not the full N×N systematic covariance.
  The released HD incorporates the full covariance + BEAMS probability weighting, so
  the two will differ by ~30-50 mmag rms. This is expected and documented.
- The tool does NOT use any tier-3 inputs (published chains or posteriors).
