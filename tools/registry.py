"""Single source of truth for the datasets this server exposes.

Keys are short dataset ids used by every tool's ``dataset`` argument.
"access" is "bundled" (file in data/), "download" (fetched on first use
into the data cache) or "remote-query" (live service). All URLs below were
verified by download during development (2026-11-04).
"""

from .sn.kernels import SOURCES

_TIER2_TOOLS = ["fetch_sn_sample", "bin_hubble_diagram", "plot_hubble_diagram", "fit_sn_cosmology",
                "plot_sn_contours", "compare_sn_samples", "compute_sn_chi2", "plot_sn_covariance"]
_TIER1_TOOLS = ["fetch_sn_salt2_params", "standardize_sn_salt2",
                "fetch_sn_lightcurve", "plot_sn_lightcurve"]

DATASETS: dict[str, dict] = {
    "pantheonplus": {
        "title": "Pantheon+ SN Ia compilation (1701 light curves of 1550 SNe, 18 surveys)",
        "products": {
            "tier1": ["per-SN SALT2 fit parameters (x1, c, mB, errors, covariances, HOST_LOGMASS) for 1701 LCs",
                      "photometric light curves for ~1220 SNe (SNANA text format, per-survey folders)"],
            "tier2": ["per-SN Hubble diagram (zHD, zHEL, m_b_corr, MU_SH0ES, IS_CALIBRATOR)",
                      "STAT+SYS covariance 1701x1701 (mu)"],
            "tier3": ["published MCMC chains (flat LCDM, flat wCDM, flat w0waCDM, SN-only and SN+SH0ES)"],
        },
        "access": "Hubble diagram + SALT2 params bundled (data/pantheonplus/Pantheon+SH0ES.dat, 579 kB); "
                  "covariance downloaded on first use (33 MB text, cached as .npy); "
                  "light curves downloaded per-SN on demand from PantheonPlusSH0ES/DataRelease/1_DATA/photometry",
        "urls": [SOURCES["pantheonplus"]["hd"], SOURCES["pantheonplus"]["cov"],
                 SOURCES["pantheonplus"]["salt2_params"],
                 "https://github.com/PantheonPlusSH0ES/DataRelease"],
        "citation": "Scolnic et al. 2022, ApJ 938, 113 (arXiv:2112.03863); "
                    "Brout et al. 2022, ApJ 938, 110 (arXiv:2202.04077)",
        "default_zmin": 0.01,
        "conventions": "mu = MU_SH0ES = m_b_corr + 19.253 (SH0ES M_B). Cosmology fits use zHD > 0.01 "
                       "(Brout+22 / cobaya PantheonPlus). 77 rows are Cepheid-host calibrators "
                       "(IS_CALIBRATOR=1; mostly z < 0.01). Covariance is the TOTAL stat+sys on mu.",
        "reference_values": {"flat LCDM Om (SN only, Brout+22)": "0.334 +- 0.018",
                             "flat wCDM (SN only, Brout+22)": "Om = 0.309 +0.063/-0.069, w = -0.90 +- 0.14",
                             "Tripp alpha (Brout+22)": "~0.15", "Tripp beta (Brout+22)": "~3.0",
                             "mass step gamma (Brout+22)": "~0.05 mag"},
        "tools": _TIER2_TOOLS + _TIER1_TOOLS,
    },
    "des5yr": {
        "title": "DES-SN5YR (Vincenzi+24 / DES Collaboration 2024): 1635 DES SNe + 194 low-z = 1829",
        "products": {
            "tier1": ["per-SN SALT3 fit parameters (x1, c, mB via x0, errors, covariances, HOST_LOGMASS) "
                      "from DES-SN5YR_HD+MetaData.csv",
                      "photometric light curves (DES5YR_SALT3_LCFIT.LCPLOT.gz, data+model, 7.4 MB)"],
            "tier2": ["per-SN Hubble diagram (zHD, zHEL, MU, MUERR_FINAL)",
                      "STAT+SYS systematic covariance 1829x1829 (stat added on the diagonal from MUERR_FINAL)"],
            "tier3": ["published polychord/emcee chains (flat LCDM, flat wCDM, flat w0waCDM)"],
        },
        "access": "Hubble diagram bundled (data/des5yr/DES-SN5YR_HD.csv, 92 kB); covariance downloaded "
                  "(16 MB gzip) from the tagged release 1.3; "
                  "SALT3 params downloaded on first use (DES-SN5YR_HD+MetaData.csv, 662 kB); "
                  "light curves from DES5YR_SALT3_LCFIT.LCPLOT.gz (7.4 MB, downloaded on first use)",
        "urls": [SOURCES["des5yr"]["hd"], SOURCES["des5yr"]["cov"],
                 SOURCES["des5yr"]["salt2_params"], SOURCES["des5yr"]["lcplot"],
                 "https://github.com/des-science/DES-SN5YR/releases/tag/1.3",
                 "https://doi.org/10.5281/zenodo.12720778"],
        "citation": "DES Collaboration 2024, ApJL 973, L14 (arXiv:2401.02929); Vincenzi et al. 2024 (arXiv:2401.02945)",
        "default_zmin": 0.0,
        "conventions": "mu = MU (H0 = 70 convention, offset marginalised). Total covariance = "
                       "STAT+SYS.txt.gz + diag(MUERR_FINAL^2), as in the DES 5_COSMOLOGY likelihood. "
                       "SALT3 used (not SALT2). Published: alpha=0.161, beta=3.118, gamma=0.038, M0_avg=-29.958.",
        "reference_values": {"flat LCDM Om (SN only, DES 2024)": "0.352 +- 0.017",
                             "flat wCDM (SN only, DES 2024)": "Om = 0.264 +0.074/-0.096, w = -0.80 +0.14/-0.16",
                             "Tripp alpha (Vincenzi+24)": "0.161 +- 0.002",
                             "Tripp beta (Vincenzi+24)": "3.118 +- 0.035",
                             "mass step gamma (Vincenzi+24)": "0.038 +- 0.008 mag"},
        "tools": _TIER2_TOOLS + _TIER1_TOOLS,
    },
    "union3": {
        "title": "Union3 + UNITY1.5 (Rubin+2023): 2087 SNe compressed to 22 spline-node distance moduli",
        "products": {
            "tier1": [],
            "tier2": ["binned distance moduli at 22 redshift nodes", "22x22 inverse covariance"],
            "tier3": [],
        },
        "access": "bundled (data/union3/mu_mat_union3_cosmo=2_mu.fits, 8.6 kB)",
        "urls": [SOURCES["union3"]["hd"], "https://github.com/rubind/union3_release",
                 "https://github.com/CobayaSampler/sn_data/tree/master/Union3"],
        "citation": "Rubin et al. 2023, arXiv:2311.12098 (Union3 & UNITY1.5)",
        "default_zmin": 0.0,
        "conventions": "FITS matrix: row 0 = z nodes, column 0 = mu, rest = INVERSE covariance "
                       "(rubind/union3_release README). Identical to cobaya's lcparam_full.txt + mag_covmat.txt. "
                       "The absolute offset is arbitrary; marginalised. There is NO public per-SN Union3 "
                       "Hubble diagram: the compressed nodes are the UNITY posterior. Calibrator selection "
                       "does not apply (no per-SN rows).",
        "reference_values": {"flat LCDM Om (SN only, Rubin+23)": "0.356 +0.028/-0.026"},
        "tools": _TIER2_TOOLS,
    },
    "des_dovekie": {
        "title": "DES-Dovekie (Popovic+2026) recalibrated DES-SN5YR: 1820 SNe",
        "products": {
            "tier1": [],
            "tier2": ["per-SN Hubble diagram (zHD, zHEL, MU, MUERR)", "STAT+SYS INVERSE covariance (upper-triangle npz)"],
            "tier3": [],
        },
        "access": "Hubble diagram bundled (data/des_dovekie/DES-Dovekie_HD.csv, 148 kB); inverse covariance "
                  "downloaded (6 MB npz)",
        "urls": [SOURCES["des_dovekie"]["hd"], SOURCES["des_dovekie"]["cov"]],
        "citation": "Popovic et al. 2026, arXiv:2511.07517",
        "default_zmin": 0.0,
        "conventions": "The npz stores the inverse of the TOTAL stat+sys covariance (DES-Dovekie-SN_Likelihood.py). "
                       "Optional extra sample; not one of the three headline datasets.",
        "reference_values": {},
        "tools": _TIER2_TOOLS,
    },
}
