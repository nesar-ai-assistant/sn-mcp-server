# Bundled data (offline core, ~1 MB total)

Each file is byte-identical to the upstream URL (md5-checked 2026-11-04). The large
covariances are NOT bundled: `tools.common.fetch_file` downloads them on first use into
`~/.cache/cosmo-mcp/sn-mcp-server/<dataset>/` (override with `$COSMO_MCP_DATA`).

## Data products
| kind | what | examples in this server |
|---|---|---|
| catalog | object-level data | per-SN SALT2/SALT3 light-curve fit parameters (x1, c, mB, HOST_LOGMASS), photometric light curves |
| summary statistic | measurements + covariance | Hubble diagrams mu(z) with full stat+sys covariance |
| published result | published inference | posterior chains (for comparison only, not used as analysis input) |

## Bundled files

| file | kind | size | upstream URL | citation |
|---|---|---|---|---|
| pantheonplus/Pantheon+SH0ES.dat | catalog + summary statistic | 579 kB | PantheonPlusSH0ES/DataRelease .../Pantheon%2BSH0ES.dat | Scolnic+2022 ApJ 938 113; Brout+2022 ApJ 938 110 |
| pantheonplus/lc_index.csv | catalog | 40 kB | built from PantheonPlusSH0ES/DataRelease/1_DATA/photometry tree | — |
| des5yr/DES-SN5YR_HD.csv | summary statistic | 92 kB | des-science/DES-SN5YR/1.3/4_DISTANCES_COVMAT/DES-SN5YR_HD.csv | DES 2024 ApJL 973 L14; Vincenzi+2024 |
| union3/mu_mat_union3_cosmo=2_mu.fits | summary statistic | 8.6 kB | rubind/union3_release/mu_mat_union3_cosmo=2_mu.fits | Rubin+2023 arXiv:2311.12098 |
| des_dovekie/DES-Dovekie_HD.csv | summary statistic | 148 kB | des-science/DES-SN5YR/main/.../DES-Dovekie_HD.csv | Popovic+2026 arXiv:2511.07517 |

## Downloaded on demand (catalogs and summary statistics, verified)

| dataset | file | kind | size | what |
|---|---|---|---|---|
| pantheonplus | Pantheon+SH0ES_STAT+SYS.cov (.npy cached) | summary statistic | 33 MB | 1701×1701 stat+sys covariance of mu |
| des5yr | STAT+SYS.txt.gz | summary statistic | 16 MB | 1829×1829 systematic covariance (stat added from MUERR_FINAL) |
| des5yr | DES-SN5YR_HD+MetaData.csv | catalog | 662 kB | per-SN SALT3 fit parameters + host properties |
| des5yr | DES5YR_SALT3_LCFIT.LCPLOT.gz | catalog | 7.4 MB | photometric light curves (data + SALT3 model) |
| des_dovekie | STAT+SYS.npz | summary statistic | 6 MB | inverse covariance, upper triangle |
| pantheonplus | per-SN .dat files | catalog | ~1-10 kB each | SNANA-format light curves (downloaded per-SN) |

Union3 format (rubind/union3_release README): FITS image, row 0 = redshift nodes,
column 0 = distance moduli, remaining 22x22 block = INVERSE covariance. The same numbers are
in cobaya's `sn_data/Union3/lcparam_full.txt` + `mag_covmat.txt`. No per-SN Union3 table is public.
