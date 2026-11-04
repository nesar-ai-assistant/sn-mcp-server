# Bundled data (offline core, ~1 MB total)

Each file is byte-identical to the upstream URL (md5-checked 2026-11-04). The large
covariances are NOT bundled: `tools.common.fetch_file` downloads them on first use into
`~/.cache/cosmo-mcp/sn-mcp-server/<dataset>/` (override with `$COSMO_MCP_DATA`).

## Data tiers
| tier | what | examples in this server |
|---|---|---|
| 1 | raw / catalog-level | per-SN SALT2/SALT3 fit params (x1, c, mB, HOST_LOGMASS), photometric light curves |
| 2 | summary statistics + covariance | Hubble diagrams mu(z) with full stat+sys covariance |
| 3 | published inference | posterior chains (comparison only, not used as analysis input) |

## Bundled files

| file | tier | size | upstream URL | citation |
|---|---|---|---|---|
| pantheonplus/Pantheon+SH0ES.dat | 1+2 | 579 kB | PantheonPlusSH0ES/DataRelease .../Pantheon%2BSH0ES.dat | Scolnic+2022 ApJ 938 113; Brout+2022 ApJ 938 110 |
| pantheonplus/lc_index.csv | 1 | 40 kB | built from PantheonPlusSH0ES/DataRelease/1_DATA/photometry tree | — |
| des5yr/DES-SN5YR_HD.csv | 2 | 92 kB | des-science/DES-SN5YR/1.3/4_DISTANCES_COVMAT/DES-SN5YR_HD.csv | DES 2024 ApJL 973 L14; Vincenzi+2024 |
| union3/mu_mat_union3_cosmo=2_mu.fits | 2 | 8.6 kB | rubind/union3_release/mu_mat_union3_cosmo=2_mu.fits | Rubin+2023 arXiv:2311.12098 |
| des_dovekie/DES-Dovekie_HD.csv | 2 | 148 kB | des-science/DES-SN5YR/main/.../DES-Dovekie_HD.csv | Popovic+2026 arXiv:2511.07517 |

## Downloaded on demand (tier 1-2, verified)

| dataset | file | tier | size | what |
|---|---|---|---|---|
| pantheonplus | Pantheon+SH0ES_STAT+SYS.cov (.npy cached) | 2 | 33 MB | 1701×1701 stat+sys covariance of mu |
| des5yr | STAT+SYS.txt.gz | 2 | 16 MB | 1829×1829 systematic covariance (stat added from MUERR_FINAL) |
| des5yr | DES-SN5YR_HD+MetaData.csv | 1 | 662 kB | per-SN SALT3 fit parameters + host properties |
| des5yr | DES5YR_SALT3_LCFIT.LCPLOT.gz | 1 | 7.4 MB | photometric light curves (data + SALT3 model) |
| des_dovekie | STAT+SYS.npz | 2 | 6 MB | inverse covariance, upper triangle |
| pantheonplus | per-SN .dat files | 1 | ~1-10 kB each | SNANA-format light curves (downloaded per-SN) |

Union3 format (rubind/union3_release README): FITS image, row 0 = redshift nodes,
column 0 = distance moduli, remaining 22x22 block = INVERSE covariance. The same numbers are
in cobaya's `sn_data/Union3/lcparam_full.txt` + `mag_covmat.txt`. No per-SN Union3 table is public.
