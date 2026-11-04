# Bundled data (offline core, 0.8 MB total)

Each file is byte-identical to the upstream URL (md5-checked 2026-11-04). The large
covariances are NOT bundled: `tools.common.fetch_file` downloads them on first use into
`~/.cache/cosmo-mcp/sn-mcp-server/<dataset>/` (override with `$COSMO_MCP_DATA`).

| file | size | upstream URL | citation |
|---|---|---|---|
| pantheonplus/Pantheon+SH0ES.dat | 579 kB | https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat | Scolnic+2022 ApJ 938 113; Brout+2022 ApJ 938 110 |
| des5yr/DES-SN5YR_HD.csv | 92 kB | https://raw.githubusercontent.com/des-science/DES-SN5YR/1.3/4_DISTANCES_COVMAT/DES-SN5YR_HD.csv (tag 1.3 = original 2024 release) | DES Collaboration 2024 ApJL 973 L14; Vincenzi+2024 |
| union3/mu_mat_union3_cosmo=2_mu.fits | 8.6 kB | https://raw.githubusercontent.com/rubind/union3_release/main/mu_mat_union3_cosmo%3D2_mu.fits | Rubin+2023 arXiv:2311.12098 |
| des_dovekie/DES-Dovekie_HD.csv | 148 kB | https://raw.githubusercontent.com/des-science/DES-SN5YR/main/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv | Popovic+2026 arXiv:2511.07517 |

Downloaded on demand (verified):

| dataset | file | size | URL |
|---|---|---|---|
| pantheonplus | Pantheon+SH0ES_STAT+SYS.cov (1701x1701 text; cached also as .npy) | 33 MB | https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov |
| des5yr | STAT+SYS.txt.gz (1829x1829 systematics; stat added from MUERR_FINAL) | 16 MB | https://raw.githubusercontent.com/des-science/DES-SN5YR/1.3/4_DISTANCES_COVMAT/STAT%2BSYS.txt.gz |
| des_dovekie | STAT+SYS.npz (inverse covariance, upper triangle) | 6 MB | https://raw.githubusercontent.com/des-science/DES-SN5YR/main/4_DISTANCES_COVMAT/STAT%2BSYS.npz |

Union3 format (rubind/union3_release README): FITS image, row 0 = redshift nodes,
column 0 = distance moduli, remaining 22x22 block = INVERSE covariance. The same numbers are
in cobaya's `sn_data/Union3/lcparam_full.txt` + `mag_covmat.txt`. No per-SN Union3 table is public.
