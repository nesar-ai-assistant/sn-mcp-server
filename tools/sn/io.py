"""Readers for the SN releases (local side): resolve files via tools.common.fetch_file
(bundled data/ first, then the cache, then download) and parse with the kernel parsers
so local and dispatched runs see byte-identical data."""

from __future__ import annotations

import numpy as np

from ..common import data_cache_dir, fetch_file, get_cached
from . import kernels

DATASET_NAMES = ("pantheonplus", "des5yr", "union3", "des_dovekie")

_LOCAL_NAMES = {
    "pantheonplus": ("Pantheon+SH0ES.dat", "Pantheon+SH0ES_STAT+SYS.cov"),
    "des5yr": ("DES-SN5YR_HD.csv", "STAT+SYS.txt.gz"),
    "union3": ("mu_mat_union3_cosmo=2_mu.fits", None),
    "des_dovekie": ("DES-Dovekie_HD.csv", "STAT+SYS.npz"),
}


def dataset_files(dataset: str) -> tuple:
    """Local paths (hd, cov|None), downloading on first use."""
    if dataset not in kernels.SOURCES:
        raise ValueError(f"Unknown dataset {dataset!r}. Available: {', '.join(DATASET_NAMES)}")
    src = kernels.SOURCES[dataset]
    hd_name, cov_name = _LOCAL_NAMES[dataset]
    hd = fetch_file(src["hd"], hd_name, subdir=dataset)
    cov = fetch_file(src["cov"], cov_name, subdir=dataset) if cov_name else None
    return hd, cov


def _load_full(dataset: str) -> dict:
    hd, cov = dataset_files(dataset)
    mat = None
    if cov is not None and str(cov).endswith((".cov", ".gz")):
        # text covariances are slow to parse (33 MB for Pantheon+): cache as .npy once
        npy = data_cache_dir() / dataset / (cov.name + ".npy")
        if npy.is_file():
            mat = np.load(npy)
        else:
            mat = kernels._read_matrix_txt(cov)
            np.save(npy, mat)
    return kernels.parse_dataset(dataset, hd, cov, cov_matrix=mat)


def load_full(dataset: str) -> dict:
    return get_cached(f"sn:{dataset}", lambda: _load_full(dataset))


def load_sample(dataset: str, zmin: float | None = None, zmax: float = 10.0,
                calibrators: str = "include") -> dict:
    """Selected sample; zmin=None uses the release's standard cosmology cut (Pantheon+ z > 0.01)."""
    from ..registry import DATASETS

    if zmin is None:
        zmin = DATASETS[dataset].get("default_zmin", 0.0)
    return kernels.select(load_full(dataset), zmin, zmax, calibrators)
