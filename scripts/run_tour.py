"""Run the sn-tour skill end to end in plain Python (all local) into output/tour/.

    MPLBACKEND=Agg .venv/bin/python scripts/run_tour.py
"""

import json
import os
import shutil
import time
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

from tools.meta import describe_sn_dataset, list_sn_datasets  # noqa: E402
from tools.sn import (  # noqa: E402
    bin_hubble_diagram,
    compare_sn_samples,
    compute_distance_modulus,
    fetch_sn_sample,
    fit_sn_cosmology,
    plot_hubble_diagram,
    plot_sn_contours,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = str(ROOT / "output" / "tour")
t0 = time.time()
summary = {}
print(list_sn_datasets().message)
print(describe_sn_dataset("pantheonplus").message)
f = fetch_sn_sample(OUT, dataset="pantheonplus")
print(f.message)
summary["pantheonplus_n"] = f.metadata["n"]
print(compute_distance_modulus(OUT, Om=0.3).metadata["mu_at_z1"])
h = plot_hubble_diagram(OUT, dataset="pantheonplus")
print(h.message)
print(bin_hubble_diagram(OUT, dataset="pantheonplus", n_bins=20).message)
for ds in ("pantheonplus", "des5yr", "union3"):
    r = fit_sn_cosmology(OUT, dataset=ds, model="lcdm")
    print(r.message)
    summary[f"{ds}_lcdm_Om"] = r.metadata["marginals"]["Om"]
w = fit_sn_cosmology(OUT, dataset="pantheonplus", model="wcdm")
print(w.message)
summary["pantheonplus_wcdm"] = w.metadata["marginals"]
c = plot_sn_contours(OUT, contour_csvs=[x for x in w.files if "_2d_Om_w" in x])
cmp = compare_sn_samples(OUT, mode="residuals")
print(cmp.message)
summary["elapsed_s"] = round(time.time() - t0, 1)
(Path(OUT) / "tour_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

figs = ROOT / "docs" / "figures"
figs.mkdir(parents=True, exist_ok=True)
shutil.copy(h.files[0], figs / "hubble_pantheonplus.png")
shutil.copy(c.files[0], figs / "contours_Om_w.png")
shutil.copy(next(x for x in cmp.files if x.endswith(".png")), figs / "compare_residuals.png")
