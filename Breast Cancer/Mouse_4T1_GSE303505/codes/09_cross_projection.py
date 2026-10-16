#!/usr/bin/env python3
"""Step 9: project every library into one library's frozen simplex.

Question: are the corners found inside one library (the host) shared by the others?
Guests are centred on the HOST's gene means and projected on the host's PCs, with no
rescaling (so each library keeps its own offset, on purpose). Each corner is described by
its signature in the host (closest 10% of host cells minus the rest), its top genes, and
the KS_EMT / cycling / depth profile of the host cells that sit there.

Outputs (results/cross_projection/<geneset>_<host><suffix>_k<k>/):
  cells.csv.gz         host-space PCs, weights, corner/interior, inside flag
  summary.md           per library: % at each host corner, % interior, % inside the simplex
                       per library x EMT state: same
Usage (from Project_1/):
  .venv/bin/python "Breast Cancer/Mouse_4T1_GSE303505"/codes/09_cross_projection.py --host M1-PT --k 4 --clean
"""
import argparse
import numpy as np, pandas as pd
import scc_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--host", required=True)
ap.add_argument("--geneset", default="all")
ap.add_argument("--k", type=int, default=4)
ap.add_argument("--clean", action="store_true")
args = ap.parse_args()
k = args.k
sfx = "_clean" if args.clean else ""
d = C.RESULTS / "within_unit" / f"{args.geneset}_{args.host}{sfx}"
out = C.RESULTS / "cross_projection" / f"{args.geneset}_{args.host}{sfx}_k{k}"
out.mkdir(parents=True, exist_ok=True)

z = np.load(d / "pca.npz")
arcs = np.load(d / f"archetypes_k{k}.npy")
a = C.load_adata()
obs = a.obs
gidx = [a.var_names.get_loc(g) for g in z["genes"]]
V = z["components"][: k - 1]
pts = np.asarray(a.X[:, gidx] @ V.T) - z["mean"] @ V.T
hr = z["rows"]
assert np.abs(pts[hr] - z["scores"][:, : k - 1]).max() < 1e-4 * np.abs(z["scores"]).max(), "host mismatch"

w = C.barycentric(pts, arcs)
pur, near = C.purity_and_nearest(w)
Vn = [f"V{i + 1}" for i in range(k)]
group = np.where(pur >= C.PURITY_CUT, np.array(Vn)[near], "interior")
dist = np.linalg.norm(pts[hr][:, None, :] - arcs[None, :, : k - 1], axis=2)
st = obs["cluster_annotation"].astype(str).values
ph = obs["phase"].astype(str).values
names = {}
lines = [f"# Cross-projection into {args.host} (genes {args.geneset}{sfx}, k={k})", "",
         "## Host corners (closest 10% of the host's cells)", ""]
for v in range(k):
    close = hr[dist[:, v] <= np.quantile(dist[:, v], 0.10)]
    close_n = hr[dist[:, v] > np.quantile(dist[:, v], 0.10)]
    sig = np.asarray(a.X[close].mean(0)).ravel() - np.asarray(a.X[close_n].mean(0)).ravel()
    top = ", ".join(a.var_names[np.argsort(sig)[::-1][:6]])
    cyc = 100 * np.isin(ph[close], ["S", "G2M"]).mean()
    ks = obs["KS_EMT"].values[close].mean()
    names[Vn[v]] = f"{Vn[v]} ({top.split(', ')[0]}...)"
    lines.append(f"- {Vn[v]}: top genes {top}; cycling {cyc:.0f}%; KS_EMT {ks:+.2f}; "
                 f"UMI {obs['nCount_RNA'].values[close].mean():.0f}")
cells = pd.DataFrame({"sample": obs["sample"].astype(str).values, "state": st, "group": group,
                      "inside": (w >= -1e-6).all(1), "is_host": np.isin(np.arange(len(obs)), hr)},
                     index=obs.index)
for i in range(k - 1):
    cells[f"PC{i + 1}"] = pts[:, i]
cells.to_csv(out / "cells.csv.gz")
order = ["E", "EM1", "EM2", "EM3", "M", "4T1"]
tab = []
for lib in order:
    m = cells["sample"] == lib
    g = cells.loc[m, "group"].value_counts(normalize=True)
    tab.append(dict(library=lib + (" (host)" if lib == args.host else ""), n=int(m.sum()),
                    **{v: round(100 * g.get(v, 0), 1) for v in Vn}, interior=round(100 * g.get("interior", 0), 1),
                    inside=round(100 * cells.loc[m, "inside"].mean(), 1)))
tab = pd.DataFrame(tab)
tab.to_csv(out / "library_in_host.csv", index=False)
lines += ["", "## % of each library's cells at each host corner (purity gap >= 0.35), interior, and inside", "",
          tab.to_string(index=False)]
sub = cells[~cells["is_host"] & True]
pool = (100 * pd.crosstab(sub["state"], sub["group"], normalize="index")).round(1)
pool = pool.reindex(C.STATE_ORDER)
lines += ["", "## All cells of the other clones, by EMT state: % at each host corner / interior", "", pool.to_string()]
(out / "summary.md").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
