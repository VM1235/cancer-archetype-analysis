#!/usr/bin/env python3
"""Step 2b: how much of each PC is tumour identity beyond EMT state?

Tumours differ in which EMT states they contain, so a pure EMT axis would
also look "tumour-driven" in step 2's table. Here each PC is first centred
within each EMT state (state averages removed), then we ask what share of
the remaining spread is explained by tumour, and by assay.

  state             share of the PC explained by EMT state (as in step 2)
  tumour_raw        share explained by tumour (as in step 2, confounded)
  tumour_in_state   share of the within-state spread explained by tumour
  assay_in_state    same for assay

High tumour_in_state = tumour/batch identity that is not EMT.

Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/Mouse_PDAC_GSE173958/codes/02b_tumour_within_state.py --geneset all --tag main
"""
import argparse

import numpy as np
import pandas as pd

import scc_common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geneset", choices=["all", "ks"], default="all")
    ap.add_argument("--tag", default="main")
    args = ap.parse_args()
    run = C.run_dir(args.geneset, args.tag)
    z = np.load(run / "pca.npz")
    meta = pd.read_csv(C.PROC / "cell_metadata.csv").set_index("barcode")
    assert list(meta.index) == list(z["barcodes"])
    o = meta.iloc[z["fit_idx"]]
    s = z["scores_all"][z["fit_idx"]]
    state = o["cluster_annotation"].astype(str).values
    rows = []
    for j in range(10):
        v = s[:, j]
        r = v - pd.Series(v).groupby(state).transform("mean").values
        rows.append(dict(PC=j + 1, pct_variance=100 * z["evr"][j],
                         state=C.eta_squared(v, state),
                         tumour_raw=C.eta_squared(v, o["sample"]),
                         tumour_in_state=C.eta_squared(r, o["sample"]),
                         assay_in_state=C.eta_squared(r, o["assay"])))
    df = pd.DataFrame(rows).round(3)
    df.to_csv(run / "pc_drivers_within_state.csv", index=False)
    print(f"{args.geneset}/{args.tag}\n{df.to_string(index=False)}")


if __name__ == "__main__":
    main()
