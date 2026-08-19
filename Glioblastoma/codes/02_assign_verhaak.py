#!/usr/bin/env python3
"""Independent Verhaak-style labels for GBM cell lines (Panel B only).

Does not use PCHA. Each line is scored on compact marker sets for Classical,
Mesenchymal, and Proneural (Wang 2017 three-type scheme; Neural is omitted).
Score = mean of gene-wise z-scores on the DepMap Panel A matrix. Assign the
max score. This stays on one platform, unlike nearest-centroid to TCGA.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
GBM = HERE.parent
ROOT = GBM.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from src.io import load_expression_csv

CL_MATRIX = GBM / "data" / "processed" / "input_panelA_glioblastoma_ccle_logtpm_filtered.csv"
OUT = GBM / "results" / "panel_b"

MARKERS = {
    "Classical": ("EGFR", "NES", "NOTCH3", "SMO", "GAS1", "GLI2", "NFKBIA"),
    "Mesenchymal": ("CD44", "CHI3L1", "RELB", "TRADD", "TNFRSF1A", "NFKB1", "STAT3"),
    "Proneural": ("OLIG2", "PDGFRA", "DLL3", "NKX2-2", "SOX10", "ASCL1", "NKX2-1"),
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cl = load_expression_csv(CL_MATRIX)
    z = cl.sub(cl.mean(axis=1), axis=0).div(cl.std(axis=1).replace(0, np.nan), axis=0)

    rows = []
    for sid in cl.columns.astype(str):
        scores = {}
        used = {}
        for subtype, genes in MARKERS.items():
            present = [g for g in genes if g in z.index]
            used[subtype] = ",".join(present)
            scores[subtype] = float(z.loc[present, sid].mean()) if present else np.nan
        best = max(scores, key=lambda s: (scores[s] if np.isfinite(scores[s]) else -np.inf))
        rows.append(
            {
                "cell_line": sid,
                "verhaak_subtype": best,
                "confidence_score": scores[best],
                **{f"score_{s}": scores[s] for s in MARKERS},
                **{f"genes_{s}": used[s] for s in MARKERS},
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "verhaak_labels_panelA.csv", index=False)
    print("Marker genes present:")
    for subtype, genes in MARKERS.items():
        present = [g for g in genes if g in cl.index]
        print(f"  {subtype}: {len(present)}/{len(genes)} {present}")
    print("\nCell-line assignments:")
    print(out["verhaak_subtype"].value_counts().to_string())
    print("Wrote", OUT / "verhaak_labels_panelA.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
