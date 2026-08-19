"""Shared helpers for hypothesis-driven GBM k=2 (PN–MES axis).

Not used by the Groves k=3–7 pipeline. Marker lists match 02_assign_verhaak.py.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

MARKERS = {
    "Classical": ("EGFR", "NES", "NOTCH3", "SMO", "GAS1", "GLI2", "NFKBIA"),
    "Mesenchymal": ("CD44", "CHI3L1", "RELB", "TRADD", "TNFRSF1A", "NFKB1", "STAT3"),
    "Proneural": ("OLIG2", "PDGFRA", "DLL3", "NKX2-2", "SOX10", "ASCL1", "NKX2-1"),
}

SUBTYPE_COLOR = {
    "Classical": "#3B6FA0",
    "Mesenchymal": "#E07A3D",
    "Proneural": "#54A24B",
    "Neural": "#D989B5",
}


def signature_scores(expr):
    """Gene-wise z-score on this matrix, then mean of each marker set.

    expr: genes × samples. Returns samples × {score_Classical, ...}.
    """
    z = expr.sub(expr.mean(axis=1), axis=0).div(
        expr.std(axis=1).replace(0, np.nan), axis=0
    )
    rows = []
    present = {}
    for subtype, genes in MARKERS.items():
        hit = [g for g in genes if g in z.index]
        present[subtype] = hit
        if hit:
            rows.append(z.loc[hit].mean(axis=0).rename(f"score_{subtype}"))
        else:
            rows.append(pd.Series(np.nan, index=expr.columns, name=f"score_{subtype}"))
    out = pd.concat(rows, axis=1)
    out.index = expr.columns.astype(str)
    out.index.name = "sample"
    return out, present


def assign_poles(weights, scores):
    """Map k=2 mixture columns to PN vs MES using Pearson of S vs signatures.

    weights: (n, 2). scores: DataFrame with score_Proneural / score_Mesenchymal.
    Vertex i is the MES pole if corr(S_i, MES − PN) is larger than for the other
    vertex (equivalently: the column more aligned with the MES–PN axis).
    """
    s = np.asarray(weights, dtype=float)
    pn = scores["score_Proneural"].to_numpy(dtype=float)
    mes = scores["score_Mesenchymal"].to_numpy(dtype=float)
    axis = mes - pn
    mask = np.isfinite(axis)
    corrs = []
    for j in range(s.shape[1]):
        if mask.sum() < 5:
            corrs.append(np.nan)
        else:
            corrs.append(float(np.corrcoef(s[mask, j], axis[mask])[0, 1]))
    mes_idx = int(np.nanargmax(corrs))
    pn_idx = 1 - mes_idx
    return {
        "mes_idx": mes_idx,
        "pn_idx": pn_idx,
        "corr_S_with_MES_minus_PN": corrs,
        "corr_S0_PN": _safe_pearson(s[:, 0], pn),
        "corr_S0_MES": _safe_pearson(s[:, 0], mes),
        "corr_S1_PN": _safe_pearson(s[:, 1], pn),
        "corr_S1_MES": _safe_pearson(s[:, 1], mes),
    }


def _safe_pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        return np.nan
    r, p = pearsonr(a[m], b[m])
    return float(r), float(p)


def corr_table(weight, scores, label="w"):
    rows = []
    w = np.asarray(weight, dtype=float)
    for col in scores.columns:
        y = scores[col].to_numpy(dtype=float)
        m = np.isfinite(w) & np.isfinite(y)
        if m.sum() < 5:
            rows.append(
                {
                    "weight": label,
                    "signature": col,
                    "n": int(m.sum()),
                    "pearson_r": np.nan,
                    "pearson_p": np.nan,
                    "spearman_r": np.nan,
                    "spearman_p": np.nan,
                }
            )
            continue
        pr, pp = pearsonr(w[m], y[m])
        sr, sp = spearmanr(w[m], y[m])
        rows.append(
            {
                "weight": label,
                "signature": col,
                "n": int(m.sum()),
                "pearson_r": float(pr),
                "pearson_p": float(pp),
                "spearman_r": float(sr),
                "spearman_p": float(sp),
            }
        )
    return pd.DataFrame(rows)


def line_barycentric(points, v0, v1):
    """Mixture weights on the line through two vertices (any dimension).

    w1 = projection coordinate of v1; w0 = 1 - w1.
    A point is between the poles iff both weights are >= 0 (equiv. t in [0, 1]).
    """
    pts = np.asarray(points, dtype=float)
    a = np.asarray(v0, dtype=float).ravel()
    b = np.asarray(v1, dtype=float).ravel()
    d = b - a
    denom = float(np.dot(d, d))
    if denom <= 0:
        w1 = np.full(pts.shape[0], np.nan)
    else:
        w1 = ((pts - a) @ d) / denom
    w0 = 1.0 - w1
    return np.column_stack([w0, w1])


def style_mpl(plt):
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.bbox": "tight",
            "savefig.dpi": 240,
        }
    )
