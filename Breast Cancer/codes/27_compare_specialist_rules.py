#!/usr/bin/env python3
"""Compare three specialist / generalist rules on the frozen Pal k=3 triangle.

Rules
-----
1. Groves et al. 2022, Figure 3E. Within each dataset, the 5% of samples
   closest in PC1–PC2 to an archetype are specialists for that archetype.
   A sample in more than one such tail is assigned to the nearer archetype.
   Everyone else is a generalist.
2. Leading weight > 0.5. AAnet (Cancer Discovery) calls a cell committed
   when its affinity for one archetype exceeds 0.5, i.e. exceeds the sum of
   the other affinities. The Groves 2021 preprint used the same majority cut
   on archetype scores. Weights are the clipped, renormalized barycentric
   weights already stored by script 26.
3. Purity gap >= 0.35. Largest weight beats the second by at least 0.35.
   This is the rule in scripts 24 and 26. It is not taken from those papers.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/27_compare_specialist_rules.py"
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
SCORES = BREAST / "results" / "specialist_interior"
PAL_FIT = BREAST / "results" / "panel_a_ks_genelist_gse161529"
OUT = SCORES / "rule_comparison"

SOURCES = {
    "GSE161529": (
        BREAST / "results" / "panel_b_gse161529" / "emt_hybrid" / "pal_cells_with_program_state.csv",
        "barcode",
    ),
    "GSE176078": (
        BREAST / "results" / "gse176078_on_gse161529" / "gse176078_in_gse161529_k3.csv",
        "barcode",
    ),
    "GSE173634": (
        BREAST / "results" / "gse173634_on_gse161529" / "gse173634_in_gse161529_k3.csv",
        "barcode",
    ),
    "DepMap": (
        BREAST / "results" / "depmap_on_gse161529" / "depmap_in_gse161529_k3.csv",
        "cell_line",
    ),
}
PROGRAMS = ("epithelial", "mesenchymal", "hybrid", "intermediate")


def closest_five_percent(dist: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    """dist is n x 3 Euclidean distances to the three archetypes."""
    n = dist.shape[0]
    n_take = int(np.ceil(0.05 * n))
    chosen = np.zeros(n, dtype=bool)
    best = np.full(n, np.inf)
    arc = np.zeros(n, dtype=int)
    for j in range(3):
        idx = np.argpartition(dist[:, j], n_take - 1)[:n_take]
        chosen[idx] = True
        closer = dist[idx, j] < best[idx]
        pick = idx[closer]
        best[pick] = dist[pick, j]
        arc[pick] = j + 1
    return chosen, arc, n_take


def summarize(name: str, df: pd.DataFrame, rule: str, specialist: np.ndarray, corner: np.ndarray) -> list[dict]:
    rows = []
    n = len(df)
    hybrid = df["program_state"].values == "hybrid"
    generalist = ~specialist
    row = {
        "dataset": name,
        "rule": rule,
        "n": n,
        "n_specialist": int(specialist.sum()),
        "pct_specialist": 100.0 * specialist.mean(),
        "n_hybrid": int(hybrid.sum()),
        "pct_hybrid_that_are_specialist": 100.0 * specialist[hybrid].mean() if hybrid.any() else np.nan,
        "pct_hybrid_at_arc1": 100.0 * ((corner == 1) & hybrid).sum() / hybrid.sum() if hybrid.any() else np.nan,
        "pct_hybrid_at_arc2": 100.0 * ((corner == 2) & hybrid).sum() / hybrid.sum() if hybrid.any() else np.nan,
        "pct_hybrid_at_arc3": 100.0 * ((corner == 3) & hybrid).sum() / hybrid.sum() if hybrid.any() else np.nan,
        "pct_generalist_that_are_hybrid": 100.0 * hybrid[generalist].mean() if generalist.any() else np.nan,
    }
    for program in ("epithelial", "mesenchymal", "hybrid"):
        mask = df["program_state"].values == program
        row[f"pct_{program}_specialist"] = 100.0 * specialist[mask].mean() if mask.any() else np.nan
    if "patient" in df.columns and hybrid.any() and ((corner == 1) & hybrid).any():
        sub = df.loc[(corner == 1) & hybrid]
        top = sub["patient"].value_counts(normalize=True)
        row["top_patient_among_hybrid_arc1"] = str(top.index[0])
        row["pct_hybrid_arc1_from_top_patient"] = 100.0 * float(top.iloc[0])
    else:
        row["top_patient_among_hybrid_arc1"] = ""
        row["pct_hybrid_arc1_from_top_patient"] = np.nan
    rows.append(row)
    return rows


def main():
    arcs = np.load(PAL_FIT / "archetypes_k3_parti.npy")
    if arcs.shape != (3, 2):
        raise SystemExit(f"expected Pal k=3 archetypes as 3x2, got {arcs.shape}")
    OUT.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    enrich_rows = []

    for name, (path, key) in SOURCES.items():
        scores = pd.read_csv(SCORES / f"scores_{name}.csv")
        xy = pd.read_csv(path, usecols=[key, "PC1", "PC2"])
        df = scores.merge(xy, left_on="sample_id", right_on=key, how="left", validate="one_to_one")
        if df[["PC1", "PC2"]].isna().any().any():
            raise SystemExit(f"{name}: PC coordinates failed to join")
        pts = df[["PC1", "PC2"]].to_numpy(float)
        dist = np.linalg.norm(pts[:, None, :] - arcs[None, :, :], axis=2)
        # Sanity: weight on an archetype should fall as distance to it rises.
        for j, col in enumerate(("w_arc1", "w_arc2", "w_arc3")):
            rho = np.corrcoef(df[col], dist[:, j])[0, 1]
            print(f"{name} corr({col}, distance to arc {j + 1}) = {rho:.3f}")

        spec5, arc5, n_take = closest_five_percent(dist)
        # overlap: a sample in two tails before the nearest-arc assignment
        in_tail = np.zeros(len(df), dtype=int)
        for j in range(3):
            idx = np.argpartition(dist[:, j], n_take - 1)[:n_take]
            in_tail[idx] += 1
        n_overlap = int((in_tail > 1).sum())

        w = df[["w_arc1", "w_arc2", "w_arc3"]].to_numpy(float)
        dominant = w.argmax(axis=1) + 1
        spec_half = w.max(axis=1) > 0.5
        arc_half = np.where(spec_half, dominant, 0)

        spec_gap = df["region"].to_numpy() == "specialist"
        arc_gap = df["corner_arc"].to_numpy(int)

        rules = {
            "groves_closest_5pct": (spec5, arc5),
            "leading_weight_gt_0.5": (spec_half, arc_half),
            "purity_gap_ge_0.35": (spec_gap, arc_gap),
        }
        for rule, (spec, arc) in rules.items():
            summary_rows.extend(summarize(name, df, rule, spec, arc))
            for j in range(3):
                if rule == "groves_closest_5pct":
                    tail = np.argpartition(dist[:, j], n_take - 1)[:n_take]
                    members = df.iloc[tail]
                else:
                    members = df.loc[arc == j + 1]
                n_m = len(members)
                enrich_rows.append(
                    {
                        "dataset": name,
                        "rule": rule,
                        "arc": j + 1,
                        "n_specialist_at_arc": n_m,
                        "pct_of_those_hybrid": 100.0 * (members["program_state"] == "hybrid").mean() if n_m else np.nan,
                        "pct_of_dataset_hybrid": 100.0 * (df["program_state"] == "hybrid").mean(),
                        "n_take_per_arc": n_take if rule == "groves_closest_5pct" else n_m,
                        "n_samples_in_two_tails": n_overlap if rule == "groves_closest_5pct" else 0,
                    }
                )
        print(
            f"{name}: n={len(df)}  5% tail per arc={n_take}  "
            f"samples in two tails={n_overlap}  "
            f"specialists 5%={spec5.sum()}  >0.5={spec_half.sum()}  gap={spec_gap.sum()}"
        )

    summary = pd.DataFrame(summary_rows)
    enrich = pd.DataFrame(enrich_rows)
    summary.to_csv(OUT / "three_rules_summary.csv", index=False)
    enrich.to_csv(OUT / "three_rules_by_arc.csv", index=False)
    print("\nWrote", OUT / "three_rules_summary.csv")
    cols = [
        "dataset",
        "rule",
        "pct_specialist",
        "pct_epithelial_specialist",
        "pct_mesenchymal_specialist",
        "pct_hybrid_that_are_specialist",
        "pct_hybrid_at_arc1",
        "pct_generalist_that_are_hybrid",
        "top_patient_among_hybrid_arc1",
        "pct_hybrid_arc1_from_top_patient",
    ]
    print(summary[cols].to_string(index=False, float_format=lambda v: f"{v:.1f}"))


if __name__ == "__main__":
    main()
