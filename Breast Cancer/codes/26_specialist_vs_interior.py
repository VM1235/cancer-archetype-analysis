#!/usr/bin/env python3
"""Specialist vs interior from k=3 archetype weights, crossed with Tan KS hybrid state.

Priority from the Pal-triangle note: score GSE161529, GSE176078, GSE173634, and
the DepMap lines the same way, and ask whether hybrid E/M cells are interior
generalists or sit at a corner.

Geometry (matches the existing Pal barycentric vertex rule in
24_panelB_emt_hybrid_gse161529.py, which gave "66% of hybrid cells at a corner,
53% at Arc 1"):

  purity = largest weight − second-largest weight
  specialist at the dominant archetype if purity >= 0.35
  interior otherwise (weights shared by two or three archetypes)

Weights that fall slightly outside the simplex are clipped at 0 and
renormalized before this cut. Program labels use the same within-dataset
rule as script 24: hybrid if both epi and mes scores are at or above the
75th percentile; otherwise epithelial / mesenchymal by EMT tertiles.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/26_specialist_vs_interior.py"
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GENELIST = BREAST / "data" / "processed" / "ks_cellline_signature_tan2014.csv"
DEPMAP = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
PAL_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse161529_magic" / "ks_genes_magic_imputed.npz"
)
WU_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse176078_magic" / "ks_genes_magic_imputed.npz"
)
G173_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse173634_magic" / "ks_genes_magic_imputed.npz"
)
PAL_CELLS = BREAST / "results" / "panel_b_gse161529" / "emt_hybrid" / "pal_cells_with_program_state.csv"
WU_PROJ = BREAST / "results" / "gse176078_on_gse161529" / "gse176078_in_gse161529_k3.csv"
G173_PROJ = BREAST / "results" / "gse173634_on_gse161529" / "gse173634_in_gse161529_k3.csv"
DEPMAP_PROJ = BREAST / "results" / "depmap_on_gse161529" / "depmap_in_gse161529_k3.csv"

OUT = BREAST / "results" / "specialist_interior"
FIG = BREAST / "figures"

PURITY_CUT = 0.35
PROGRAM_LEVELS = ("epithelial", "mesenchymal", "hybrid", "intermediate")
ARC_COLORS = {1: "#4C78A8", 2: "#F58518", 3: "#E45756"}
INTERIOR_COLOR = "#B0B0B0"
PROGRAM_COLORS = {
    "epithelial": "#4C78A8",
    "mesenchymal": "#E45756",
    "hybrid": "#9C755F",
    "intermediate": "#B0B0B0",
}
DATASETS = ("GSE161529", "GSE176078", "GSE173634", "DepMap")


def style():
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


def zscore_genes(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).replace(0, np.nan).fillna(1.0).clip(lower=1e-8)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def ks_cell_scores(expr: pd.DataFrame, epi: list[str], mes: list[str]) -> pd.DataFrame:
    z = zscore_genes(expr)
    epi_g = [g for g in epi if g in z.index]
    mes_g = [g for g in mes if g in z.index]
    epi_s = z.loc[epi_g].mean(axis=0)
    mes_s = z.loc[mes_g].mean(axis=0)
    out = pd.DataFrame({"epi_score": epi_s, "mes_score": mes_s})
    out["emt_score"] = out["mes_score"] - out["epi_score"]
    out.attrs["n_epi"] = len(epi_g)
    out.attrs["n_mes"] = len(mes_g)
    return out


def classify_program(epi: np.ndarray, mes: np.ndarray, emt: np.ndarray) -> np.ndarray:
    epi_q75, mes_q75 = np.quantile(epi, 0.75), np.quantile(mes, 0.75)
    emt_lo, emt_hi = np.quantile(emt, [1 / 3, 2 / 3])
    hybrid = (epi >= epi_q75) & (mes >= mes_q75)
    out = np.full(len(emt), "intermediate", dtype=object)
    out[hybrid] = "hybrid"
    rest = ~hybrid
    out[rest & (emt <= emt_lo)] = "epithelial"
    out[rest & (emt >= emt_hi)] = "mesenchymal"
    return out


def composition(w: np.ndarray) -> np.ndarray:
    w = np.clip(np.asarray(w, dtype=float), 0, None)
    total = w.sum(axis=1, keepdims=True)
    total[total == 0] = 1.0
    return w / total


def geometry_from_weights(w: np.ndarray, purity_cut: float = PURITY_CUT) -> pd.DataFrame:
    w = composition(w)
    order = np.argsort(-w, axis=1)
    top = order[:, 0] + 1
    sorted_w = np.sort(w, axis=1)
    max_w = sorted_w[:, -1]
    purity = sorted_w[:, -1] - sorted_w[:, -2]
    specialist = purity >= purity_cut
    # Shannon evenness: 0 = one archetype, 1 = equal weights.
    wc = np.clip(w, 1e-12, None)
    evenness = -(wc * np.log(wc)).sum(axis=1) / np.log(w.shape[1])
    # 0 at the centroid (max=1/k), 1 at a vertex (max=1).
    specialization = (w.shape[1] * max_w - 1.0) / (w.shape[1] - 1.0)
    region = np.where(specialist, "specialist", "interior")
    corner = np.where(specialist, top, 0)
    return pd.DataFrame(
        {
            "w_arc1": w[:, 0],
            "w_arc2": w[:, 1],
            "w_arc3": w[:, 2],
            "max_weight": max_w,
            "purity": purity,
            "evenness": evenness,
            "specialization": specialization,
            "region": region,
            "corner_arc": corner,
            "dominant_arc": top,
        }
    )


def load_magic(path: Path) -> pd.DataFrame:
    z = np.load(path, allow_pickle=True)
    return pd.DataFrame(
        z["imputed"],
        index=[str(g) for g in z["genes"]],
        columns=[str(b) for b in z["barcodes"]],
    )


def build_frames(epi_genes: list[str], mes_genes: list[str]) -> dict[str, pd.DataFrame]:
    frames = {}

    pal = pd.read_csv(
        PAL_CELLS,
        usecols=["barcode", "patient", "subtype", "w_arc1", "w_arc2", "w_arc3", "program_state"],
    )
    pal_expr = load_magic(PAL_MAGIC)
    pal_ks = ks_cell_scores(pal_expr, epi_genes, mes_genes)
    pal = pal.merge(pal_ks, left_on="barcode", right_index=True, how="left", validate="one_to_one")
    if pal[["epi_score", "mes_score"]].isna().any().any():
        raise ValueError("Pal KS scores failed to join on barcode")
    recomputed = classify_program(
        pal["epi_score"].values, pal["mes_score"].values, pal["emt_score"].values
    )
    agree = float((recomputed == pal["program_state"].values).mean())
    print(f"Pal program labels recomputed vs saved: {agree:.1%} agreement")
    pal["program_state"] = recomputed
    pal["dataset"] = "GSE161529"
    pal["sample_id"] = pal["barcode"]
    frames["GSE161529"] = pal

    wu = pd.read_csv(WU_PROJ)
    wu_ks = ks_cell_scores(load_magic(WU_MAGIC), epi_genes, mes_genes)
    wu = wu.merge(wu_ks, left_on="barcode", right_index=True, how="left", validate="one_to_one")
    if wu[["epi_score", "mes_score"]].isna().any().any():
        raise ValueError("Wu KS scores failed to join on barcode")
    wu["program_state"] = classify_program(
        wu["epi_score"].values, wu["mes_score"].values, wu["emt_score"].values
    )
    wu["dataset"] = "GSE176078"
    wu["sample_id"] = wu["barcode"]
    frames["GSE176078"] = wu

    g173 = pd.read_csv(G173_PROJ)
    g173_ks = ks_cell_scores(load_magic(G173_MAGIC), epi_genes, mes_genes)
    g173 = g173.merge(g173_ks, left_on="barcode", right_index=True, how="left", validate="one_to_one")
    if g173[["epi_score", "mes_score"]].isna().any().any():
        raise ValueError("GSE173634 KS scores failed to join on barcode")
    g173["program_state"] = classify_program(
        g173["epi_score"].values, g173["mes_score"].values, g173["emt_score"].values
    )
    g173["dataset"] = "GSE173634"
    g173["sample_id"] = g173["barcode"]
    g173["patient"] = g173["cell_line"]
    frames["GSE173634"] = g173

    dep = pd.read_csv(DEPMAP_PROJ)
    dep_expr = pd.read_csv(DEPMAP, index_col=0)
    dep_expr.index = dep_expr.index.astype(str)
    dep_ks = ks_cell_scores(dep_expr, epi_genes, mes_genes)
    dep = dep.merge(dep_ks, left_on="cell_line", right_index=True, how="left", validate="one_to_one")
    if dep[["epi_score", "mes_score"]].isna().any().any():
        raise ValueError("DepMap KS scores failed to join on cell_line")
    dep["program_state"] = classify_program(
        dep["epi_score"].values, dep["mes_score"].values, dep["emt_score"].values
    )
    dep["dataset"] = "DepMap"
    dep["sample_id"] = dep["cell_line"]
    dep["subtype"] = dep["pam50"]
    frames["DepMap"] = dep
    return frames


def attach_geometry(frame: pd.DataFrame) -> pd.DataFrame:
    geo = geometry_from_weights(frame[["w_arc1", "w_arc2", "w_arc3"]].values)
    out = frame.drop(columns=["w_arc1", "w_arc2", "w_arc3"]).reset_index(drop=True)
    return pd.concat([out, geo], axis=1)


def conditional_tables(scored: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    by_program_rows = []
    p_program_given_region = []
    p_region_given_program = []
    for name, df in scored.items():
        n = len(df)
        for program in PROGRAM_LEVELS:
            sub = df[df["program_state"] == program]
            n_p = len(sub)
            if n_p == 0:
                continue
            by_program_rows.append(
                {
                    "dataset": name,
                    "program": program,
                    "n": n_p,
                    "frac_of_dataset": n_p / n,
                    "p_specialist": float((sub["region"] == "specialist").mean()),
                    "p_interior": float((sub["region"] == "interior").mean()),
                    "p_corner_arc1": float((sub["corner_arc"] == 1).mean()),
                    "p_corner_arc2": float((sub["corner_arc"] == 2).mean()),
                    "p_corner_arc3": float((sub["corner_arc"] == 3).mean()),
                    "median_max_weight": float(sub["max_weight"].median()),
                    "median_specialization": float(sub["specialization"].median()),
                    "median_evenness": float(sub["evenness"].median()),
                }
            )
            p_region_given_program.append(
                {
                    "dataset": name,
                    "given": program,
                    "n": n_p,
                    "p_specialist": float((sub["region"] == "specialist").mean()),
                    "p_interior": float((sub["region"] == "interior").mean()),
                }
            )
        for region in ("specialist", "interior"):
            sub = df[df["region"] == region]
            n_r = len(sub)
            row = {"dataset": name, "given": region, "n": n_r}
            for program in PROGRAM_LEVELS:
                row[f"p_{program}"] = float((sub["program_state"] == program).mean()) if n_r else np.nan
            p_program_given_region.append(row)
    return (
        pd.DataFrame(by_program_rows),
        pd.DataFrame(p_program_given_region),
        pd.DataFrame(p_region_given_program),
    )


def sensitivity(scored: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    cuts = (0.20, 0.35, 0.50)
    max_cuts = (0.50, 0.60, 0.70)
    for name, df in scored.items():
        w = df[["w_arc1", "w_arc2", "w_arc3"]].values
        hybrid = df["program_state"].values == "hybrid"
        for cut in cuts:
            geo = geometry_from_weights(w, purity_cut=cut)
            spec = geo["region"].values == "specialist"
            arc1 = geo["corner_arc"].values == 1
            rows.append(
                {
                    "dataset": name,
                    "rule": f"purity>={cut:.2f}",
                    "p_specialist_all": float(spec.mean()),
                    "p_specialist_given_hybrid": float(spec[hybrid].mean()) if hybrid.any() else np.nan,
                    "p_arc1_given_hybrid": float(arc1[hybrid].mean()) if hybrid.any() else np.nan,
                    "n_hybrid": int(hybrid.sum()),
                }
            )
        for cut in max_cuts:
            spec = df["max_weight"].values >= cut
            # corner = dominant arc among those passing the max-weight cut
            arc1 = spec & (df["dominant_arc"].values == 1)
            rows.append(
                {
                    "dataset": name,
                    "rule": f"max_weight>={cut:.2f}",
                    "p_specialist_all": float(spec.mean()),
                    "p_specialist_given_hybrid": float(spec[hybrid].mean()) if hybrid.any() else np.nan,
                    "p_arc1_given_hybrid": float(arc1[hybrid].mean()) if hybrid.any() else np.nan,
                    "n_hybrid": int(hybrid.sum()),
                }
            )
    return pd.DataFrame(rows)


def pal_arc1_patients(pal: pd.DataFrame) -> pd.DataFrame:
    hybrid_arc1 = pal[(pal["program_state"] == "hybrid") & (pal["corner_arc"] == 1)]
    if hybrid_arc1.empty:
        return pd.DataFrame()
    ct = hybrid_arc1.groupby(["patient", "subtype"], observed=True).size().reset_index(name="n")
    ct["frac_of_hybrid_arc1"] = ct["n"] / ct["n"].sum()
    ct = ct.sort_values("n", ascending=False)
    return ct


def plot(by_program: pd.DataFrame, given_region: pd.DataFrame, path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))

    ax = axes[0]
    hybrid = by_program[by_program["program"] == "hybrid"].set_index("dataset").reindex(DATASETS)
    for col, _, _ in [
        ("p_corner_arc1", "", ""),
        ("p_corner_arc2", "", ""),
        ("p_corner_arc3", "", ""),
        ("p_interior", "", ""),
    ]:
        hybrid[col] = hybrid[col].fillna(0.0)
    x = np.arange(len(DATASETS))
    pieces = [
        ("p_corner_arc1", "specialist Arc 1 (hybrid pole)", ARC_COLORS[1]),
        ("p_corner_arc2", "specialist Arc 2 (epi pole)", ARC_COLORS[2]),
        ("p_corner_arc3", "specialist Arc 3 (mes pole)", ARC_COLORS[3]),
        ("p_interior", "interior", INTERIOR_COLOR),
    ]
    bottom = np.zeros(len(DATASETS))
    for col, label, color in pieces:
        vals = hybrid[col].values.astype(float)
        ax.bar(x, vals, bottom=bottom, color=color, width=0.72, label=label)
        bottom += vals
    ax.set_xticks(x)
    labels = []
    for n in DATASETS:
        n_h = hybrid.loc[n, "n"] if n in hybrid.index else np.nan
        labels.append(f"{n}\n(n={0 if pd.isna(n_h) else int(n_h)})")
    ax.set_xticklabels(labels, fontsize=8)
    for i, n in enumerate(DATASETS):
        if n not in hybrid.index or pd.isna(hybrid.loc[n, "n"]):
            ax.text(i, 0.04, "no hybrid\ncells", ha="center", va="bottom", fontsize=7, color="0.35")
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of hybrid cells")
    ax.set_title("Where hybrid E/M cells sit")
    ax.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)

    ax = axes[1]
    width = 0.18
    offsets = np.linspace(-1.5, 1.5, len(DATASETS)) * width
    for name, off in zip(DATASETS, offsets):
        sub = by_program[by_program["dataset"] == name].set_index("program").reindex(PROGRAM_LEVELS)
        vals = sub["p_specialist"].fillna(0).values
        ax.bar(
            np.arange(len(PROGRAM_LEVELS)) + off,
            vals,
            width=width,
            label=name,
        )
    ax.axhline(0.5, color="0.75", lw=1, ls="--")
    ax.set_xticks(np.arange(len(PROGRAM_LEVELS)))
    ax.set_xticklabels(PROGRAM_LEVELS, rotation=15, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("P(specialist | program)")
    ax.set_title("Corner specialists by KS program")
    ax.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2)

    fig.suptitle(
        "k=3 Pal triangle — specialist if largest archetype weight exceeds the next by ≥ 0.35",
        y=1.02,
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)

    # Second figure: P(program | region) for each dataset.
    fig, axes = plt.subplots(1, 4, figsize=(12.2, 3.6), sharey=True)
    for ax, name in zip(axes, DATASETS):
        sub = given_region[given_region["dataset"] == name].set_index("given")
        x = np.arange(2)
        bottom = np.zeros(2)
        for program in PROGRAM_LEVELS:
            vals = np.array(
                [sub.loc["specialist", f"p_{program}"], sub.loc["interior", f"p_{program}"]],
                dtype=float,
            )
            ax.bar(x, vals, bottom=bottom, color=PROGRAM_COLORS[program], width=0.7, label=program)
            bottom += vals
        ax.set_xticks(x)
        labels = []
        for region in ("specialist", "interior"):
            labels.append(f"{region}\n(n={int(sub.loc[region, 'n'])})")
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(name)
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("fraction of cells")
    axes[-1].legend(
        fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0)
    )
    fig.suptitle("P(KS program | specialist or interior)", y=1.05, fontsize=11)
    fig.tight_layout()
    fig.savefig(path.with_name("Figure_program_given_specialist.png"))
    plt.close(fig)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    style()

    sig = pd.read_csv(GENELIST)
    epi_genes = sig.loc[sig["category"] == "Epi", "gene_symbol"].astype(str).tolist()
    mes_genes = sig.loc[sig["category"] == "Mes", "gene_symbol"].astype(str).tolist()

    frames = build_frames(epi_genes, mes_genes)
    scored = {name: attach_geometry(df) for name, df in frames.items()}

    keep = [
        "dataset",
        "sample_id",
        "subtype",
        "program_state",
        "epi_score",
        "mes_score",
        "emt_score",
        "w_arc1",
        "w_arc2",
        "w_arc3",
        "max_weight",
        "purity",
        "evenness",
        "specialization",
        "region",
        "corner_arc",
        "dominant_arc",
    ]
    for name, df in scored.items():
        cols = [c for c in keep if c in df.columns]
        if "patient" in df.columns:
            cols = cols[:2] + ["patient"] + cols[2:]
        df[cols].to_csv(OUT / f"scores_{name}.csv", index=False)
        n_spec = int((df["region"] == "specialist").sum())
        corr = float(np.corrcoef(df["epi_score"], df["mes_score"])[0, 1])
        print(
            f"{name}: n={len(df)} specialist {n_spec}/{len(df)} "
            f"({n_spec / len(df):.1%}); hybrid {(df.program_state == 'hybrid').mean():.1%}; "
            f"epi–mes r={corr:.2f}"
        )

    by_program, given_region, given_program = conditional_tables(scored)
    by_program.to_csv(OUT / "by_program.csv", index=False)
    given_region.to_csv(OUT / "p_program_given_region.csv", index=False)
    given_program.to_csv(OUT / "p_region_given_program.csv", index=False)

    sens = sensitivity(scored)
    sens.to_csv(OUT / "threshold_sensitivity.csv", index=False)

    patients = pal_arc1_patients(scored["GSE161529"])
    patients.to_csv(OUT / "pal_hybrid_arc1_by_patient.csv", index=False)

    print("\nHybrid cells — P(specialist), P(Arc 1 corner):")
    hyb = by_program[by_program["program"] == "hybrid"]
    print(
        hyb[
            [
                "dataset",
                "n",
                "p_specialist",
                "p_interior",
                "p_corner_arc1",
                "p_corner_arc2",
                "p_corner_arc3",
                "median_max_weight",
            ]
        ].to_string(index=False)
    )
    print("\nP(program | region):")
    print(given_region.to_string(index=False))
    if len(patients):
        top = patients.iloc[0]
        print(
            f"\nPal hybrid specialists at Arc 1: n={int(patients['n'].sum())}; "
            f"largest patient {top['patient']} ({top['subtype']}) "
            f"{top['frac_of_hybrid_arc1']:.1%}"
        )

    plot(by_program, given_region, FIG / "Figure_specialist_vs_interior_hybrid.png")

    summary = {
        "geometry_rule": (
            "specialist if (largest archetype weight − second) >= 0.35 "
            "after clipping negative barycentric weights to 0 and renormalizing; "
            "else interior"
        ),
        "program_rule": (
            "within each dataset: hybrid if epi and mes KS scores both >= 75th percentile; "
            "else epithelial if EMT <= 33rd percentile; mesenchymal if EMT >= 67th percentile"
        ),
        "pole_names": {"1": "hybrid", "2": "epithelial", "3": "mesenchymal"},
        "hybrid_placement": hyb.to_dict(orient="records"),
        "pal_program_note": (
            "Program labels were recomputed from MAGIC KS scores with the same quantiles "
            "as script 24. Agreement with the saved Pal labels is printed at runtime."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("\nWrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
