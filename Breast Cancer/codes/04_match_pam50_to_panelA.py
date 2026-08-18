#!/usr/bin/env python3
"""Step 0d: match genefu PAM50 output to the exact 63 Panel A sample IDs.

Does not assume row order. Matches by ModelID first, then stripped aliases
from input_panelA_models_used.csv.
"""

from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
OUT = BREAST / "results" / "panel_b"

MATRIX = BREAST / "data" / "processed" / "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv"
MODELS = BREAST / "data" / "processed" / "input_panelA_models_used.csv"
PAM50_RAW = OUT / "pam50_genefu_raw.csv"
OUT_LABELS = OUT / "pam50_labels_panelA.csv"
OUT_REPORT = OUT / "step0d_name_match_report.txt"


def norm_name(x):
    s = str(x).strip().upper()
    s = s.replace("_BREAST", "")
    s = re.sub(r"[^A-Z0-9]", "", s)
    return s


def main():
    import pandas as pd

    expr = pd.read_csv(MATRIX, index_col=0)
    panel_ids = list(expr.columns.astype(str))
    print(f"Panel A sample IDs: {len(panel_ids)}")

    raw = pd.read_csv(PAM50_RAW)
    raw["cell_line_genefu"] = raw["cell_line_genefu"].astype(str)
    print(f"genefu rows: {len(raw)}")
    print("genefu ID examples:", list(raw["cell_line_genefu"].head(5)))

    models = pd.read_csv(MODELS)
    id_col = "ModelID" if "ModelID" in models.columns else models.columns[0]
    alias_cols = [
        c
        for c in (
            "ModelID",
            "CellLineName",
            "StrippedCellLineName",
            "CCLEName",
            "RRID",
        )
        if c in models.columns
    ]

    alias_to_id = {}
    for _, row in models.iterrows():
        mid = str(row[id_col])
        for c in alias_cols:
            alias_to_id[str(row[c])] = mid
            alias_to_id[norm_name(row[c])] = mid
        alias_to_id[mid] = mid
        alias_to_id[norm_name(mid)] = mid

    matched = {}
    unmatched_genefu = []
    for _, row in raw.iterrows():
        gid = str(row["cell_line_genefu"])
        mid = None
        if gid in panel_ids:
            mid = gid
        elif gid in alias_to_id and alias_to_id[gid] in panel_ids:
            mid = alias_to_id[gid]
        elif norm_name(gid) in alias_to_id and alias_to_id[norm_name(gid)] in panel_ids:
            mid = alias_to_id[norm_name(gid)]
        if mid is None:
            unmatched_genefu.append(gid)
            continue
        matched[mid] = row

    missing_panel = [s for s in panel_ids if s not in matched]
    n_ok = len(panel_ids) - len(missing_panel)
    coverage = n_ok / len(panel_ids)

    print(f"\nMatched {n_ok}/{len(panel_ids)} Panel A lines ({100 * coverage:.1f}%)")
    if missing_panel:
        print("Panel A IDs with no PAM50 match:")
        for s in missing_panel:
            print(" ", s)
    else:
        print("All 63 Panel A IDs matched a PAM50 row.")
    if unmatched_genefu:
        print("genefu IDs that did not map onto the Panel A set:")
        for s in unmatched_genefu:
            print(" ", s)

    rows = []
    for sid in panel_ids:
        if sid not in matched:
            rows.append(
                {
                    "cell_line": sid,
                    "pam50_subtype": pd.NA,
                    "confidence_score": pd.NA,
                }
            )
            continue
        r = matched[sid]
        rows.append(
            {
                "cell_line": sid,
                "pam50_subtype": r["pam50_subtype"],
                "confidence_score": r["confidence_score"],
            }
        )
    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_LABELS, index=False)
    report = [
        f"n_panelA={len(panel_ids)}",
        f"n_matched={n_ok}",
        f"coverage={coverage:.4f}",
        f"unmatched_panelA={','.join(missing_panel) if missing_panel else ''}",
        f"unmatched_genefu={','.join(unmatched_genefu) if unmatched_genefu else ''}",
        f"output={OUT_LABELS}",
        "",
    ]
    OUT_REPORT.write_text("\n".join(report))
    print("Wrote", OUT_LABELS)
    print(out["pam50_subtype"].value_counts(dropna=False).to_string())

    if coverage < 0.80:
        print(
            "\nSTOP: match coverage < 80%. Not running Panel B. "
            "This is likely a naming-convention bug, not unclassifiable lines."
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
