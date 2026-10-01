#!/usr/bin/env python3
"""Assemble a combined genes x cells raw count matrix for GSE161529.

Unlike GSE176078 (one pre-merged count_matrix_sparse.mtx for the whole
atlas), GSE161529's GEO deposit is Cell Ranger's raw per-sample output:
one barcodes.tsv.gz + matrix.mtx.gz pair per GSM (69 samples), sharing a
single features.tsv gene reference across all of them (the standalone
features.tsv you already downloaded).

This script loops over the 27 samples actually used by the companion
TNBC/HER2/ER analyses (github.com/yunshun/HumanBreast10X), applies the
same per-sample QC thresholds Chen et al. published in SampleStats.txt,
and concatenates into one matrix -- playing the same role
count_matrix_sparse.mtx / count_matrix_genes.tsv / count_matrix_barcodes.tsv
played for GSE176078.

BEFORE RUNNING:
  1. Extract GSE161529_RAW.tar somewhere (RAW_DIR below).
  2. Run 00_extract_metadata_gse161529.R first, so META_CSV exists.
  3. Matching is keyed on GSM ID (from SampleStats.txt), not on the raw
     filename's sample-name text, since GEO's on-disk names carry extra
     site-code prefixes the clean names don't. The __main__ block still
     sanity-checks the first sample before processing all 27.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/00_build_matrix_gse161529.py"
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
BASE = BREAST / "data" / "gse161529"
RAW_DIR = BASE / "GSE161529_RAW"          # you'll create this by extracting the tar (see chat)
FEATURES_PATH = BASE / "features.tsv"
SAMPLE_STATS = BASE / "SampleStats_GSE161529.txt"
META_CSV = BASE / "metadata_gse161529.csv"  # written by the R script

OUT_DIR = BASE / "processed"

# The 27 samples used across TNBC.R + HER2.R + ER.R in
# github.com/yunshun/HumanBreast10X -- the same set the Figshare
# SeuratObject_*Tum.rds files were built from.
TNBC_SAMPLES = ["TN-B1-0554", "TN-B1-0177", "TN-0135", "TN-B1-4031",
                "TN-B1-0131", "TN-0126", "TN-0106", "TN-0114-T2"]
HER2_SAMPLES = ["HER2-0031", "HER2-0337", "HER2-0308", "HER2-0069",
                "HER2-0161", "HER2-0176"]
ER_SAMPLES = ["ER-0114-T3", "ER-0360", "ER-0167-T", "ER-0151", "ER-0032",
              "ER-0125", "ER-0043-T", "ER-0025", "ER-0001", "ER-0042",
              "ER-0319", "ER-0040-T", "ER-0163"]
ALL_SAMPLES = TNBC_SAMPLES + HER2_SAMPLES + ER_SAMPLES


def read_features(path: Path) -> list[str]:
    """Return gene symbols: col 2 (0-indexed 1) if present, else col 1."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        rows = [l.rstrip("\n").split("\t") for l in f]
    col = 1 if len(rows[0]) > 1 else 0
    return [r[col] for r in rows]


def find_sample_files(raw_dir: Path, gsm_id: str) -> dict[str, Path] | None:
    """Locate the barcodes/matrix pair for one GEO sample, keyed by GSM ID.

    Actual GSE161529 filenames look like:
        GSM4909287_TN-B1-0554-barcodes.tsv.gz
        GSM4909287_TN-B1-0554-matrix.mtx.gz
    The "sample name" chunk in the middle sometimes carries extra site-code
    prefixes GEO doesn't expose anywhere else (e.g. normal samples show up
    as "N-PM0092-Total" on disk vs "N-0092-total" in SampleStats.txt), so
    matching on GSM ID -- which SampleStats.txt gives unambiguously -- is
    the reliable key, not the human-readable name.
    """
    cand = [f for f in raw_dir.glob(f"{gsm_id}_*") if f.is_file()]
    if not cand:
        return None
    hits = {}
    bc = [f for f in cand if "barcode" in f.name.lower()]
    mx = [f for f in cand if "matrix" in f.name.lower()]
    if bc:
        hits["barcodes"] = bc[0]
    if mx:
        hits["matrix"] = mx[0]
    return hits if len(hits) == 2 else None


def load_one_sample(files: dict[str, Path]) -> tuple[sp.csr_matrix, list[str]]:
    """Return (genes x cells sparse counts, barcodes) for one sample."""
    mopen = gzip.open if files["matrix"].suffix == ".gz" else open
    with mopen(files["matrix"], "rt") as f:
        mat = sio.mmread(f).tocsr()
    bopen = gzip.open if files["barcodes"].suffix == ".gz" else open
    with bopen(files["barcodes"], "rt") as f:
        barcodes = [l.strip() for l in f]
    return mat, barcodes


def main() -> int:
    if not FEATURES_PATH.is_file():
        print(f"Missing {FEATURES_PATH}")
        return 1
    genes_ref = read_features(FEATURES_PATH)
    n_genes_ref = len(genes_ref)
    print(f"Shared gene reference: {n_genes_ref} genes")

    if not SAMPLE_STATS.is_file():
        print(f"Missing {SAMPLE_STATS}")
        return 1
    stats = pd.read_csv(SAMPLE_STATS, sep="\t")
    stats.columns = [c.strip() for c in stats.columns]
    stats = stats.set_index("SampleName")

    # --- sanity check on the first sample before committing to all 27 ---
    probe_name = ALL_SAMPLES[0]
    probe_gsm = str(stats.loc[probe_name, "GEO_ID"]).strip()
    probe = find_sample_files(RAW_DIR, probe_gsm)
    if probe is None:
        print(f"!! Could not locate files for {probe_name!r} (GSM {probe_gsm}) under {RAW_DIR}.")
        print("   List what's actually there and adjust find_sample_files():")
        print("   ", [p.name for p in list(RAW_DIR.glob('*'))[:20]])
        return 1
    print(f"Found for {probe_name} ({probe_gsm}): {probe}")

    mats, all_barcodes, sample_of = [], [], []
    for name in ALL_SAMPLES:
        gsm_id = str(stats.loc[name, "GEO_ID"]).strip()
        files = find_sample_files(RAW_DIR, gsm_id)
        if files is None:
            print(f"!! no files found for {name} (GSM {gsm_id}), skipping")
            continue
        mat, barcodes = load_one_sample(files)
        if mat.shape[0] != n_genes_ref:
            print(f"!! {name}: {mat.shape[0]} genes != shared reference {n_genes_ref}; skipping")
            continue

        row = stats.loc[name]
        lib = np.asarray(mat.sum(axis=0)).ravel()
        n_genes_detected = np.asarray((mat > 0).sum(axis=0)).ravel()
        mito_idx = [i for i, g in enumerate(genes_ref) if str(g).startswith("MT-")]
        mito_frac = (
            np.asarray(mat[mito_idx, :].sum(axis=0)).ravel() / np.maximum(lib, 1)
            if mito_idx else np.zeros(mat.shape[1])
        )
        keep = (
            (mito_frac < float(row["Mito"]))
            & (n_genes_detected > float(row["GeneLower"]))
            & (n_genes_detected < float(row["GeneUpper"]))
            & (lib < float(row["LibSize"]))
        )
        print(f"{name} (GSM {gsm_id}): {int(keep.sum())}/{mat.shape[1]} cells pass QC "
              f"(mito<{row['Mito']}, {row['GeneLower']}<genes<{row['GeneUpper']}, lib<{row['LibSize']})")

        mat_keep = mat[:, keep]
        sample_comb = name.replace("-", "_")
        bc_keep = [f"{sample_comb}_{barcodes[i]}" for i in np.where(keep)[0]]

        mats.append(mat_keep)
        all_barcodes.extend(bc_keep)
        sample_of.extend([name] * mat_keep.shape[1])

    if not mats:
        print("No samples loaded -- nothing to write.")
        return 1

    combined = sp.hstack(mats).tocsr()  # genes x cells
    print(f"Combined matrix: {combined.shape[0]} genes x {combined.shape[1]} cells "
          f"across {len(mats)} samples")

    if META_CSV.is_file():
        meta = pd.read_csv(META_CSV)
        have = set(meta["full_barcode"])
        matched = sum(bc in have for bc in all_barcodes)
        print(f"{matched}/{len(all_barcodes)} raw cells have a Figshare cell-type match "
              f"(run 00_extract_metadata_gse161529.R first if this is 0)")
    else:
        print(f"Note: {META_CSV} not found yet -- run the R script first to get cell-type labels.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sp.save_npz(OUT_DIR / "gse161529_combined_counts.npz", combined)
    pd.Series(genes_ref).to_csv(OUT_DIR / "gse161529_genes.csv", index=False, header=["gene"])
    pd.DataFrame({"barcode": all_barcodes, "sample": sample_of}).to_csv(
        OUT_DIR / "gse161529_barcodes.csv", index=False
    )
    print("Wrote", OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
