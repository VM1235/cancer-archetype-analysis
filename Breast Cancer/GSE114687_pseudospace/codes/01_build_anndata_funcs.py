#!/usr/bin/env python3
"""Step 1: read GSE114687 pseudospace data (McFaline-Figueroa et al.), QC, normalise, save.

Two human mammary epithelial systems, each grown as a colony so that EMT progresses with distance from the
colony edge ("pseudospace"), with and without TGF-beta:
  MCF10A  (immortalised)          matrix 'pseudospace'   9,268 cells kept by the authors (of 23,951 droplets)
  HuMEC   (primary)               matrix 'HuMEC'        22,352 cells
Arms: Mock (spontaneous EMT) and TGFB (4 ng/mL). Position: inner or outer colony.
The CROP-seq perturbation screen (also in GSE114687) is not used.

Only cells present in the authors' metadata are kept (the matrix also holds low-quality droplets).
Labels carried per cell:
  line      = MCF10A / HuMEC          (also stored as `assay`)
  arm       = Mock / TGFB
  position  = inner / outer           (ground-truth side of the EMT axis, like time in Cook)
  sample    = <line>_<arm>            (batch unit, the analogue of the tumour / library)
Light QC (as for 4T1): >= MIN_GENES genes, mito fraction <= MAX_MT.  Genes detected in >= 1% of kept cells; mito
and ribosomal genes dropped. Expression = log1p of counts per 10,000. Output is named expr_tumour.h5ad so the shared
code finds it, although these are cultured cells, not tumours.
"""
import gzip
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, scanpy as sc

HERE = Path(__file__).resolve().parent.parent
RAW = HERE / "data" / "raw"
OUT = HERE / "data" / "processed"
MIN_GENES, MAX_MT = 500, 15.0
S_GENES = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RAP2C RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8".split()
G2M_GENES = "HMGB2 DDX11 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPK TACC3 SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA".split()


def read_mtx(path):
    """Fast MatrixMarket reader: returns (genes x cells) COO parts."""
    with gzip.open(path, "rt") as f:
        l = f.readline()
        while l.startswith("%"):
            l = f.readline()
        ng, nc, nnz = map(int, l.split())
    df = pd.read_csv(path, sep=" ", header=None, skiprows=3, names=["i", "j", "x"],
                     dtype={"i": np.int32, "j": np.int32, "x": np.float32}, compression="gzip")
    assert len(df) == nnz, (len(df), nnz)
    return df, ng, nc


def load_line(line, prefix, meta_file, header):
    genes = pd.read_csv(RAW / f"GSE114687_{prefix}_genes.tsv.gz", sep="\t", header=None, names=["ens", "sym"])
    bcs = pd.read_csv(RAW / f"GSE114687_{prefix}_barcodes.tsv.gz", header=None)[0].astype(str).values
    m = pd.read_csv(RAW / meta_file, sep="\t", header=0 if header else None)
    if not header:
        m.columns = ["cell", "total_umis", "sample", "TSNE.1", "TSNE.2", "Size_Factor", "treatment_id", "spatial_id"]
    m["barcode"] = m["cell"].astype(str).str.split("_").str[0]
    pos = {b: i for i, b in enumerate(bcs)}
    col = m["barcode"].map(pos)
    print(f"{line}: {len(m)} metadata cells, {col.isna().sum()} not found among {len(bcs)} matrix barcodes", flush=True)
    m = m[col.notna()].copy()
    col = col[col.notna()].astype(int).values
    df, ng, nc = read_mtx(RAW / f"GSE114687_{prefix}_matrix.mtx.gz")
    new = np.full(nc, -1, dtype=np.int64)
    new[col] = np.arange(len(col))
    keep = new[df["j"].values - 1] >= 0
    X = sp.coo_matrix((df["x"].values[keep], (new[df["j"].values[keep] - 1], df["i"].values[keep] - 1)),
                      shape=(len(col), ng)).tocsr()
    del df
    sym = genes["sym"].astype(str).values
    dup = pd.Index(sym).duplicated(keep=False)
    sym = np.array([f"{s}__{e}" if d else s for s, e, d in zip(sym, genes["ens"], dup)])
    mt = np.array([s.startswith("MT-") for s in sym])
    n = np.asarray(X.sum(1)).ravel(); ngn = np.diff(X.indptr)
    pm = 100 * np.asarray(X[:, np.where(mt)[0]].sum(1)).ravel() / np.maximum(n, 1)
    obs = pd.DataFrame({"line": line, "assay": line, "arm": m["treatment_id"].values, "position": m["spatial_id"].values,
                        "sample": [f"{line}_{a}" for a in m["treatment_id"].values],
                        "nCount_RNA": n, "nGene": ngn, "pct_counts_mt": pm}, index=[f"{line}_{c}" for c in m["cell"].values])
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=sym))


