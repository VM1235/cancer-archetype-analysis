#!/usr/bin/env python3
"""Robust ENA FASTQ download + Salmon quantification for one study.

Resumes partial downloads. Deletes FASTQs after each sample.
Usage: python3 codes/04b_quant_ena_resumable.py GSE102511
"""
from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes")
MAP = ROOT / "data/sra/gsm_srr_fastq_map.csv"
INDEX = ROOT / "data/ref/salmon_index"
TX2GENE = ROOT / "data/ref/tx2gene.tsv"
THREADS = "4"


def log(msg: str, logf: Path) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with logf.open("a") as f:
        f.write(line + "\n")


def download(url: str, dest: Path, logf: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not url.startswith("http"):
        url = "https://" + url
    # curl resume; retries
    cmd = [
        "curl",
        "-fL",
        "--retry",
        "20",
        "--retry-delay",
        "5",
        "--retry-all-errors",
        "-C",
        "-",
        "-o",
        str(dest),
        url,
    ]
    log(f"DOWNLOAD {url} -> {dest.name}", logf)
    subprocess.run(cmd, check=True)


def aggregate_gene_counts(outdir: Path, out_csv: Path) -> None:
    import pandas as pd

    tx2gene = pd.read_csv(TX2GENE, sep="\t")
    tx_map = dict(zip(tx2gene["TXNAME"].astype(str), tx2gene["GENEID"].astype(str)))
    samples = sorted([p for p in outdir.iterdir() if (p / "quant.sf").exists()])
    cols = {}
    for sdir in samples:
        q = pd.read_csv(sdir / "quant.sf", sep="\t")
        q["TX"] = q["Name"].astype(str).str.replace(r"\..*$", "", regex=True)
        q["GENE"] = q["TX"].map(tx_map)
        q = q.dropna(subset=["GENE"])
        cols[sdir.name] = q.groupby("GENE")["NumReads"].sum()
    mat = pd.DataFrame(cols).fillna(0.0)
    mat.index.name = "gene"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    mat.to_csv(out_csv)
    print(f"wrote {out_csv} {mat.shape}")


def main() -> None:
    study = sys.argv[1]
    work = ROOT / f"data/sra/work_{study}"
    outdir = ROOT / f"data/sra/salmon_{study}"
    logf = ROOT / f"data/sra/quant_{study}.log"
    work.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    rows = [r for r in csv.DictReader(MAP.open()) if r["study"] == study]
    log(f"START {study} n={len(rows)}", logf)

    for r in rows:
        sid, layout = r["sample_id"], r["layout"]
        qdir = outdir / sid
        if (qdir / "quant.sf").exists():
            log(f"SKIP {sid}", logf)
            continue

        sdir = work / sid
        sdir.mkdir(parents=True, exist_ok=True)

        urls = [u for u in r["fastq_ftp"].split(";") if u]
        fqs = []
        try:
            for i, u in enumerate(urls, 1):
                dest = sdir / f"r{i}.fastq.gz"
                download(u, dest, logf)  # curl -C - resumes partial
                fqs.append(dest)

            qdir.mkdir(parents=True, exist_ok=True)
            if layout == "PAIRED":
                cmd = [
                    "salmon",
                    "quant",
                    "-i",
                    str(INDEX),
                    "-l",
                    "A",
                    "-1",
                    str(fqs[0]),
                    "-2",
                    str(fqs[1]),
                    "-p",
                    THREADS,
                    "-o",
                    str(qdir),
                    "--validateMappings",
                ]
            else:
                cmd = [
                    "salmon",
                    "quant",
                    "-i",
                    str(INDEX),
                    "-l",
                    "A",
                    "-r",
                    str(fqs[0]),
                    "-p",
                    THREADS,
                    "-o",
                    str(qdir),
                    "--validateMappings",
                ]
            log(f"SALMON {sid}", logf)
            subprocess.run(cmd, check=True)
            log(f"DONE {sid}", logf)
        except Exception as e:
            log(f"FAIL {sid}: {e}", logf)
        finally:
            shutil.rmtree(sdir, ignore_errors=True)

    out_csv = ROOT / "data/processed" / f"counts_{study}.csv"
    n_done = sum(1 for p in outdir.iterdir() if (p / "quant.sf").exists())
    if n_done:
        aggregate_gene_counts(outdir, out_csv)
    log(f"COMPLETE {study} done={n_done}/{len(rows)}", logf)


if __name__ == "__main__":
    main()
