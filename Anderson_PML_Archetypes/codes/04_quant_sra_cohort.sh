#!/usr/bin/env bash
# Quantify one PML cohort: NCBI prefetch -> fasterq-dump -> Salmon -> gene counts.
# Deletes intermediate SRA/FASTQ after each sample to save disk.
# Usage: bash codes/04_quant_sra_cohort.sh GSE102511
set -euo pipefail
export KMP_DUPLICATE_LIB_OK=TRUE
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

STUDY="${1:?study accession e.g. GSE102511}"
ROOT="/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes"
MAP="$ROOT/data/sra/gsm_srr_fastq_map.csv"
INDEX="$ROOT/data/ref/salmon_index"
TX2GENE="$ROOT/data/ref/tx2gene.tsv"
WORKDIR="$ROOT/data/sra/work_${STUDY}"
OUTDIR="$ROOT/data/sra/salmon_${STUDY}"
LOG="$ROOT/data/sra/quant_${STUDY}.log"
mkdir -p "$WORKDIR" "$OUTDIR"

echo "[$(date)] START $STUDY" | tee "$LOG"

python3 - "$MAP" "$STUDY" "$WORKDIR/samples.tsv" <<'PY'
import csv, sys
inp, study, out = sys.argv[1:4]
rows=[r for r in csv.DictReader(open(inp)) if r["study"]==study]
with open(out,"w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=["sample_id","srr","layout"], delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow({k:r[k] for k in w.fieldnames})
print(len(rows), "samples", file=sys.stderr)
PY

while IFS=$'\t' read -r sample_id srr layout || [[ -n "${sample_id:-}" ]]; do
  sample_id="${sample_id//$'\r'/}"
  srr="${srr//$'\r'/}"
  layout="${layout//$'\r'/}"
  [[ "$sample_id" == "sample_id" || -z "$sample_id" ]] && continue

  qdir="$OUTDIR/$sample_id"
  if [[ -f "$qdir/quant.sf" ]]; then
    echo "[$(date)] SKIP $sample_id" | tee -a "$LOG"
    continue
  fi

  echo "[$(date)] QUANT $sample_id $srr $layout" | tee -a "$LOG"
  sdir="$WORKDIR/$sample_id"
  rm -rf "$sdir"
  mkdir -p "$sdir"
  cd "$sdir"

  # 1) Download SRA via NCBI (HTTPS / cloud)
  prefetch "$srr" -O "$sdir" -v 2>>"$LOG"

  # 2) Convert to FASTQ
  if [[ "$layout" == "PAIRED" ]]; then
    fasterq-dump "$srr" -O "$sdir" -e "$OMP_NUM_THREADS" --split-files 2>>"$LOG"
    # compress for salmon streaming friendliness / disk
    gzip -f "${srr}_1.fastq" "${srr}_2.fastq"
    fq1="${sdir}/${srr}_1.fastq.gz"
    fq2="${sdir}/${srr}_2.fastq.gz"
  else
    fasterq-dump "$srr" -O "$sdir" -e "$OMP_NUM_THREADS" 2>>"$LOG"
    gzip -f "${srr}.fastq"
    fq1="${sdir}/${srr}.fastq.gz"
  fi

  # free SRA container
  rm -rf "${sdir}/${srr}" "${sdir}/${srr}.sra" 2>/dev/null || true

  mkdir -p "$qdir"
  if [[ "$layout" == "PAIRED" ]]; then
    salmon quant -i "$INDEX" -l A -1 "$fq1" -2 "$fq2" \
      -p "$OMP_NUM_THREADS" -o "$qdir" --validateMappings 2>>"$LOG"
  else
    salmon quant -i "$INDEX" -l A -r "$fq1" \
      -p "$OMP_NUM_THREADS" -o "$qdir" --validateMappings 2>>"$LOG"
  fi

  rm -rf "$sdir"
  echo "[$(date)] DONE $sample_id" | tee -a "$LOG"
done < "$WORKDIR/samples.tsv"

python3 - "$OUTDIR" "$TX2GENE" "$ROOT/data/processed/counts_${STUDY}.csv" <<'PY'
import sys
from pathlib import Path
import pandas as pd
outdir, tx2gene_path, out_csv = map(Path, sys.argv[1:4])
tx2gene = pd.read_csv(tx2gene_path, sep="\t")
tx_map = dict(zip(tx2gene["TXNAME"].astype(str), tx2gene["GENEID"].astype(str)))
samples = sorted([p for p in outdir.iterdir() if (p/"quant.sf").exists()])
gene_cols = {}
for sdir in samples:
    q = pd.read_csv(sdir/"quant.sf", sep="\t")
    q["TX"] = q["Name"].astype(str).str.replace(r"\..*$", "", regex=True)
    q["GENE"] = q["TX"].map(tx_map)
    q = q.dropna(subset=["GENE"])
    gene_cols[sdir.name] = q.groupby("GENE")["NumReads"].sum()
mat = pd.DataFrame(gene_cols).fillna(0.0)
mat.index.name = "gene"
mat.to_csv(out_csv)
print(f"wrote {out_csv} genes={mat.shape[0]} samples={mat.shape[1]}")
PY

echo "[$(date)] COMPLETE $STUDY" | tee -a "$LOG"
