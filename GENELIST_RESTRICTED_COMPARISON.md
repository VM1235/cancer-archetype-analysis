# Gene-list-restricted archetype analysis vs full transcriptome

**Headline: restricting PCA/PCHA to curated subtype markers did not produce a significant simplex for either cancer.** No gene-list \(k\) has \(p < 0.05\). The full-transcriptome runs were also non-significant. The previous nulls are **not** explained by “irrelevant genes dominating PCA,” under this protocol.

Official full-transcriptome outputs in `results/panel_a/`, `panel_b/`, and `panel_c_tcga/` were **not** overwritten. This run lives in `panel_a_genelist/`, `panel_b_genelist/`, `panel_c_genelist/` (plus `*_genelist.png` figures).

---

## Gene lists (Step 1)

### Breast — PAM50 (Parker / genefu)

- **Source:** `genefu` 2.42.0, `data(pam50)$centroids.map` (identical 50 probes to `pam50.robust`, the object used for labeling).
- **Citation:** Parker et al., *J Clin Oncol* 2009; centroids shipped with genefu.
- **File:** `Breast Cancer/data/processed/pam50_genefu_centroids_map.csv`
- **n:** 50 unique symbols. Probe column *is* the gene symbol (`ACTR3B`, `ESR1`, …).
- **In the official 16,500 × 63 matrix:** 47 exact matches; 3 HGNC aliases (`CDCA1→NUF2`, `KNTC2→NDC80`, `ORC6L→ORC6`); **0 missing**.
- **Restricted matrix:** **50 genes × 63 lines**  
  `Breast Cancer/data/processed/input_panelA_pam50_genelist.csv`

### GBM — Wang et al. 2017 (not the 7-gene handmade set)

- **Source:** Wang et al., *Cancer Cell* **32**, 42–56 (2017), Table S1 sheet **Subtype Signatures** (`mmc2.xlsx` from the article). doi:[10.1016/j.ccell.2017.06.003](https://doi.org/10.1016/j.ccell.2017.06.003).
- **Copy in repo:** `Glioblastoma/data/processed/Wang2017_CancerCell_TableS1_mmc2.xlsx` and `wang2017_tableS1_signatures.csv`.
- **What the table is:** 50 **upregulated** genes each for Mesenchymal, Proneural, and Classical (Neural is not tumor-intrinsic in that paper). Headers on this mmc2 file are the **corrected** MES/PN labels (the original supplement had those two headers swapped; this file matches the published correction).
- **Counts:** MES 50, PN 50, CL 50; **150 unique** after union (no cross-subtype duplicates).
- **In the official 15,833 × 54 matrix (explicit missing list — not silent):**
  - 135 present under the paper symbol
  - 8 present only after alias: `ACPP→ACP3`, `HN1→JPT1`, `HRASLS→PLAAT1`, `ZNF643→ZSCAN31`, `LHFP→LHFPL6`, `SEPT11→SEPTIN11`, `C14orf159→DGLUCY`, `KIAA0494→EFCAB14`
  - **7 missing (no row in DepMap matrix):**  
    Proneural: `NPPA`, `PAK7`, `PCDHA9`, `SLC17A6`, `NHLH1`  
    Classical: `ARNTL`, `GRIK1`  
    (`PAK7`/`PAK5` both absent.)
- **Restricted matrix:** **143 genes × 54 lines**  
  `Glioblastoma/data/processed/input_panelA_wang2017_genelist.csv`

No extra log and no low-expression refilter: both restricted matrices are **row subsets** of the official Panel A CSVs.

---

## Significance flag

| Cancer | Full transcriptome: any \(p<0.05\)? | Gene-list: any \(p<0.05\)? | Changed conclusion? |
|---|---|---|---|
| Breast | No (best \(p=0.082\) at \(k=4\)) | No (best \(p=0.762\) at \(k=3\) and \(k=4\)) | **No** — still NS; \(p\) **increased** |
| GBM | No (smallest \(p=0.068\) at \(k=7\)) | No (smallest \(p=0.066\) at \(k=6\)) | **No** — still NS |

---

## Breast

| Metric | Full transcriptome (existing) | Gene-list restricted (new) |
|---|---|---|
| n_genes used for PCA/PCHA | 16,500 | 50 (PAM50) |
| n_PCs for ≥50% variance | 8 (50.0%) | **3 (59.2%)** |
| PCHA `numIter` (official for this cancer) | 5 (15/5 inits) | 5 (same) |
| k grid actually fit | 3–7 | 3–4 (k≥5 needs >3 PCs) |
| Suggested / “best” k | 4 (DimensionFinder; no \(p<0.05\)) | 3 (DimensionFinder; no \(p<0.05\)) |
| t-ratio at that k | 0.354 | 0.635 |
| p at that k (500 shuffles) | 0.082 | 0.762 |
| Subtype–archetype match count (bin-0 peak, q<0.1) | **2/3** (LumB, Basal; Her2 unmatched) | **2/3** (LumB, Basal; Her2 unmatched) |
| Tumor containment | **0/1097** | **526/1097** |

Breast Panel B still uses genefu PAM50 labels (LumA/Normal dropped from the test only). Vertex **indices** are not comparable across k=4 vs k=3.

**Containment caveat:** 526/1097 is a 2-D triangle (\(k=3\)). Points fall inside a low-D simplex much more easily than inside the official 3-D tetrahedron. Combined with \(p=0.762\), this is **not** a rescued Pareto claim.

Panel C ComBat: 49 shared genes (Xena HiSeqV2 has no `ORC6`/`ORC6L` in the prepared tumor table).

---

## GBM

| Metric | Full transcriptome (existing) | Gene-list restricted (new) |
|---|---|---|
| n_genes used for PCA/PCHA | 15,833 | 143 (Wang 2017 CL/MES/PN; 7 paper genes absent) |
| n_PCs for ≥50% variance | 12 (51.4%) | **7 (53.6%)** |
| PCHA `numIter` (official for this cancer) | 50 (150/50 inits) | 50 (same) |
| k grid actually fit | 3–7 | 3–7 |
| Suggested / “best” k | 7 (DimensionFinder; no \(p<0.05\)) | 5 (DimensionFinder; no \(p<0.05\)) |
| t-ratio at that k | 0.027 | 0.145 |
| p at that k (500 shuffles) | 0.068 | 0.104 |
| Smallest p in the k grid | 0.068 (\(k=7\)) | 0.066 (\(k=6\), t=0.060) |
| Subtype–archetype match count | **0/3** (CL, MES, PN marker labels) | **0/3** |
| Tumor containment | **0/154** | **2/154** |

Panel B still uses the existing cell-line marker z-score labels (`verhaak_labels_panelA.csv`), not a new ssGSEA. Labels are independent of PCHA either way.

Panel C ComBat: 135 genes. The 8 alias-renamed DepMap symbols (`JPT1`, `PLAAT1`, …) are absent from the already-prepared Xena table (that table was intersected with **full** Panel A symbols, which also use the modern names). Re-extracting HiSeq under 2017 symbols was not done.

---

## How to rerun (does not touch official folders)

From repo root, `.venv`, NumPy 1.x:

```bash
.venv/bin/python -u "Breast Cancer/codes/00_build_genelist_input.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelA_genelist.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelB_genelist.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelC_genelist.py" --n-shuffle 20

.venv/bin/python -u "Glioblastoma/codes/00_build_genelist_input.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelA_genelist.py"   # slow; numIter=50
.venv/bin/python -u "Glioblastoma/codes/run_panelB_genelist.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelC_genelist.py" --n-shuffle 20
```

---

## Reading for the thesis

The advisor’s hypothesis was: the cloud looks round because housekeeping/irrelevant genes dominate the leading PCs. After restricting to the **same gene sets that define the subtypes**, breast PAM50 space is even **less** simplex-like by t-ratio permutation, and GBM Wang-2017 space remains non-significant (closest \(p=0.066\)). A 1-to-1 vertex–subtype map still fails for GBM; breast still maps LumB and Basal only. Tumor-in-simplex counts should not be compared across different \(k\) and dimension without that caveat.
