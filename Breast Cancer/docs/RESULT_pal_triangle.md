# Pal (GSE161529) — where cells sit, and corner versus middle

Host dataset for the k=3 triangle. 84,602 cells. Archetype weights are the Pal PCHA weights. Corners were named from the Tan 2014 KS genes on the cells nearest each vertex: Arc 2 epithelial, Arc 3 mesenchymal, Arc 1 hybrid (both gene sets high).

Two labels are used, and they are not the same thing.

**Nearest corner.** The archetype with the largest weight.

**Gene score.** Computed inside this dataset. Hybrid means both the epithelial score and the mesenchymal score are at or above the 75th percentile. Among the remaining cells, epithelial means the mesenchymal-minus-epithelial score is at or below the 33rd percentile, and mesenchymal means it is at or above the 67th percentile. The rest are intermediate.

**Corner versus middle.** A cell is a specialist (at a corner) when its largest weight beats the second-largest by at least 0.35. Otherwise it is interior (in the middle, including cells on an edge). This is the rule behind the 66% and 53% figures. In the comparison note, interior is the generalist class.

Script: `Breast Cancer/codes/26_specialist_vs_interior.py`. Per-cell file: `Breast Cancer/results/specialist_interior/scores_GSE161529.csv`.

Epithelial and mesenchymal scores in Pal are positively correlated (r = 0.29), so cells with both scores high are common: 8,719 hybrid cells, 10.3% of the dataset.

## 1. Where each group sits

The three “nearest” columns add to 100%. The three gene-score columns do not, because some cells are intermediate.

| Group | Cells | Nearest hybrid (Arc 1) | Nearest epithelial (Arc 2) | Nearest mesenchymal (Arc 3) | Gene score hybrid | Gene score epithelial | Gene score mesenchymal |
|---|---:|---:|---:|---:|---:|---:|---:|
| ER+ | 36,857 | 0% | 86.3% | 13.7% | 5.0% | 31.3% | 37.3% |
| HER2+ | 17,557 | 0% | 99.7% | 0.2% | 2.3% | 72.2% | 5.9% |
| TNBC | 30,188 | 66.1% | 1.3% | 32.6% | 21.5% | 5.8% | 36.6% |
| TN-0135 | 14,322 | 100% | 0% | 0% | 29.8% | 11.8% | 4.2% |
| TN-B1-0131 | 5,491 | 0% | 0% | 100% | 0% | 0% | 95.5% |
| TN-0126 | 1,223 | 8.3% | 0.7% | 91.0% | 21.8% | 0% | 77.8% |
| TN-B1-4031 | 5,128 | 39.5% | 0.7% | 59.8% | 12.4% | 0.3% | 32.9% |
| TN-B1-0554 | 2,302 | 96.7% | 2.6% | 0.8% | 20.3% | 0.4% | 78.4% |

ER+ and HER2+ sit on the epithelial corner. HER2+ has no corner of its own. TNBC uses two corners, and that split is mostly between patients.

TN-0135 (14,322 cells) is entirely nearest Arc 1, and every one of those cells passes the corner cut. Median weights are 0.88, 0.01, 0.11. Of its cells, 29.8% are gene-score hybrid. This is the hybrid specialist.

TN-B1-0131 (5,491 cells) is entirely nearest Arc 3. 95.5% are gene-score mesenchymal and none are hybrid. This is the mesenchymal TNBC patient.

TN-B1-4031 contains both (39.5% nearest Arc 1, 59.8% nearest Arc 3).

TN-B1-0554 is a warning about the two labels. 96.7% of its cells are nearest the hybrid corner, but 78.4% are called mesenchymal by the gene score. Nearest Arc 1 does not mean the cell turned both gene sets on.

## 2. Corner or middle

Given the gene score, where the cell sits:

| Gene score | Cells | At a corner | In the middle |
|---|---:|---:|---:|
| Epithelial | 25,964 | 69.3% | 30.7% |
| Mesenchymal | 25,868 | 73.0% | 27.0% |
| Hybrid | 8,719 | 65.6% | 34.4% |
| Intermediate | 24,051 | 81.5% | 18.5% |

Given the position, what the gene score is:

| Position | Cells | Epithelial | Mesenchymal | Hybrid | Intermediate |
|---|---:|---:|---:|---:|---:|
| At a corner | 62,200 | 28.9% | 30.4% | 9.2% | 31.5% |
| In the middle | 22,402 | 35.6% | 31.1% | 13.4% | 19.9% |

Of the 8,719 hybrid cells, 65.6% are specialists. 53.2% are specialists at Arc 1, 10.3% at Arc 2, and 2.1% at Arc 3. Median weights of hybrid cells are 0.68, 0.02, 0.11. Median purity is 0.52, so the leading weight is well clear of the second.

Of the 4,635 hybrid cells that are specialists at Arc 1, 92.2% are TN-0135. The hybrid corner in Pal is that patient. TN-0135 hybrid cells have median weights 0.91, 0.00, 0.08, and all of them pass the corner cut.

A hybrid cell is about as likely to sit at a corner as an epithelial cell (69.3%) or a mesenchymal cell (73.0%). A cell in the middle is hybrid only 13.4% of the time. It is epithelial 35.6% of the time and mesenchymal 31.1% of the time.

## Reading

On Pal, hybrid is a third corner, not the middle of the triangle. That corner is one TNBC patient. Other TNBC patients sit at the mesenchymal corner. The interior hybrid cells (34.4% of hybrids) are a real subset, but epithelial and mesenchymal cells are interior at similar rates, so the middle is not a hybrid compartment.
