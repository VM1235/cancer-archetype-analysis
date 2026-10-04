# GSE173634 — same triangle, cell lines rather than patients

35,271 cells from breast cancer cell lines (Gambardella et al.), projected into the frozen Pal k=3 triangle. Weights and the corner cut are the same as for Wu. Gene-score percentiles are computed inside GSE173634.

Study labels: LA luminal A, LB luminal B, H HER2, TNA and TNB two triple-negative groups, Basal-like.

Script: `Breast Cancer/codes/26_specialist_vs_interior.py`. Per-cell file: `Breast Cancer/results/specialist_interior/scores_GSE173634.csv`.

Epithelial and mesenchymal scores are strongly anti-correlated (r = −0.77). A cell almost never has both scores in the top quartile. There are 65 hybrid cells, 0.2% of the dataset, and all 65 come from one line, CAL851 (TNB).

38.4% of cells are specialists (13,529 / 35,271). Every mesenchymal cell is interior (0 / 11,757 specialists).

## 1. Where each group sits

| Group | Cells | Nearest hybrid (Arc 1) | Nearest epithelial (Arc 2) | Nearest mesenchymal (Arc 3) | Gene score hybrid | Gene score epithelial | Gene score mesenchymal | At a corner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| TNA | 11,716 | 55.6% | 44.4% | 0% | 0% | 24.0% | 35.8% | 18.6% |
| LA | 7,565 | 0% | 100% | 0% | 0% | 60.9% | 0% | 88.8% |
| TNB | 6,988 | 5.8% | 18.2% | 76.0% | 0.9% | 0.3% | 85.7% | 0% |
| H | 4,347 | 50.3% | 49.7% | 0% | 0% | 32.4% | 8.1% | 39.4% |
| LB | 2,911 | 0.2% | 99.8% | 0% | 0% | 99.8% | 0.1% | 99.8% |
| Basal-like | 1,744 | 86.4% | 0% | 13.6% | 0% | 0% | 69.9% | 0.7% |

Luminal lines occupy the epithelial corner and pass the corner cut. LB is 99.8% epithelial by gene score and 99.8% specialists. LA is 100% nearest Arc 2 and 88.8% specialists, with no mesenchymal gene-score cells.

TNB is 76.0% nearest the mesenchymal corner and 85.7% mesenchymal by gene score, and 0% of TNB cells are specialists. Basal-like (this cohort’s Basal-like cells are the MCF12A line) is 86.4% nearest Arc 1 and 69.9% mesenchymal by gene score, and 99.3% interior. Nearest the hybrid corner, here, is a mesenchymal gene-score line that does not reach the vertex.

HER2 (H) is split evenly between Arc 1 and Arc 2, with no hybrid cells.

Lines with at least 800 cells:

| Line | Subtype | Cells | Nearest Arc 1 | Nearest Arc 2 | Nearest Arc 3 | Mesenchymal | Epithelial | At a corner |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| HCC70 | TNA | 2,818 | 100% | 0% | 0% | 0% | 62.4% | 75.5% |
| MDAMB436 | TNA | 2,817 | 0.1% | 99.9% | 0% | 100% | 0% | 0% |
| BT549 | TNB | 2,295 | 0% | 0% | 100% | 100% | 0% | 0% |
| BT474 | LB | 2,158 | 0.2% | 99.8% | 0% | 0.1% | 99.8% | 99.8% |
| MCF12A | Basal-like | 1,744 | 86.4% | 0% | 13.6% | 69.9% | 0% | 0.7% |
| HCC1954 | H | 1,623 | 100% | 0% | 0% | 0% | 17.1% | 0% |
| MDAMB468 | TNA | 1,573 | 15.6% | 84.4% | 0% | 0% | 3.6% | 0% |
| EFM19 | LA | 1,316 | 0% | 100% | 0% | 0% | 41.0% | 91.3% |
| HCC1143 | TNA | 1,280 | 99.9% | 0.1% | 0% | 99.8% | 0% | 0% |
| CAL51 | TNB | 1,260 | 0% | 100% | 0% | 100% | 0% | 0% |
| HCC1187 | TNA | 1,147 | 100% | 0% | 0% | 7.8% | 0% | 0% |
| MDAMB453 | H | 1,119 | 0% | 100% | 0% | 0% | 8.6% | 99.9% |
| HCC1937 | TNA | 1,028 | 99.0% | 1.0% | 0% | 0% | 97.2% | 4.8% |
| CAL851 | TNB | 995 | 33.9% | 0.2% | 65.9% | 0.1% | 2.0% | 0% |
| ZR751 | LA | 943 | 0% | 100% | 0% | 0% | 80.9% | 22.9% |
| KPL1 | LA | 942 | 0% | 100% | 0% | 0% | 91.1% | 100% |
| HS578T | TNB | 879 | 0% | 0% | 100% | 100% | 0% | 0% |
| MDAMB415 | LA | 877 | 0% | 100% | 0% | 0% | 81.4% | 100% |
| MCF7 | LA | 839 | 0% | 100% | 0% | 0% | 70.4% | 100% |
| T47D | LA | 825 | 0% | 100% | 0% | 0% | 91.5% | 99.6% |
| CAMA1 | LA | 823 | 0% | 100% | 0% | 0% | 0.1% | 100% |

TNA is not one position. HCC70 is an Arc 1 specialist and mostly epithelial by gene score. MDAMB436 is nearest Arc 2, 100% mesenchymal by gene score, and entirely interior. HCC1143 is nearest Arc 1, 99.8% mesenchymal, and entirely interior. The subtype label does not place the line.

BT549 and HS578T are the clean mesenchymal lines: 100% nearest Arc 3, 100% mesenchymal by gene score, and 0% specialists. They point at the mesenchymal corner without reaching the cut. CAL51 is the opposite mismatch: 100% nearest the epithelial corner, 100% mesenchymal by gene score, 0% specialists.

All 65 hybrid cells are CAL851. That line is 65.9% nearest Arc 3, 33.9% nearest Arc 1, and 0% specialists.

## 2. Corner or middle

| Gene score | Cells | At a corner | In the middle |
|---|---:|---:|---:|
| Epithelial | 11,757 | 80.0% | 20.0% |
| Mesenchymal | 11,757 | 0% | 100% |
| Hybrid | 65 | 0% | 100% |
| Intermediate | 11,692 | 35.3% | 64.7% |

| Position | Cells | Epithelial | Mesenchymal | Hybrid | Intermediate |
|---|---:|---:|---:|---:|---:|
| At a corner | 13,529 | 69.5% | 0% | 0% | 30.5% |
| In the middle | 21,742 | 10.8% | 54.1% | 0.3% | 34.8% |

The 65 hybrid cells have median weights 0.36, 0.28, 0.37. Median purity is 0.01. That is the centroid: equal weight on all three corners. All of them are interior. A cell at a corner is never hybrid and never mesenchymal. Corners in this dataset are epithelial (69.5% of specialists) or intermediate.

## Reading

Joint-high hybrid cells are almost absent, because the two gene scores move against each other. The 65 that exist are one TNB line sitting in the middle of the triangle. Luminal lines are epithelial specialists at Arc 2. Mesenchymal cells, including the TNB lines that point at Arc 3, are interior. The hybrid corner is not a hybrid gene-score compartment here.
