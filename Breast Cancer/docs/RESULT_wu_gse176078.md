# Wu (GSE176078) — same triangle, same two questions

24,162 Wu cells, projected into the frozen Pal k=3 triangle. No new archetypes were fit. Weights are barycentric coordinates in that triangle, with negative weights clipped at 0 and renormalized. Gene scores use the Tan KS lists, with percentiles computed inside Wu, not borrowed from Pal. Corner rule is the same: largest weight beats the second by at least 0.35, or the cell is interior.

Script: `Breast Cancer/codes/26_specialist_vs_interior.py`. Per-cell file: `Breast Cancer/results/specialist_interior/scores_GSE176078.csv`.

Epithelial and mesenchymal scores in Wu are negatively correlated (r = −0.29). Hybrid cells still exist: 2,413 cells, 10.0% of the dataset. Only 36.0% of all Wu cells pass the corner cut (8,689 / 24,162), against 73.5% in Pal.

## 1. Where each group sits

| Group | Cells | Nearest hybrid (Arc 1) | Nearest epithelial (Arc 2) | Nearest mesenchymal (Arc 3) | Gene score hybrid | Gene score epithelial | Gene score mesenchymal | At a corner |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ER+ | 11,860 | 11.3% | 88.7% | 0% | 3.4% | 40.8% | 15.7% | 66.6% |
| HER2+ | 1,756 | 63.3% | 36.7% | 0% | 23.6% | 46.4% | 17.3% | 29.5% |
| TNBC | 10,546 | 34.6% | 4.4% | 61.0% | 15.2% | 14.8% | 46.8% | 2.5% |

ER+ sits nearest the epithelial corner and is the group that actually reaches it: 66.6% of ER+ cells are specialists. TNBC is mostly nearest the mesenchymal corner (61.0%) or the hybrid corner (34.6%), but 97.5% of TNBC cells fail the corner cut. They have a largest weight, and that weight is not far enough ahead of the next one to count as a corner.

HER2+ is split between Arc 1 (63.3%) and Arc 2 (36.7%), with a high hybrid gene-score rate (23.6%). Only 29.5% of HER2+ cells are specialists, so most of that Arc 1 mass is inside the triangle.

Patients with at least 200 cells:

| Patient | Subtype | Cells | Nearest Arc 1 | Nearest Arc 2 | Nearest Arc 3 | Hybrid | Epithelial | Mesenchymal | At a corner |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CID4290A | ER+ | 4,052 | 0% | 100% | 0% | 0% | 80.1% | 0% | 93.9% |
| CID44991 | TNBC | 4,016 | 12.7% | 5.8% | 81.5% | 2.1% | 38.8% | 13.3% | 0% |
| CID4067 | ER+ | 2,351 | 0% | 100% | 0% | 0% | 0.2% | 23.4% | 31.2% |
| CID4535 | ER+ | 2,213 | 0% | 100% | 0% | 0.1% | 0% | 40.6% | 99.9% |
| CID4515 | TNBC | 2,166 | 99.6% | 0% | 0.4% | 55.7% | 0% | 22.0% | 0% |
| CID4530N | ER+ | 1,715 | 78.0% | 22.0% | 0% | 14.8% | 58.3% | 0.1% | 0% |
| CID4523 | TNBC | 1,158 | 56.2% | 5.4% | 38.4% | 26.5% | 0.3% | 61.7% | 2.8% |
| CID4495 | TNBC | 1,062 | 9.2% | 0.8% | 89.9% | 0% | 0% | 100% | 0% |
| CID4513 | TNBC | 926 | 11.9% | 0% | 88.1% | 0% | 0% | 100% | 0% |
| CID44971 | TNBC | 883 | 0.1% | 0.3% | 99.5% | 0% | 0% | 99.9% | 26.5% |
| CID45171 | HER2+ | 799 | 99.1% | 0.9% | 0% | 0% | 47.6% | 37.4% | 0% |
| CID4463 | ER+ | 659 | 0% | 100% | 0% | 0% | 89.1% | 0.3% | 100% |
| CID4066 | HER2+ | 520 | 0% | 100% | 0% | 26.7% | 62.7% | 0.8% | 99.6% |
| CID3921 | HER2+ | 437 | 73.2% | 26.8% | 0% | 63.2% | 24.9% | 0% | 0% |

The Pal-like patient is CID4515: TNBC, 99.6% nearest Arc 1, 55.7% gene-score hybrid. It is not a Pal-like corner. Median weights for CID4515 are 0.52, 0.03, 0.43. Arc 1 leads Arc 3 by about 0.10, and the corner cut is 0.35, so none of its 2,166 cells are specialists. CID4515 is the Arc 1–Arc 3 edge.

The other large TNBC patients (CID4495, CID4513, CID44971, CID44991) are nearest Arc 3 and almost entirely interior. CID4495 and CID4513 are 100% mesenchymal by gene score and 0% at a corner. Their largest weight is Arc 3, but the other weights are close enough that they do not count as a mesenchymal specialist.

ER+ corners that do lock in are CID4290A, CID4535, and CID4463, all nearest Arc 2 and at least 94% specialists.

## 2. Corner or middle

| Gene score | Cells | At a corner | In the middle |
|---|---:|---:|---:|
| Epithelial | 7,217 | 55.3% | 44.7% |
| Mesenchymal | 7,101 | 20.1% | 79.9% |
| Hybrid | 2,413 | 6.3% | 93.7% |
| Intermediate | 7,431 | 41.9% | 58.1% |

| Position | Cells | Epithelial | Mesenchymal | Hybrid | Intermediate |
|---|---:|---:|---:|---:|---:|
| At a corner | 8,689 | 45.9% | 16.5% | 1.8% | 35.9% |
| In the middle | 15,473 | 20.9% | 36.7% | 14.6% | 27.9% |

Hybrid cells: median weights 0.51, 0.11, 0.38. Median purity is 0.12. Arc 1 is the largest weight for 83.3% of them, but only 0.5% are specialists at Arc 1 and 5.8% are specialists at Arc 2. None are specialists at Arc 3. Half of all Wu hybrid cells (1,207 / 2,413) are CID4515, and none of those are specialists.

A cell in the middle is hybrid 14.6% of the time, against 1.8% for a cell at a corner. Mesenchymal cells are interior 79.9% of the time. Epithelial cells are the ones that reach a corner (55.3%).

## Reading

Wu hybrid cells are interior. They are not in the centroid. The typical hybrid weight vector is about half Arc 1 and 0.38 Arc 3, which is the edge between the hybrid pole and the mesenchymal pole. The one TNBC patient that looks like Pal’s TN-0135 (CID4515, nearest Arc 1, mostly hybrid) sits on that edge, with Arc 1 only slightly ahead of Arc 3. The epithelial corner is an ER+ corner.
