# Four datasets on the Pal triangle — the three questions

Pal (GSE161529) is the host. The k=3 archetypes were fit there and the corners were named there: Arc 2 epithelial, Arc 3 mesenchymal, Arc 1 hybrid, from the Tan KS genes on cells nearest each vertex. Wu (GSE176078), GSE173634, and the 63 DepMap lines are projected into that frozen triangle. Nothing was refit.

Shared rules, applied inside each dataset:

- **Gene score.** Hybrid: both epithelial and mesenchymal KS scores at or above that dataset’s 75th percentile. Otherwise epithelial or mesenchymal by the mesenchymal-minus-epithelial score, split at the 33rd and 67th percentiles. The rest are intermediate.
- **Specialist (corner).** Largest archetype weight beats the second by at least 0.35, after clipping negative weights at 0 and renormalizing.
- **Interior, called generalist below.** Every other sample. This includes edges, where two weights are close, and the centroid, where all three weights are near one-third.

Per-sample files are in `Breast Cancer/results/specialist_interior/`. Dataset notes: `RESULT_pal_triangle.md`, `RESULT_wu_gse176078.md`, `RESULT_gse173634.md`, `RESULT_depmap_cell_lines.md`.

| Dataset | Samples | Hybrid | Epi–mes correlation | Specialists |
|---|---:|---:|---:|---:|
| Pal GSE161529 | 84,602 cells | 8,719 (10.3%) | +0.29 | 73.5% |
| Wu GSE176078 | 24,162 cells | 2,413 (10.0%) | −0.29 | 36.0% |
| GSE173634 | 35,271 cells | 65 (0.2%) | −0.77 | 38.4% |
| DepMap | 63 lines | 0 | −0.69 | 28.6% |

## 1. If a sample is a generalist, what is its gene score?

This is P(gene score | interior).

| Dataset | Generalists | Epithelial | Mesenchymal | Hybrid | Intermediate |
|---|---:|---:|---:|---:|---:|
| Pal | 22,402 | 35.6% | 31.1% | 13.4% | 19.9% |
| Wu | 15,473 | 20.9% | 36.7% | 14.6% | 27.9% |
| GSE173634 | 21,742 | 10.8% | 54.1% | 0.3% | 34.8% |
| DepMap | 45 | 28.9% | 46.7% | 0% | 24.4% |

A generalist is not a hybrid. In every dataset the largest shares are epithelial, mesenchymal, or intermediate. Hybrid is 13–15% of generalists in the two tumor single-cell sets, and essentially absent in the two cell-line sets. In GSE173634 and DepMap, the generalist class is where the mesenchymal samples are: 54.1% and 46.7%.

For comparison, the same split among specialists:

| Dataset | Specialists | Epithelial | Mesenchymal | Hybrid | Intermediate |
|---|---:|---:|---:|---:|---:|
| Pal | 62,200 | 28.9% | 30.4% | 9.2% | 31.5% |
| Wu | 8,689 | 45.9% | 16.5% | 1.8% | 35.9% |
| GSE173634 | 13,529 | 69.5% | 0% | 0% | 30.5% |
| DepMap | 18 | 44.4% | 0% | 0% | 55.6% |

Moving from specialist to generalist raises the hybrid share only modestly in Pal (9.2% to 13.4%) and more clearly in Wu (1.8% to 14.6%). In the cell-line datasets the specialist class contains no mesenchymal sample and no hybrid sample.

## 2. If a sample is epithelial, mesenchymal, or hybrid, is it a generalist?

This is P(interior | gene score). The complement is the specialist rate.

| Dataset | Epithelial | Mesenchymal | Hybrid |
|---|---:|---:|---:|
| Pal | 30.7% generalist (69.3% corner) | 27.0% (73.0% corner) | 34.4% (65.6% corner) |
| Wu | 44.7% (55.3% corner) | 79.9% (20.1% corner) | 93.7% (6.3% corner) |
| GSE173634 | 20.0% (80.0% corner) | 100% (0% corner) | 100% (0 of 65 at a corner) |
| DepMap | 61.9% (38.1% corner) | 100% (0 of 21 at a corner) | no hybrid line |

On Pal, hybrid cells are generalists about as often as epithelial and mesenchymal cells are. The hybrid label does not pick out the middle.

On Wu, it does. 93.7% of hybrid cells are generalists, against 44.7% of epithelial cells. Mesenchymal cells are also mostly generalists (79.9%).

On GSE173634 and DepMap, every mesenchymal sample is a generalist, and every hybrid sample (65 cells, all the line CAL851) is a generalist. Epithelial samples are the ones that can be specialists, especially in GSE173634 (80.0%).

## 3. Do hybrid cells sit at a corner, or are they interior generalists?

Each sample was scored from its three weights. Hybrid cells were then counted by that score.

| Dataset | Hybrid samples | At a corner | Of which Arc 1 | Of which Arc 2 | Of which Arc 3 | In the middle | Median weights (Arc 1, 2, 3) |
|---|---:|---:|---:|---:|---:|---:|---|
| Pal | 8,719 | 65.6% | 53.2% | 10.3% | 2.1% | 34.4% | 0.68, 0.02, 0.11 |
| Wu | 2,413 | 6.3% | 0.5% | 5.8% | 0% | 93.7% | 0.51, 0.11, 0.38 |
| GSE173634 | 65 | 0% | 0% | 0% | 0% | 100% | 0.36, 0.28, 0.37 |
| DepMap | 0 | — | — | — | — | — | no line is jointly high |

Pal hybrid cells sit at a corner. The corner is Arc 1. Median purity is 0.52, so this is not a near miss on the cut. Of the hybrid cells that are Arc 1 specialists, 92.2% are one patient, TN-0135, whose median weights are 0.91, 0.00, 0.08.

Wu hybrid cells fail the corner cut. They are not at the centroid either. Arc 1 is the largest weight for 83.3% of them, and Arc 3 is close behind (median 0.38 against 0.51). Median purity is 0.12. That is the edge between the hybrid pole and the mesenchymal pole. Half of Wu’s hybrid cells are one TNBC patient, CID4515, with median weights 0.52, 0.04, 0.43 and a specialist rate of 0%.

GSE173634 hybrid cells are the centroid: weights about one-third each, median purity 0.01, all interior, all from CAL851.

DepMap has no hybrid sample to place. The related observation is about basal and mesenchymal lines: all 21 mesenchymal lines are interior, and 26 of 27 basal lines are interior, spread between Arc 1 and Arc 3. The lines that sit at a corner are LumB and Her2 at Arc 2.

## What holds across the four

Hybrid-as-its-own-corner is a Pal result, and inside Pal it is TN-0135. It does not recur. Where both gene sets are high in another dataset, those cells are interior: an Arc 1–Arc 3 edge in Wu, the centroid in GSE173634, and absent in DepMap.

Generalist does not mean hybrid in any dataset. The middle is mostly epithelial or mesenchymal in Pal, mostly mesenchymal in the cell-line datasets, and mixed in Wu, with hybrid the smallest named share except for intermediate.

The consistent pattern outside Pal is which program reaches a corner. Epithelial samples can be specialists at Arc 2. Mesenchymal samples are interior in Wu (79.9%), GSE173634 (100%), and DepMap (100%). Basal and TNBC samples are the ones that leave the epithelial corner, and in those three datasets they usually do not lock onto a vertex.
