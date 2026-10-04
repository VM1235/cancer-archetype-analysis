# DepMap cell lines — same triangle, 63 bulk lines

The 63 DepMap invasive breast lines, projected into the frozen Pal k=3 triangle. Each row is one cell line, not one cell. Gene-score percentiles are computed across these 63 lines. Labels are PAM50 from the Panel B call (Basal 27, LumB 17, Her2 14, LumA 4, Normal 1).

Script: `Breast Cancer/codes/26_specialist_vs_interior.py`. Per-line file: `Breast Cancer/results/specialist_interior/scores_DepMap.csv`.

Epithelial and mesenchymal scores are anti-correlated (r = −0.69). No line is in the top quartile of both, so the hybrid count is 0. The closest line still fails the joint-high cut: the higher of its two percentile ranks, taking the worse of the two, peaks at the 75th percentile and does not clear both.

18 of 63 lines are specialists (28.6%). All 21 mesenchymal lines are interior.

## 1. Where each group sits

| PAM50 | Lines | Nearest hybrid (Arc 1) | Nearest epithelial (Arc 2) | Nearest mesenchymal (Arc 3) | Gene score epithelial | Gene score mesenchymal | At a corner |
|---|---:|---:|---:|---:|---:|---:|---:|
| Basal | 27 | 59.3% | 7.4% | 33.3% | 3.7% | 66.7% | 3.7% |
| LumB | 17 | 5.9% | 88.2% | 5.9% | 41.2% | 11.8% | 64.7% |
| Her2 | 14 | 7.1% | 92.9% | 0% | 71.4% | 7.1% | 35.7% |
| LumA | 4 | 25.0% | 75.0% | 0% | 75.0% | 0% | 0% |
| Normal | 1 | 100% | 0% | 0% | 0% | 0% | 100% |

Basal lines leave the epithelial corner: 59.3% nearest Arc 1, 33.3% nearest Arc 3, and 66.7% mesenchymal by gene score. Only 1 of 27 basal lines passes the corner cut. LumB is the epithelial specialist group (88.2% nearest Arc 2, 64.7% at a corner). Her2 is also nearest Arc 2 (92.9%) and is epithelial by gene score (71.4%), but only 5 of 14 Her2 lines clear the corner cut. Her2 does not have its own corner. LumA is nearest Arc 2 and entirely interior (4 lines; the weight gap stays under 0.35).

The single Normal line is a specialist at Arc 1 with an intermediate gene score. One line is not a corner assignment.

Basal lines, all 27:

| Line | Gene score | Position | Dominant arc | Weights (Arc 1, 2, 3) | Purity |
|---|---|---|---:|---|---:|
| ACH-000223 | intermediate | corner | 1 | 0.65, 0.20, 0.14 | 0.45 |
| ACH-000668 | epithelial | interior | 1 | 0.53, 0.19, 0.28 | 0.25 |
| ACH-000643 | intermediate | interior | 1 | 0.43, 0.16, 0.40 | 0.03 |
| ACH-000573 | mesenchymal | interior | 1 | 0.44, 0.28, 0.29 | 0.15 |
| ACH-000711 | mesenchymal | interior | 1 | 0.43, 0.36, 0.21 | 0.07 |
| ACH-000857 | intermediate | interior | 1 | 0.41, 0.21, 0.38 | 0.04 |
| ACH-000930 | mesenchymal | interior | 1 | 0.48, 0.19, 0.33 | 0.15 |
| ACH-000374 | mesenchymal | interior | 1 | 0.46, 0.16, 0.38 | 0.08 |
| ACH-001705 | intermediate | interior | 1 | 0.44, 0.35, 0.22 | 0.09 |
| ACH-000699 | mesenchymal | interior | 1 | 0.45, 0.26, 0.29 | 0.17 |
| ACH-001390 | mesenchymal | interior | 1 | 0.48, 0.20, 0.32 | 0.16 |
| ACH-001389 | mesenchymal | interior | 1 | 0.39, 0.30, 0.31 | 0.08 |
| ACH-000621 | mesenchymal | interior | 1 | 0.36, 0.36, 0.28 | 0.00 |
| ACH-000849 | intermediate | interior | 1 | 0.43, 0.26, 0.31 | 0.12 |
| ACH-001394 | mesenchymal | interior | 1 | 0.48, 0.17, 0.35 | 0.13 |
| ACH-001662 | mesenchymal | interior | 1 | 0.49, 0.15, 0.36 | 0.12 |
| ACH-000148 | mesenchymal | interior | 2 | 0.33, 0.37, 0.30 | 0.04 |
| ACH-001391 | mesenchymal | interior | 2 | 0.30, 0.42, 0.29 | 0.12 |
| ACH-000276 | mesenchymal | interior | 3 | 0.40, 0.03, 0.58 | 0.18 |
| ACH-002499 | mesenchymal | interior | 3 | 0.37, 0.22, 0.40 | 0.03 |
| ACH-000196 | intermediate | interior | 3 | 0.40, 0.11, 0.49 | 0.08 |
| ACH-000691 | intermediate | interior | 3 | 0.38, 0.08, 0.54 | 0.15 |
| ACH-000536 | intermediate | interior | 3 | 0.35, 0.28, 0.37 | 0.02 |
| ACH-000111 | mesenchymal | interior | 3 | 0.35, 0.24, 0.42 | 0.07 |
| ACH-000212 | mesenchymal | interior | 3 | 0.21, 0.32, 0.48 | 0.16 |
| ACH-000721 | mesenchymal | interior | 3 | 0.16, 0.35, 0.49 | 0.14 |
| ACH-000288 | mesenchymal | interior | 3 | 0.26, 0.27, 0.47 | 0.19 |

ACH-000223 is the only basal specialist, at Arc 1, and its gene score is intermediate rather than hybrid. The other basal lines split between a near-tie of Arc 1 and Arc 3 (purity often under 0.20) and a milder Arc 3 lead that still misses 0.35. ACH-000276 is the strongest Arc 3 basal line (weights 0.40, 0.03, 0.58) and is still interior.

Specialists that do exist are LumB and Her2 at Arc 2. Eleven LumB lines and five Her2 lines pass the cut, all at Arc 2. Examples at the extreme: ACH-002921 (LumB) has weights 0.03, 0.97, 0.00. ACH-000828 (Her2) has weights 0.24, 0.76, 0.00.

## 2. Corner or middle

There is no hybrid row.

| Gene score | Lines | At a corner | In the middle |
|---|---:|---:|---:|
| Epithelial | 21 | 38.1% | 61.9% |
| Mesenchymal | 21 | 0% | 100% |
| Intermediate | 21 | 47.6% | 52.4% |

| Position | Lines | Epithelial | Mesenchymal | Hybrid | Intermediate |
|---|---:|---:|---:|---:|---:|
| At a corner | 18 | 44.4% | 0% | 0% | 55.6% |
| In the middle | 45 | 28.9% | 46.7% | 0% | 24.4% |

A line in the middle is mesenchymal 46.7% of the time and epithelial 28.9% of the time. A line at a corner is never mesenchymal. Eight of the 21 epithelial lines are specialists, all at Arc 2.

## Reading

DepMap has no hybrid line under the top-quartile rule, so there is no hybrid corner and no hybrid generalist. Basal lines are the ones pulled off the epithelial corner, toward Arc 1 and Arc 3, and they stay interior. The lines that lock onto a corner are LumB and Her2 at the epithelial pole. That matches the single-cell cohorts in one respect: mesenchymal state and the basal/TNBC label sit inside the triangle, and the epithelial corner is luminal.
