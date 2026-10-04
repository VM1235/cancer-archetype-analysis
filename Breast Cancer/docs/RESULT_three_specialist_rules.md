# Three specialist rules on the Pal triangle

Same cells, same clipped barycentric weights, same KS program labels as `26_specialist_vs_interior.py`. Distances for the 5% rule are Euclidean distance in Pal PC1–PC2 to the frozen k=3 vertices. Weight on a vertex and distance to that vertex correlate between −0.87 and −0.99, so the two descriptions agree about which corner a cell is near.

| Rule | Source | A sample is a specialist when |
|---|---|---|
| Closest 5% | Groves et al. 2022, Figure 3E | It is among the closest 5% of its own dataset to one vertex. Ties across vertices did not occur. |
| Leading weight > 0.5 | AAnet; Groves 2021 preprint | One archetype weight is above one half. |
| Purity gap ≥ 0.35 | Scripts 24 and 26 | Largest weight beats the second by at least 0.35. Not from those papers. |

The 5% rule labels 15% of each single-cell dataset as specialists (19% of the 63 DepMap lines, because 5% of 63 rounds up to 4 lines per vertex).

## Are hybrid cells specialists?

Percent of hybrid samples called specialist:

| Dataset | Closest 5% | Weight > 0.5 | Gap ≥ 0.35 |
|---|---:|---:|---:|
| Pal | 29.8 | 96.4 | 65.6 |
| Wu | 11.3 | 86.9 | 6.3 |
| GSE173634 | 0 | 0 | 0 |
| DepMap | no hybrid line | no hybrid line | no hybrid line |

Pal hybrids are enriched in the specialist class under the 5% rule (29.8% of hybrids, against 15% of all cells). That enrichment is Arc 1: 61% of the cells closest to Arc 1 are hybrid, against 10.3% of Pal overall (about 6-fold). Every one of those Arc-1 hybrid cells is TN-0135. The tails at Arc 2 and Arc 3 are 0.4% and 0% hybrid.

Wu hybrids are not enriched among specialists under the 5% rule (11.3% versus 15% of all Wu cells). They are enriched inside the Arc 1 tail only: 22.5% of that tail is hybrid, against 10.0% of Wu (about 2-fold). 44% of those Arc-1 hybrid cells are CID4515. The weight > 0.5 rule calls 87% of Wu hybrids specialists, mostly at Arc 1, because their leading weight is about 0.51. The gap rule calls 6% of them specialists. Those two paper-backed rules disagree with each other on Wu, and the gap rule agrees with the 5% rule that Wu hybrids are not, as a class, the extreme cells.

GSE173634 hybrid cells (65 cells, all CAL851) are specialists under none of the three rules.

## Is a generalist a hybrid?

Percent of generalists that are hybrid. Pal’s base rate is 10.3%. Wu’s is 10.0%.

| Dataset | Closest 5% | Weight > 0.5 | Gap ≥ 0.35 |
|---|---:|---:|---:|
| Pal | 8.5 | 7.3 | 13.4 |
| Wu | 10.4 | 4.1 | 14.6 |
| GSE173634 | 0.2 | 0.4 | 0.3 |
| DepMap | 0 | 0 | 0 |

Under all three rules the generalist class is not a hybrid class.

## What is stable

The closest-5% rule and the gap rule agree on the sentence that matters: hybrid cells are not the generalists, and the Pal hybrid extreme is TN-0135 at Arc 1. The majority-weight rule agrees that generalists are not hybrids, but on Pal it calls 95% of all cells specialists, so it barely leaves an interior, and on Wu it turns the Arc 1–Arc 3 edge into a specialist call.

Tables: `results/specialist_interior/rule_comparison/`.
