# Cleaning audit: old vs new

## Cook: cleaned, before/after dropping human technical genes (MALAT1, NEAT1, haemoglobin)
Before = old cleaned fit (filter matched nothing).

| Fit | Before | After | Passes (p<0.05) before -> after |
|---|---|---|---|
| all A549 k3 | 0.010 (2646 cells) | 0.010 (2646 cells) | yes -> yes |
| all A549 k4 | 1.000 (2646 cells) | 1.000 (2646 cells) | no -> no |
| all DU145 k3 | 0.249 (1839 cells) | 0.254 (1839 cells) | no -> no |
| all DU145 k4 | 0.338 (1839 cells) | 0.338 (1839 cells) | no -> no |
| all MCF7 k3 | 0.443 (1378 cells) | 0.463 (1378 cells) | no -> no |
| all MCF7 k4 | 0.930 (1378 cells) | 0.930 (1378 cells) | no -> no |
| all OVCA420 k3 | 0.726 (355 cells) | 0.756 (355 cells) | no -> no |
| all OVCA420 k4 | 1.000 (355 cells) | 1.000 (355 cells) | no -> no |
| ks A549 k3 | 0.741 (2646 cells) | 0.741 (2646 cells) | no -> no |
| ks A549 k4 | 0.716 (2646 cells) | 0.716 (2646 cells) | no -> no |
| ks DU145 k3 | 0.950 (1839 cells) | 0.950 (1839 cells) | no -> no |
| ks DU145 k4 | 1.000 (1839 cells) | 1.000 (1839 cells) | no -> no |
| ks MCF7 k3 | 0.030 (1378 cells) | 0.030 (1378 cells) | yes -> yes |
| ks MCF7 k4 | 0.194 (1378 cells) | 0.194 (1378 cells) | no -> no |
| ks OVCA420 k3 | 0.537 (355 cells) | 0.537 (355 cells) | no -> no |
| ks OVCA420 k4 | 0.448 (355 cells) | 0.448 (355 cells) | no -> no |
| all A549 time-balanced k3 | 0.224 (1303 cells) | 0.269 (1303 cells) | no -> no |
| ks A549 time-balanced k3 | 0.726 (1303 cells) | 0.726 (1303 cells) | no -> no |

## GSE114687: cleaned, before/after dropping human technical genes
Before = old cleaned fit (filter matched nothing).

| Fit | Before | After | Passes (p<0.05) before -> after |
|---|---|---|---|
| all MCF10A_Mock k3 | 0.005 (2536 cells) | 0.005 (2536 cells) | yes -> yes |
| all HuMEC_Mock k3 | 0.005 (3889 cells) | 0.005 (3889 cells) | yes -> yes |
| all MCF10A_TGFB k3 | 0.005 (3127 cells) | 0.005 (3127 cells) | yes -> yes |
| all HuMEC_TGFB k3 | 0.025 (4000 cells) | 0.025 (4000 cells) | yes -> yes |
| all MCF10A k3 | 0.005 (4000 cells) | 0.005 (4000 cells) | yes -> yes |
| all HuMEC k3 | 0.925 (4000 cells) | 0.886 (4000 cells) | no -> no |
| ks MCF10A_Mock k3 | 0.005 (2536 cells) | 0.005 (2536 cells) | yes -> yes |
| ks HuMEC_Mock k3 | 0.005 (3889 cells) | 0.005 (3889 cells) | yes -> yes |
| ks MCF10A_TGFB k3 | 0.562 (3127 cells) | 0.562 (3127 cells) | no -> no |
| ks HuMEC_TGFB k3 | 0.005 (4000 cells) | 0.005 (4000 cells) | yes -> yes |
| ks MCF10A k3 | 0.726 (4000 cells) | 0.726 (4000 cells) | no -> no |
| ks HuMEC k3 | 0.075 (4000 cells) | 0.075 (4000 cells) | no -> no |
| all MCF10A_Mock k4 | 0.005 (2536 cells) | 0.005 (2536 cells) | yes -> yes |
| all HuMEC_Mock k4 | 0.005 (3889 cells) | 0.005 (3889 cells) | yes -> yes |
| all MCF10A_TGFB k4 | 0.005 (3127 cells) | 0.005 (3127 cells) | yes -> yes |
| all HuMEC_TGFB k4 | 0.005 (4000 cells) | 0.005 (4000 cells) | yes -> yes |
| all MCF10A k4 | 0.005 (4000 cells) | 0.005 (4000 cells) | yes -> yes |
| all HuMEC k4 | 0.005 (4000 cells) | 0.005 (4000 cells) | yes -> yes |
| ks MCF10A_Mock k4 | 0.005 (2536 cells) | 0.005 (2536 cells) | yes -> yes |
| ks HuMEC_Mock k4 | 0.005 (3889 cells) | 0.005 (3889 cells) | yes -> yes |
| ks MCF10A_TGFB k4 | 0.393 (3127 cells) | 0.393 (3127 cells) | no -> no |
| ks HuMEC_TGFB k4 | 0.005 (4000 cells) | 0.005 (4000 cells) | yes -> yes |
| ks MCF10A k4 | 0.005 (4000 cells) | 0.005 (4000 cells) | yes -> yes |
| ks HuMEC k4 | 0.169 (4000 cells) | 0.169 (4000 cells) | no -> no |

## GSE114687: UMI cutoff 10,000 -> 5,000
Before = 10,000 with human technical genes dropped.

| Fit | Before | After | Passes (p<0.05) before -> after |
|---|---|---|---|
| all MCF10A_Mock k3 | 0.005 (2536 cells) | 0.005 (4000 cells) | yes -> yes |
| all HuMEC_Mock k3 | 0.005 (3889 cells) | 0.005 (4000 cells) | yes -> yes |
| all MCF10A_TGFB k3 | 0.005 (3127 cells) | 0.005 (4000 cells) | yes -> yes |
| all HuMEC_TGFB k3 | 0.025 (4000 cells) | 1.000 (4000 cells) | yes -> no |
| all MCF10A k3 | 0.005 (4000 cells) | 0.005 (4000 cells) | yes -> yes |
| all HuMEC k3 | 0.886 (4000 cells) | 0.005 (4000 cells) | no -> yes |
| ks MCF10A_Mock k3 | 0.005 (2536 cells) | 1.000 (4000 cells) | yes -> no |
| ks HuMEC_Mock k3 | 0.005 (3889 cells) | 0.085 (4000 cells) | yes -> no |
| ks MCF10A_TGFB k3 | 0.562 (3127 cells) | 0.751 (4000 cells) | no -> no |
| ks HuMEC_TGFB k3 | 0.005 (4000 cells) | 0.080 (4000 cells) | yes -> no |
| ks MCF10A k3 | 0.726 (4000 cells) | 0.557 (4000 cells) | no -> no |
| ks HuMEC k3 | 0.075 (4000 cells) | 0.005 (4000 cells) | no -> yes |

## GSE114687: UMI cutoff 10,000 -> 15,000
Before = 10,000 with human technical genes dropped.

| Fit | Before | After | Passes (p<0.05) before -> after |
|---|---|---|---|
| all MCF10A_Mock k3 | 0.005 (2536 cells) | 0.005 (542 cells) | yes -> yes |
| all HuMEC_Mock k3 | 0.005 (3889 cells) | 0.015 (573 cells) | yes -> yes |
| all MCF10A_TGFB k3 | 0.005 (3127 cells) | 0.005 (1296 cells) | yes -> yes |
| all HuMEC_TGFB k3 | 0.025 (4000 cells) | 0.005 (284 cells) | yes -> yes |
| all MCF10A k3 | 0.005 (4000 cells) | 0.005 (1838 cells) | yes -> yes |
| all HuMEC k3 | 0.886 (4000 cells) | 0.005 (857 cells) | no -> yes |
| ks MCF10A_Mock k3 | 0.005 (2536 cells) | 0.005 (542 cells) | yes -> yes |
| ks HuMEC_Mock k3 | 0.005 (3889 cells) | 0.005 (573 cells) | yes -> yes |
| ks MCF10A_TGFB k3 | 0.562 (3127 cells) | 0.353 (1296 cells) | no -> no |
| ks HuMEC_TGFB k3 | 0.005 (4000 cells) | 0.005 (284 cells) | yes -> yes |
| ks MCF10A k3 | 0.726 (4000 cells) | 0.622 (1838 cells) | no -> no |
| ks HuMEC k3 | 0.075 (4000 cells) | 0.005 (857 cells) | no -> yes |

## 4T1: contamination cutoff 0.3 -> 0.2
Within-clone, cleaned (>=10,000 UMI).

| Fit | Before | After | Passes (p<0.05) before -> after |
|---|---|---|---|
| all EM1 k3 | 0.622 (1365 cells) | 0.657 (1332 cells) | no -> no |
| all EM2 k3 | 0.627 (2701 cells) | 0.682 (2562 cells) | no -> no |
| all EM3 k3 | 0.512 (1507 cells) | 0.259 (1208 cells) | no -> no |
| all M k3 | 0.498 (3324 cells) | 0.602 (2763 cells) | no -> no |
| all 4T1 k3 | 0.284 (4000 cells) | 0.199 (3821 cells) | no -> no |
| ks EM1 k3 | 0.826 (1365 cells) | 0.617 (1332 cells) | no -> no |
| ks EM2 k3 | 0.015 (2701 cells) | 0.035 (2562 cells) | yes -> yes |
| ks EM3 k3 | 0.856 (1507 cells) | 0.945 (1208 cells) | no -> no |
| ks M k3 | 0.164 (3324 cells) | 0.294 (2763 cells) | no -> no |
| ks 4T1 k3 | 0.100 (4000 cells) | 0.100 (3821 cells) | no -> no |

## 4T1: contamination cutoff 0.3 -> 0.5
Within-clone, cleaned (>=10,000 UMI).

| Fit | Before | After | Passes (p<0.05) before -> after |
|---|---|---|---|
| all EM1 k3 | 0.622 (1365 cells) | 0.701 (1373 cells) | no -> no |
| all EM2 k3 | 0.627 (2701 cells) | 0.572 (2738 cells) | no -> no |
| all EM3 k3 | 0.512 (1507 cells) | 0.612 (1582 cells) | no -> no |
| all M k3 | 0.498 (3324 cells) | 0.373 (3525 cells) | no -> no |
| all 4T1 k3 | 0.284 (4000 cells) | 0.433 (4000 cells) | no -> no |
| ks EM1 k3 | 0.826 (1365 cells) | 0.846 (1373 cells) | no -> no |
| ks EM2 k3 | 0.015 (2701 cells) | 0.050 (2738 cells) | yes -> yes |
| ks EM3 k3 | 0.856 (1507 cells) | 0.741 (1582 cells) | no -> no |
| ks M k3 | 0.164 (3324 cells) | 0.179 (3525 cells) | no -> no |
| ks 4T1 k3 | 0.100 (4000 cells) | 0.005 (4000 cells) | no -> yes |

## PDAC: cleaned KS-gene fits (never run before)
No earlier result to compare.

| Fit | Before | After | Passes (p<0.05) before -> after |
|---|---|---|---|
| ks M2-Met k3 | n/a | 0.124 (4000 cells) | - -> no |
| ks M2-PT k3 | n/a | 0.622 (4000 cells) | - -> no |
| ks M2-Liver k3 | n/a | 0.856 (4000 cells) | - -> no |
| ks M1-PT k3 | n/a | 0.990 (3770 cells) | - -> no |
| ks M1-Met k3 | n/a | 0.627 (4000 cells) | - -> no |
| ks M1-Liver k3 | n/a | 0.995 (3521 cells) | - -> no |
| ks M1-SS k3 | n/a | 0.333 (2397 cells) | - -> no |
| ks M1 k3 | n/a | 1.000 (4000 cells) | - -> no |
| ks M2 k3 | n/a | 0.866 (2728 cells) | - -> no |
| ks M2-Met k4 | n/a | 0.040 (4000 cells) | - -> yes |
| ks M2-PT k4 | n/a | 0.910 (4000 cells) | - -> no |
| ks M2-Liver k4 | n/a | 0.746 (4000 cells) | - -> no |
| ks M1-PT k4 | n/a | 0.154 (3770 cells) | - -> no |
| ks M1-Met k4 | n/a | 0.597 (4000 cells) | - -> no |
| ks M1-Liver k4 | n/a | 0.970 (3521 cells) | - -> no |
| ks M1-SS k4 | n/a | 0.005 (2397 cells) | - -> yes |
| ks M1 k4 | n/a | 0.985 (4000 cells) | - -> no |
| ks M2 k4 | n/a | 0.025 (2728 cells) | - -> yes |
