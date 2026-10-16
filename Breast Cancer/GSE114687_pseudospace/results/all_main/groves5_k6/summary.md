# Groves closest-5% specialists: all/main, k=6

Reference cells (5% cut-off set on these): 30590
Cells within the 5% cut-off of two vertices: 0
Specialists among reference cells: 30.0%

## state: counts

groves5               V1   V2   V3   V4   V5    V6  generalist  total
cluster_annotation                                                   
Epithelial             0  607   33  112  281  1514        5101   7648
Early-Hybrid          17  458  239  375  687    16        5855   7647
Late-Hybrid          132  283  479  605  484     0        5664   7647
Mesenchymal         1381  182  779  438   78     0        4790   7648

## state: % of each state in each class (rows sum to 100)

Read across: where do the cells of this label go?

groves5               V1   V2    V3   V4   V5    V6  generalist
cluster_annotation                                             
Epithelial           0.0  7.9   0.4  1.5  3.7  19.8        66.7
Early-Hybrid         0.2  6.0   3.1  4.9  9.0   0.2        76.6
Late-Hybrid          1.7  3.7   6.3  7.9  6.3   0.0        74.1
Mesenchymal         18.1  2.4  10.2  5.7  1.0   0.0        62.6

## Each class: % made up of each state (columns sum to 100)

Read down: what is this class made of? Compare with the last column (overall share).

                      V1    V2    V3    V4    V5    V6  generalist  all reference cells
cluster_annotation                                                                     
Epithelial           0.0  39.7   2.2   7.3  18.4  99.0        23.8                 25.0
Early-Hybrid         1.1  29.9  15.6  24.5  44.9   1.0        27.3                 25.0
Late-Hybrid          8.6  18.5  31.3  39.5  31.6   0.0        26.5                 25.0
Mesenchymal         90.3  11.9  50.9  28.6   5.1   0.0        22.4                 25.0

## clusterU: counts

groves5     V1   V2   V3   V4   V5    V6  generalist  total
clusterU                                                   
K0           5  449    0    5  844     0        8760  10063
K1           0  173    0    3    1  1358        2369   3904
K2        1505    4  981  633    0     0        1877   5000
K3          20   53  549  888    0   172        1815   3497
K4           0  851    0    1  685     0        6589   8126

## clusterU: % of each clusterU in each class (rows sum to 100)

Read across: where do the cells of this label go?

groves5     V1    V2    V3    V4   V5    V6  generalist
clusterU                                               
K0         0.0   4.5   0.0   0.0  8.4   0.0        87.1
K1         0.0   4.4   0.0   0.1  0.0  34.8        60.7
K2        30.1   0.1  19.6  12.7  0.0   0.0        37.5
K3         0.6   1.5  15.7  25.4  0.0   4.9        51.9
K4         0.0  10.5   0.0   0.0  8.4   0.0        81.1

## Each class: % made up of each clusterU (columns sum to 100)

Read down: what is this class made of? Compare with the last column (overall share).

            V1    V2    V3    V4    V5    V6  generalist  all reference cells
clusterU                                                                     
K0         0.3  29.3   0.0   0.3  55.2   0.0        40.9                 32.9
K1         0.0  11.3   0.0   0.2   0.1  88.8        11.1                 12.8
K2        98.4   0.3  64.1  41.4   0.0   0.0         8.8                 16.3
K3         1.3   3.5  35.9  58.0   0.0  11.2         8.5                 11.4
K4         0.0  55.6   0.0   0.1  44.8   0.0        30.8                 26.6

## Tumour of each class (% of class)

sample      HuMEC_Mock  HuMEC_TGFB  MCF10A_Mock  MCF10A_TGFB
groves5                                                     
V1                 0.0         0.0          0.0        100.0
V2                66.7        33.3          0.0          0.0
V3                 0.0         0.0         71.2         28.8
V4                 0.0         0.0         90.5          9.5
V5                48.2        51.8          0.0          0.0
V6                34.9        12.0         52.4          0.7
generalist        41.7        40.5          8.0          9.8

## Cell-cycle phase of each class (% of class)

phase         G1   G2M     S
groves5                     
V1          83.9   5.1  11.0
V2          67.0   6.0  27.0
V3           7.9  47.3  44.8
V4          82.1   1.6  16.3
V5          62.0   1.4  36.7
V6          97.9   0.1   2.0
generalist  60.8  11.6  27.7

## Rule comparison, reference cells (counts)

groves5 = closest 5% (Fig 3E); weight095 = weight > 0.95 (STAR Methods); gap035 = purity gap >= 0.35 (Pal work)

            groves5  weight095  gap035
V1             1530          0    1606
V2             1530          1    1458
V3             1530          3     649
V4             1530          1     698
V5             1530          3    9964
V6             1530          0     604
generalist    21410      30582   15611

## % of each state called specialist, by rule

                    groves5  weight095  gap035
cluster_annotation                            
Epithelial             33.3        0.0    42.0
Early-Hybrid           23.4        0.0    54.4
Late-Hybrid            25.9        0.0    53.5
Mesenchymal            37.4        0.0    45.9
