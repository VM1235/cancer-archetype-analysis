# Groves closest-5% specialists: all/clean, k=3

Reference cells (5% cut-off set on these): 13256
Cells within the 5% cut-off of two vertices: 0
Specialists among reference cells: 15.0%

## state: counts

groves5              V1   V2   V3  generalist  total
cluster_annotation                                  
Epithelial           72  228   49        2123   2472
Early-Hybrid        174  181  108        2977   3440
Late-Hybrid         214  147  186        3121   3668
Mesenchymal         203  107  320        3046   3676

## state: % of each state in each class (rows sum to 100)

Read across: where do the cells of this label go?

groves5              V1   V2   V3  generalist
cluster_annotation                           
Epithelial          2.9  9.2  2.0        85.9
Early-Hybrid        5.1  5.3  3.1        86.5
Late-Hybrid         5.8  4.0  5.1        85.1
Mesenchymal         5.5  2.9  8.7        82.9

## Each class: % made up of each state (columns sum to 100)

Read down: what is this class made of? Compare with the last column (overall share).

                      V1    V2    V3  generalist  all reference cells
cluster_annotation                                                   
Epithelial          10.9  34.4   7.4        18.8                 18.6
Early-Hybrid        26.2  27.3  16.3        26.4                 26.0
Late-Hybrid         32.3  22.2  28.1        27.7                 27.7
Mesenchymal         30.6  16.1  48.3        27.0                 27.7

## clusterU: counts

groves5    V1   V2   V3  generalist  total
clusterU                                  
K0        458  322    1        4076   4857
K1         11   73    2         389    475
K2        128    8  224        1246   1606
K3          1   45  366        1795   2207
K4         65  215   70        3761   4111

## clusterU: % of each clusterU in each class (rows sum to 100)

Read across: where do the cells of this label go?

groves5    V1    V2    V3  generalist
clusterU                             
K0        9.4   6.6   0.0        83.9
K1        2.3  15.4   0.4        81.9
K2        8.0   0.5  13.9        77.6
K3        0.0   2.0  16.6        81.3
K4        1.6   5.2   1.7        91.5

## Each class: % made up of each clusterU (columns sum to 100)

Read down: what is this class made of? Compare with the last column (overall share).

            V1    V2    V3  generalist  all reference cells
clusterU                                                   
K0        69.1  48.6   0.2        36.2                 36.6
K1         1.7  11.0   0.3         3.5                  3.6
K2        19.3   1.2  33.8        11.1                 12.1
K3         0.2   6.8  55.2        15.9                 16.6
K4         9.8  32.4  10.6        33.4                 31.0

## Tumour of each class (% of class)

sample       4T1    E   EM1   EM2   EM3     M
groves5                                      
V1          51.4  0.0   3.2  15.4   3.3  26.7
V2          56.0  4.2   2.9  22.6  14.3   0.0
V3           0.3  0.2   6.6   0.5   2.6  89.9
generalist  31.1  1.0  11.4  21.7  12.2  22.6

## Cell-cycle phase of each class (% of class)

phase         G1   G2M     S
groves5                     
V1           1.5  47.5  51.0
V2          99.4   0.0   0.6
V3          97.0   1.2   1.8
generalist  66.4  14.5  19.1

## Rule comparison, reference cells (counts)

groves5 = closest 5% (Fig 3E); weight095 = weight > 0.95 (STAR Methods); gap035 = purity gap >= 0.35 (Pal work)

            groves5  weight095  gap035
V1              663          5    1295
V2              663          1    2320
V3              663          2     766
generalist    11267      13248    8875

## % of each state called specialist, by rule

                    groves5  weight095  gap035
cluster_annotation                            
Epithelial             14.1        0.0    35.0
Early-Hybrid           13.5        0.0    32.3
Late-Hybrid            14.9        0.0    31.7
Mesenchymal            17.1        0.1    33.8
