# Groves closest-5% specialists: all/main, k=4

Reference cells (5% cut-off set on these): 19930
Cells within the 5% cut-off of two vertices: 0
Specialists among reference cells: 20.0%

## state: counts

groves5              V1   V2   V3   V4  generalist  total
cluster_annotation                                       
Epithelial          610  635   38  120        3580   4983
Early-Hybrid        160  136  128  303        4255   4982
Late-Hybrid         111  100  292  350        4129   4982
Mesenchymal         116  126  539  224        3978   4983

## state: % of each state in each class (rows sum to 100)

Read across: where do the cells of this label go?

groves5               V1    V2    V3   V4  generalist
cluster_annotation                                   
Epithelial          12.2  12.7   0.8  2.4        71.8
Early-Hybrid         3.2   2.7   2.6  6.1        85.4
Late-Hybrid          2.2   2.0   5.9  7.0        82.9
Mesenchymal          2.3   2.5  10.8  4.5        79.8

## Each class: % made up of each state (columns sum to 100)

Read down: what is this class made of? Compare with the last column (overall share).

                      V1    V2    V3    V4  generalist  all reference cells
cluster_annotation                                                         
Epithelial          61.2  63.7   3.8  12.0        22.5                 25.0
Early-Hybrid        16.0  13.6  12.8  30.4        26.7                 25.0
Late-Hybrid         11.1  10.0  29.3  35.1        25.9                 25.0
Mesenchymal         11.6  12.6  54.1  22.5        25.0                 25.0

## clusterU: counts

groves5    V1   V2   V3   V4  generalist  total
clusterU                                       
K0         49   26    4  656        5178   5913
K1        392  275    0    3         782   1452
K2        406  556  353   58        2270   3643
K3         49   85  443    3        3088   3668
K4        101   55  197  277        4624   5254

## clusterU: % of each clusterU in each class (rows sum to 100)

Read across: where do the cells of this label go?

groves5     V1    V2    V3    V4  generalist
clusterU                                    
K0         0.8   0.4   0.1  11.1        87.6
K1        27.0  18.9   0.0   0.2        53.9
K2        11.1  15.3   9.7   1.6        62.3
K3         1.3   2.3  12.1   0.1        84.2
K4         1.9   1.0   3.7   5.3        88.0

## Each class: % made up of each clusterU (columns sum to 100)

Read down: what is this class made of? Compare with the last column (overall share).

            V1    V2    V3    V4  generalist  all reference cells
clusterU                                                         
K0         4.9   2.6   0.4  65.8        32.5                 29.7
K1        39.3  27.6   0.0   0.3         4.9                  7.3
K2        40.7  55.8  35.4   5.8        14.2                 18.3
K3         4.9   8.5  44.4   0.3        19.4                 18.4
K4        10.1   5.5  19.8  27.8        29.0                 26.4

## Tumour of each class (% of class)

sample       4T1     E   EM1   EM2   EM3     M
groves5                                       
V1          19.8  37.0  17.4  10.9  13.9   1.0
V2          32.3  22.4  18.1  10.3   6.0  10.9
V3           0.0   0.0   1.2   0.0   0.0  98.8
V4          42.7   0.0   2.5  39.3   6.3   9.1
generalist  31.8   1.3  15.8  17.4  10.8  22.9

## Cell-cycle phase of each class (% of class)

phase         G1   G2M     S
groves5                     
V1          82.9   7.4   9.6
V2          68.9  10.7  20.4
V3          98.0   0.7   1.3
V4          10.7  38.6  50.7
generalist  72.5  11.7  15.7

## Rule comparison, reference cells (counts)

groves5 = closest 5% (Fig 3E); weight095 = weight > 0.95 (STAR Methods); gap035 = purity gap >= 0.35 (Pal work)

            groves5  weight095  gap035
V1              997          1     367
V2              997          0     261
V3              997         17    3030
V4              997         29    4865
generalist    15942      19883   11407

## % of each state called specialist, by rule

                    groves5  weight095  gap035
cluster_annotation                            
Epithelial             28.2        0.1    35.2
Early-Hybrid           14.6        0.2    43.7
Late-Hybrid            17.1        0.3    45.6
Mesenchymal            20.2        0.4    46.6
