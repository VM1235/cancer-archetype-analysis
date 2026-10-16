# Groves closest-5% specialists: all/scrna, k=3

Reference cells (5% cut-off set on these): 14059
Cells within the 5% cut-off of two vertices: 0
Specialists among reference cells: 15.0%

## state: counts

groves5              V1   V2   V3  generalist  total
cluster_annotation                                  
Epithelial          295    0  283         926   1504
Early-Hybrid        311    0  404        2756   3471
Hybrid               14   10   16        3388   3428
Late-EMT             83  693    0        4880   5656

## state: % of each state in each class (rows sum to 100)

Read across: where do the cells of this label go?

groves5               V1    V2    V3  generalist
cluster_annotation                              
Epithelial          19.6   0.0  18.8        61.6
Early-Hybrid         9.0   0.0  11.6        79.4
Hybrid               0.4   0.3   0.5        98.8
Late-EMT             1.5  12.3   0.0        86.3

## Each class: % made up of each state (columns sum to 100)

Read down: what is this class made of? Compare with the last column (overall share).

                      V1    V2    V3  generalist  all reference cells
cluster_annotation                                                   
Epithelial          42.0   0.0  40.3         7.7                 10.7
Early-Hybrid        44.2   0.0  57.5        23.1                 24.7
Hybrid               2.0   1.4   2.3        28.4                 24.4
Late-EMT            11.8  98.6   0.0        40.8                 40.2

## clusterU: counts

groves5                          V1   V2   V3  generalist  total
clusterU                                                        
Epithelial                      288    0  184         772   1244
Epithelial (Differentiated)       9    0   25         116    150
Proliferative Early Hybrid EMT   12    0  407        1006   1425
Early Hybrid EMT_1               89    0   63         786    938
Early Hybrid EMT_2              155    0    0         591    746
Early Hybrid EMT_3               51    0   10         396    457
Proliferative Hybrid EMT         15    2   14        1954   1985
Hybrid EMT                       49  126    0        3120   3295
Late EMT                         35  575    0        3209   3819

## clusterU: % of each clusterU in each class (rows sum to 100)

Read across: where do the cells of this label go?

groves5                           V1    V2    V3  generalist
clusterU                                                    
Epithelial                      23.2   0.0  14.8        62.1
Epithelial (Differentiated)      6.0   0.0  16.7        77.3
Proliferative Early Hybrid EMT   0.8   0.0  28.6        70.6
Early Hybrid EMT_1               9.5   0.0   6.7        83.8
Early Hybrid EMT_2              20.8   0.0   0.0        79.2
Early Hybrid EMT_3              11.2   0.0   2.2        86.7
Proliferative Hybrid EMT         0.8   0.1   0.7        98.4
Hybrid EMT                       1.5   3.8   0.0        94.7
Late EMT                         0.9  15.1   0.0        84.0

## Each class: % made up of each clusterU (columns sum to 100)

Read down: what is this class made of? Compare with the last column (overall share).

                                  V1    V2    V3  generalist  all reference cells
clusterU                                                                         
Epithelial                      41.0   0.0  26.2         6.5                  8.8
Epithelial (Differentiated)      1.3   0.0   3.6         1.0                  1.1
Proliferative Early Hybrid EMT   1.7   0.0  57.9         8.4                 10.1
Early Hybrid EMT_1              12.7   0.0   9.0         6.6                  6.7
Early Hybrid EMT_2              22.0   0.0   0.0         4.9                  5.3
Early Hybrid EMT_3               7.3   0.0   1.4         3.3                  3.3
Proliferative Hybrid EMT         2.1   0.3   2.0        16.4                 14.1
Hybrid EMT                       7.0  17.9   0.0        26.1                 23.4
Late EMT                         5.0  81.8   0.0        26.9                 27.2

## Tumour of each class (% of class)

sample      SC15  SC2b1  SC2b3  SC32  SC33  SC35  SC43  SC46
groves5                                                     
V1          82.4    0.0    1.1   2.7   1.1   0.9   0.6  11.2
V2           0.0   92.2    6.1   0.0   0.0   1.0   0.7   0.0
V3           2.0    1.6    4.4  41.5  48.1   0.9   0.4   1.1
generalist   8.2    3.7   24.1  13.4   8.6  21.4   8.6  12.0

## Cell-cycle phase of each class (% of class)

phase         G1   G2M     S
groves5                     
V1          66.3  16.5  17.2
V2          89.2   2.6   8.3
V3          20.1  40.1  39.8
generalist  67.2  18.0  14.9

## Rule comparison, reference cells (counts)

groves5 = closest 5% (Fig 3E); weight095 = weight > 0.95 (STAR Methods); gap035 = purity gap >= 0.35 (Pal work)

            groves5  weight095  gap035
V1              703          4     387
V2              703         80    3797
V3              703         35    3212
generalist    11950      13940    6663

## % of each state called specialist, by rule

                    groves5  weight095  gap035
cluster_annotation                            
Epithelial             38.4        1.7    71.7
Early-Hybrid           20.6        0.4    61.1
Hybrid                  1.2        0.0    17.6
Late-EMT               13.7        1.4    63.5

## Non-reference cells (e.g. multiome), same cut-offs: % of each state

groves5               V1   V2   V3  generalist
cluster_annotation                            
Epithelial          49.3  0.0  0.1        50.6
Early-Hybrid        54.1  0.0  0.0        45.9
Hybrid               8.3  0.0  0.0        91.7
Late-EMT            40.6  0.0  0.0        59.4
