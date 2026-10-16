# Simplex description: gene set all, tag main, k=3

Cells: 23959; inside simplex 73.1%; specialists 40.0%; interior 60.0%
Fit cells inside simplex: 71.8%

## Specialists: vertex x EMT state (counts)

cluster_annotation  Epithelial  Early-Hybrid  Hybrid  Late-EMT
nearest                                                       
V1                         658          1151       1        24
V2                        1030          1869     132         0
V3                           0            33     297      4383

## Specialists: vertex x EMT state (% of row)

cluster_annotation  Epithelial  Early-Hybrid  Hybrid  Late-EMT
nearest                                                       
V1                        35.9          62.8     0.1       1.3
V2                        34.0          61.7     4.4       0.0
V3                         0.0           0.7     6.3      93.0

## Specialists: vertex x clusterU (% of row)

clusterU  Epithelial  Epithelial (Differentiated)  Proliferative Early Hybrid EMT  Early Hybrid EMT_1  Early Hybrid EMT_2  Early Hybrid EMT_3  Proliferative Hybrid EMT  Hybrid EMT  Late EMT
nearest                                                                                                                                                                                      
V1              34.7                          0.9                             0.6                25.4                28.7                 8.4                       0.2         0.9       0.2
V2              25.0                          3.8                            42.0                20.4                 0.7                 3.8                       4.2         0.2       0.0
V3               0.0                          0.0                             0.0                 0.0                 0.5                 0.0                       1.9        25.8      71.8

## Specialists: vertex x phase (% of row)

phase      G1   G2M     S
nearest                  
V1       79.8   8.0  12.3
V2       42.1  32.8  25.1
V3       89.0   3.9   7.0

## Specialists: vertex x assay (% of row)

assay    multiome  scRNA
nearest                 
V1           79.3   20.7
V2           10.4   89.6
V3            0.4   99.6

## Specialists: vertex x tumour (% of row)

sample   MU3  MU4   MU7   MU8  SC15  SC2b1  SC2b3  SC32  SC33  SC35  SC43  SC46
nearest                                                                        
V1       2.3  2.5  44.7  29.8  19.6    0.0    0.0   0.1   0.0   0.0   0.0   1.0
V2       3.3  2.7   1.1   3.3   6.8    0.8    5.6  35.9  34.4   1.5   0.5   4.1
V3       0.0  0.3   0.0   0.0   1.4   18.4   35.4   1.3   0.4  28.8   6.7   7.3

## Each EMT state: % of its cells at each vertex or interior

group                 V1    V2    V3  interior
cluster_annotation                            
Epithelial          15.9  24.8   0.0      59.3
Early-Hybrid        11.5  18.7   0.3      69.4
Hybrid               0.0   3.7   8.2      88.1
Late-EMT             0.4   0.0  70.5      29.1

## Each clusterU: % of its cells at each vertex or interior

group                             V1    V2    V3  interior
clusterU                                                  
Epithelial                      18.5  22.0   0.0      59.5
Epithelial (Differentiated)      3.0  21.7   0.0      75.3
Proliferative Early Hybrid EMT   0.3  37.7   0.0      62.0
Early Hybrid EMT_1              13.9  18.4   0.0      67.8
Early Hybrid EMT_2              23.0   0.9   1.0      75.1
Early Hybrid EMT_3              13.4  10.0   0.1      76.4
Proliferative Hybrid EMT         0.1   6.1   4.2      89.6
Hybrid EMT                       0.5   0.2  33.4      65.9
Late EMT                         0.1   0.0  82.8      17.1

## Each tumour: % of its cells at each vertex or interior

group     V1    V2    V3  interior
sample                            
MU3      4.9  11.9   0.1      83.1
MU4      4.9   8.6   1.7      84.8
MU7     29.3   1.1   0.0      69.5
MU8     10.3   1.9   0.0      87.8
SC15    22.8  13.1   4.1      60.0
SC2b1    0.0   2.1  78.9      19.0
SC2b3    0.0   5.8  56.4      37.8
SC32     0.1  57.0   3.1      39.8
SC33     0.0  75.6   1.3      23.1
SC35     0.0   1.8  52.6      45.6
SC43     0.0   1.3  30.3      68.4
SC46     1.2   8.2  22.7      67.9

## Tumour check per vertex

        n_specialists top_tumour  top_share_pct  n_tumours_ge5pct
vertex                                                           
V1               1834        MU7           44.7                 3
V2               3031       SC32           35.9                 4
V3               4713      SC2b3           35.4                 5

## Mean scores per vertex (specialists) and interior

          KS_EMT  LKPR_Epithelial_andrea  LKPR_EarlyHybridEMT_andrea  LKPR_HybridEMT_andrea  LKPR_LateEMT_andrea  purity
group                                                                                                                   
V1        -0.320                   0.696                       0.804                 -0.069               -0.086   0.479
V2        -0.308                   0.962                       1.444                 -0.027               -0.049   0.573
V3         0.513                  -0.419                       0.211                  0.391                0.613   0.576
interior  -0.062                   0.438                       0.867                  0.107                0.097   0.162

## Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)

                   label_type                         subtype  k_in_bin  fold_enrichment  q_value
archetype                                                                                        
V1                     sample                             MU7       944           3.3701      0.0
V1                     sample                            SC15       421           2.6712      0.0
V1                   clusterU              Early Hybrid EMT_2       584           2.5490      0.0
V1                   clusterU                      Epithelial       788           2.2873      0.0
V1                   clusterU              Early Hybrid EMT_1       730           2.1719      0.0
V1         cluster_annotation                      Epithelial       825           1.9898      0.0
V1                      assay                        multiome      1948           1.9676      0.0
V1                     sample                             MU8       845           1.5919      0.0
V1                   clusterU              Early Hybrid EMT_3       177           1.5458      0.0
V1         cluster_annotation                    Early-Hybrid      1510           1.5122      0.0
V2                     sample                            SC33       900           6.5262      0.0
V2                     sample                            SC32       905           4.7380      0.0
V2                   clusterU  Proliferative Early Hybrid EMT      1049           3.1053      0.0
V2                      phase                             G2M       844           1.8197      0.0
V2                   clusterU     Epithelial (Differentiated)        93           1.7382      0.0
V2         cluster_annotation                      Epithelial       716           1.7269      0.0
V2                      assay                           scRNA      2347           1.6693      0.0
V3                     sample                           SC2b1       806           7.3336      0.0
V3                   clusterU                        Late EMT      1880           4.5975      0.0
V3         cluster_annotation                        Late-EMT      2292           3.6847      0.0
V3                     sample                           SC2b3       773           2.6105      0.0
V3                     sample                            SC35       590           2.2921      0.0
V3                      assay                           scRNA      2396           1.7042      0.0
V3                     sample                            SC43       171           1.6426      0.0

## Top 15 genes higher near each vertex (log-expression difference)

       V1_up          V2_up   V3_up
0     Camk1d          Krt14  Col3a1
1      Maml2           Krt5  Col1a1
2      Pcdh7          Krt6a  Col1a2
3   AY036118          Krt17     Dcn
4     Frmd4b         Krtdap     Bgn
5     Adgrl3           Rbm3   Postn
6      Pard3            Sfn     Fn1
7      Magi1         Lgals7   Aebp1
8      Mast4          Fabp5   Sparc
9     Zswim6           Gclm  Crabp1
10      Ext1  2200002D01Rik    Mmp2
11     Atxn1          H2afz  Col5a2
12     Aopep           Ldha  Cxcl12
13    Ptpn13          Fxyd3     Gsn
14      Nrg1           Hras    Fth1
