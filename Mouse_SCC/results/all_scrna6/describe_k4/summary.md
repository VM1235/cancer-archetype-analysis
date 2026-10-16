# Simplex description: gene set all, tag scrna6, k=4

Cells: 23959; inside simplex 20.6%; specialists 33.5%; interior 66.5%
Fit cells inside simplex: 23.2%

## Specialists: vertex x EMT state (counts)

cluster_annotation  Epithelial  Early-Hybrid  Hybrid  Late-EMT
nearest                                                       
V1                           0             1    1607       136
V2                           0             1     175       726
V3                        1368          2364       1         3
V4                          87           495      11      1059

## Specialists: vertex x EMT state (% of row)

cluster_annotation  Epithelial  Early-Hybrid  Hybrid  Late-EMT
nearest                                                       
V1                         0.0           0.1    92.1       7.8
V2                         0.0           0.1    19.4      80.5
V3                        36.6          63.3     0.0       0.1
V4                         5.3          30.0     0.7      64.1

## Specialists: vertex x clusterU (% of row)

clusterU  Epithelial  Epithelial (Differentiated)  Proliferative Early Hybrid EMT  Early Hybrid EMT_1  Early Hybrid EMT_2  Early Hybrid EMT_3  Proliferative Hybrid EMT  Hybrid EMT  Late EMT
nearest                                                                                                                                                                                      
V1               0.0                          0.0                             0.0                 0.0                 0.1                 0.0                      64.8        33.7       1.4
V2               0.0                          0.0                             0.0                 0.0                 0.2                 0.0                       7.2        30.2      62.4
V3              26.2                          7.6                            34.0                24.3                 1.9                 5.8                       0.0         0.2       0.0
V4               5.2                          0.0                             0.0                 2.1                24.6                 4.9                       0.2        13.6      49.4

## Specialists: vertex x phase (% of row)

phase      G1   G2M     S
nearest                  
V1       20.5  46.3  33.1
V2       79.7   9.2  11.1
V3       51.6  26.7  21.6
V4       85.5   5.0   9.6

## Specialists: vertex x assay (% of row)

assay    multiome  scRNA
nearest                 
V1            0.5   99.5
V2            0.0  100.0
V3           26.9   73.1
V4           42.1   57.9

## Specialists: vertex x tumour (% of row)

sample   MU3  MU4   MU7   MU8  SC15  SC2b1  SC2b3  SC32  SC33  SC35  SC43  SC46
nearest                                                                        
V1       0.1  0.1   0.2   0.1   2.8    0.0   13.8   3.5   1.0  37.3  35.1   6.0
V2       0.0  0.0   0.0   0.0   0.1   99.7    0.2   0.0   0.0   0.0   0.0   0.0
V3       3.8  5.5   1.4  16.2   4.6    0.0    4.8  32.6  30.2   0.2   0.3   0.4
V4       2.4  2.6  12.7  24.4  17.9    0.0   11.3   1.6   0.4  14.6   0.9  11.3

## Each EMT state: % of its cells at each vertex or interior

group                 V1    V2    V3    V4  interior
cluster_annotation                                  
Epithelial           0.0   0.0  33.0   2.1      64.9
Early-Hybrid         0.0   0.0  23.7   5.0      71.3
Hybrid              44.5   4.9   0.0   0.3      50.3
Late-EMT             2.2  11.7   0.0  17.0      69.1

## Each clusterU: % of its cells at each vertex or interior

group                             V1    V2    V3    V4  interior
clusterU                                                        
Epithelial                       0.0   0.0  28.4   2.5      69.1
Epithelial (Differentiated)      0.0   0.0  53.1   0.0      46.9
Proliferative Early Hybrid EMT   0.0   0.0  37.7   0.0      62.3
Early Hybrid EMT_1               0.0   0.0  27.0   1.0      72.0
Early Hybrid EMT_2               0.0   0.1   3.1  17.8      79.0
Early Hybrid EMT_3               0.0   0.0  18.9   7.1      74.1
Proliferative Hybrid EMT        54.3   3.1   0.0   0.2      42.3
Hybrid EMT                      16.2   7.5   0.2   6.2      70.0
Late EMT                         0.6  13.8   0.0  20.0      65.7

## Each tumour: % of its cells at each vertex or interior

group     V1    V2    V3    V4  interior
sample                                  
MU3      0.1   0.0  16.7   4.7      78.5
MU4      0.2   0.0  21.9   4.6      73.3
MU7      0.1   0.0   1.8   7.5      90.6
MU8      0.0   0.0  11.4   7.6      81.0
SC15     3.0   0.1  10.9  18.7      67.3
SC2b1    0.0  81.8   0.0   0.0      18.2
SC2b3    8.1   0.1   6.1   6.3      79.4
SC32     3.2   0.0  63.8   1.4      31.6
SC33     1.3   0.0  81.8   0.4      16.5
SC35    25.3   0.0   0.2   9.4      65.2
SC43    58.9   0.0   1.0   1.4      38.7
SC46     6.8   0.0   1.1  12.2      79.9

## Tumour check per vertex

        n_specialists top_tumour  top_share_pct  n_tumours_ge5pct
vertex                                                           
V1               1744       SC35           37.3                 4
V2                902      SC2b1           99.7                 1
V3               3736       SC32           32.6                 4
V4               1652        MU8           24.4                 6

## Mean scores per vertex (specialists) and interior

          KS_EMT  LKPR_Epithelial_andrea  LKPR_EarlyHybridEMT_andrea  LKPR_HybridEMT_andrea  LKPR_LateEMT_andrea  purity
group                                                                                                                   
V1         0.484                  -0.496                       0.192                  0.355                0.303   0.541
V2         0.464                  -0.474                       0.259                  0.256                0.549   0.625
V3        -0.418                   1.086                       1.458                 -0.094               -0.089   0.524
V4         0.260                  -0.049                       0.291                  0.346                0.431   0.468
interior  -0.008                   0.366                       0.805                  0.132                0.161   0.158

## Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)

                   label_type                         subtype  k_in_bin  fold_enrichment  q_value
archetype                                                                                        
V1                   clusterU        Proliferative Hybrid EMT      1513           7.2737      0.0
V1                     sample                            SC43       641           6.1573      0.0
V1         cluster_annotation                          Hybrid      2203           6.1056      0.0
V1                     sample                            SC35       786           3.0535      0.0
V1                   clusterU                      Hybrid EMT       838           2.3053      0.0
V1                      phase                             G2M      1048           2.2595      0.0
V1                      phase                               S       796           1.8277      0.0
V1                     sample                            SC46       266           1.7511      0.0
V1                     sample                           SC2b3       518           1.7493      0.0
V1                      assay                           scRNA      2380           1.6928      0.0
V2                     sample                           SC2b1      1018           9.2626      0.0
V2                     sample                           SC2b3      1047           3.5358      0.0
V2                   clusterU                        Late EMT      1378           3.3699      0.0
V2         cluster_annotation                        Late-EMT      1886           3.0320      0.0
V2                   clusterU                      Hybrid EMT       823           2.2640      0.0
V2                      assay                           scRNA      2396           1.7042      0.0
V3                     sample                            SC33       903           6.5480      0.0
V3                     sample                            SC32       913           4.7799      0.0
V3                   clusterU  Proliferative Early Hybrid EMT       994           2.9424      0.0
V3                   clusterU     Epithelial (Differentiated)       145           2.7102      0.0
V3         cluster_annotation                      Epithelial       841           2.0284      0.0
V3                   clusterU                      Epithelial       600           1.7416      0.0
V3                   clusterU              Early Hybrid EMT_1       549           1.6334      0.0
V3                      phase                             G2M       752           1.6213      0.0
V3         cluster_annotation                    Early-Hybrid      1555           1.5573      0.0
V4                   clusterU                        Late EMT      1283           3.1376      0.0
V4         cluster_annotation                        Late-EMT      1584           2.5465      0.0
V4                   clusterU              Early Hybrid EMT_2       572           2.4966      0.0
V4                     sample                            SC15       291           1.8464      0.0
V4                     sample                            SC35       395           1.5345      0.0

## Top 15 genes higher near each vertex (log-expression difference)

       V1_up   V2_up   V3_up    V4_up
0     Lgals1  Col1a1   Krt14   Col5a2
1      Csrp2  Col1a2    Krt5      Fn1
2     Tm4sf1  Col3a1   Krt6a   Col1a1
3        Vim     Dcn  Krtdap   Col3a1
4     S100a4     Bgn   Krt17   Malat1
5        Dcn     Fn1   Fabp5   Col1a2
6      Timp1  Crabp1     Sfn      Bgn
7       Pfn1   Postn    Gclm    Hmcn1
8      Actg1  Lgals1   Fxyd3    Postn
9      H2afz     Vim    Perp  Col12a1
10      Actb  Cxcl12  Lgals7    Fgfr1
11       Ran   Aebp1    Rbm3     Fbn1
12    Crabp1   Sparc    Aqp3   Col6a1
13  Ppp1r14b   Crip1    Ly6a       Cp
14     Stmn2  Tmsb10   Pitx1      Dcn
