# Simplex description: gene set all, tag scrna6, k=3

Cells: 23959; inside simplex 45.1%; specialists 57.7%; interior 42.3%
Fit cells inside simplex: 60.0%

## Specialists: vertex x EMT state (counts)

cluster_annotation  Epithelial  Early-Hybrid  Hybrid  Late-EMT
nearest                                                       
V1                           0            14      16      2959
V2                        2984          6512       7        15
V3                           2             0    1301        10

## Specialists: vertex x EMT state (% of row)

cluster_annotation  Epithelial  Early-Hybrid  Hybrid  Late-EMT
nearest                                                       
V1                         0.0           0.5     0.5      99.0
V2                        31.4          68.4     0.1       0.2
V3                         0.2           0.0    99.1       0.8

## Specialists: vertex x clusterU (% of row)

clusterU  Epithelial  Epithelial (Differentiated)  Proliferative Early Hybrid EMT  Early Hybrid EMT_1  Early Hybrid EMT_2  Early Hybrid EMT_3  Proliferative Hybrid EMT  Hybrid EMT  Late EMT
nearest                                                                                                                                                                                      
V1               0.0                          0.0                             0.0                 0.0                 0.3                 0.0                       0.2         9.1      90.4
V2              25.5                          5.1                            27.7                29.8                 5.4                 6.3                       0.1         0.2       0.0
V3               0.0                          0.0                             0.2                 0.0                 0.0                 0.0                      79.3        20.6       0.0

## Specialists: vertex x phase (% of row)

phase      G1   G2M     S
nearest                  
V1       90.4   2.1   7.4
V2       52.9  25.0  22.1
V3       11.9  50.8  37.3

## Specialists: vertex x assay (% of row)

assay    multiome  scRNA
nearest                 
V1            2.1   97.9
V2           68.4   31.6
V3            0.0  100.0

## Specialists: vertex x tumour (% of row)

sample   MU3  MU4   MU7   MU8  SC15  SC2b1  SC2b3  SC32  SC33  SC35  SC43  SC46
nearest                                                                        
V1       0.3  1.1   0.3   0.4   0.7   17.1   43.5   0.4   0.1  27.8   2.5   6.0
V2       7.1  7.9  13.7  39.7   5.2    0.0    1.9  12.8  10.8   0.1   0.1   0.7
V3       0.0  0.0   0.0   0.0   1.5    8.2   22.5   2.9   0.7  27.3  32.6   4.2

## Each EMT state: % of its cells at each vertex or interior

group                 V1    V2    V3  interior
cluster_annotation                            
Epithelial           0.0  72.0   0.0      28.0
Early-Hybrid         0.1  65.2   0.0      34.6
Hybrid               0.4   0.2  36.1      63.3
Late-EMT            47.6   0.2   0.2      52.0

## Each clusterU: % of its cells at each vertex or interior

group                             V1    V2    V3  interior
clusterU                                                  
Epithelial                       0.0  70.4   0.0      29.6
Epithelial (Differentiated)      0.0  90.1   0.0       9.9
Proliferative Early Hybrid EMT   0.0  78.0   0.1      22.0
Early Hybrid EMT_1               0.0  84.4   0.0      15.6
Early Hybrid EMT_2               0.4  22.4   0.0      77.2
Early Hybrid EMT_3               0.0  52.4   0.0      47.6
Proliferative Hybrid EMT         0.3   0.2  50.0      49.4
Hybrid EMT                       7.5   0.6   7.4      84.5
Late EMT                        66.1   0.0   0.0      33.9

## Each tumour: % of its cells at each vertex or interior

group     V1    V2    V3  interior
sample                            
MU3      1.2  79.0   0.0      19.8
MU4      3.5  79.8   0.0      16.7
MU7      0.3  46.4   0.0      53.3
MU8      0.2  71.3   0.0      28.5
SC15     1.3  31.5   1.3      65.9
SC2b1   46.4   0.1   9.8      43.7
SC2b3   43.9   6.2  10.0      39.9
SC32     0.6  64.0   2.0      33.4
SC33     0.1  74.3   0.7      24.9
SC35    32.2   0.3  13.9      53.5
SC43     7.1   0.8  41.1      51.0
SC46    11.8   4.4   3.6      80.2

## Tumour check per vertex

        n_specialists top_tumour  top_share_pct  n_tumours_ge5pct
vertex                                                           
V1               2989      SC2b3           43.5                 4
V2               9518        MU8           39.7                 7
V3               1313       SC43           32.6                 4

## Mean scores per vertex (specialists) and interior

          KS_EMT  LKPR_Epithelial_andrea  LKPR_EarlyHybridEMT_andrea  LKPR_HybridEMT_andrea  LKPR_LateEMT_andrea  purity
group                                                                                                                   
V1         0.524                  -0.384                       0.205                  0.405                0.691   0.529
V2        -0.390                   1.021                       1.305                 -0.096               -0.087   0.547
V3         0.475                  -0.492                       0.304                  0.334                0.267   0.507
interior   0.150                   0.058                       0.580                  0.240                0.235   0.191

## Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)

                   label_type                         subtype  k_in_bin  fold_enrichment  q_value
archetype                                                                                        
V1                   clusterU                        Late EMT      2255           5.5146      0.0
V1                     sample                           SC2b1       502           4.5676      0.0
V1         cluster_annotation                        Late-EMT      2379           3.8246      0.0
V1                     sample                           SC2b3      1095           3.6979      0.0
V1                     sample                            SC35       590           2.2921      0.0
V1                      assay                           scRNA      2345           1.6679      0.0
V2                     sample                            SC33       442           3.2051      0.0
V2                   clusterU     Epithelial (Differentiated)       169           3.1587      0.0
V2                     sample                            SC32       521           2.7276      0.0
V2                   clusterU  Proliferative Early Hybrid EMT       877           2.5961      0.0
V2                   clusterU              Early Hybrid EMT_1       719           2.1392      0.0
V2                     sample                             MU4       196           2.0872      0.0
V2                      phase                             G2M       817           1.7615      0.0
V2         cluster_annotation                      Epithelial       724           1.7462      0.0
V2         cluster_annotation                    Early-Hybrid      1672           1.6744      0.0
V2                     sample                             MU8       858           1.6164      0.0
V2                     sample                             MU3       130           1.5258      0.0
V2                   clusterU                      Epithelial       523           1.5181      0.0
V3                   clusterU        Proliferative Hybrid EMT      1613           7.7545      0.0
V3         cluster_annotation                          Hybrid      2247           6.2276      0.0
V3                     sample                            SC43       586           5.6290      0.0
V3                     sample                            SC35       669           2.5990      0.0
V3                      phase                             G2M      1123           2.4212      0.0
V3                   clusterU                      Hybrid EMT       719           1.9779      0.0
V3                      phase                               S       812           1.8644      0.0
V3                     sample                           SC2b1       204           1.8562      0.0
V3                     sample                           SC2b3       530           1.7899      0.0
V3                      assay                           scRNA      2387           1.6978      0.0
V3                     sample                            SC46       232           1.5273      0.0

## Top 15 genes higher near each vertex (log-expression difference)

     V1_up    V2_up   V3_up
0   Col3a1    Krt14  Lgals1
1   Col1a1     Krt5   Csrp2
2   Col1a2    Krt6a  Tm4sf1
3    Postn   Krtdap     Vim
4      Dcn    Fabp5  S100a4
5      Bgn      Sfn   H2afz
6      Fn1    Krt17    Pfn1
7   Col5a2     Gclm     Dcn
8     Mmp2     Perp   Actg1
9    Sparc      Dsp   Timp1
10   Aebp1     Nrg1     Ran
11  Cxcl12  Col17a1   Hmgb2
12  Crabp1    Pitx1  Crabp1
13   Fstl1    Itga6    Actb
14     Id3    Fxyd3  Ranbp1
