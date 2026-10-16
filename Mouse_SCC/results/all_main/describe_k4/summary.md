# Simplex description: gene set all, tag main, k=4

Cells: 23959; inside simplex 31.3%; specialists 34.6%; interior 65.4%
Fit cells inside simplex: 30.5%

## Specialists: vertex x EMT state (counts)

cluster_annotation  Epithelial  Early-Hybrid  Hybrid  Late-EMT
nearest                                                       
V1                         951          1838       0         0
V2                           1             7    2005       177
V3                         421           632       0        25
V4                           0            18      57      2169

## Specialists: vertex x EMT state (% of row)

cluster_annotation  Epithelial  Early-Hybrid  Hybrid  Late-EMT
nearest                                                       
V1                        34.1          65.9     0.0       0.0
V2                         0.0           0.3    91.6       8.1
V3                        39.1          58.6     0.0       2.3
V4                         0.0           0.8     2.5      96.7

## Specialists: vertex x clusterU (% of row)

clusterU  Epithelial  Epithelial (Differentiated)  Proliferative Early Hybrid EMT  Early Hybrid EMT_1  Early Hybrid EMT_2  Early Hybrid EMT_3  Proliferative Hybrid EMT  Hybrid EMT  Late EMT
nearest                                                                                                                                                                                      
V1              23.2                          8.9                            29.8                30.5                 1.7                 6.0                       0.0         0.0       0.0
V2               0.1                          0.0                             0.1                 0.0                 0.0                 0.0                      64.0        34.0       1.8
V3              38.1                          0.2                             0.5                16.3                35.9                 6.9                       0.0         1.8       0.4
V4               0.0                          0.0                             0.0                 0.0                 0.6                 0.0                       0.4        15.6      83.4

## Specialists: vertex x phase (% of row)

phase      G1   G2M     S
nearest                  
V1       58.6  23.3  18.1
V2       20.7  45.7  33.6
V3       81.4   7.2  11.4
V4       89.0   2.6   8.4

## Specialists: vertex x assay (% of row)

assay    multiome  scRNA
nearest                 
V1           23.9   76.1
V2            1.1   98.9
V3           66.7   33.3
V4            0.0  100.0

## Specialists: vertex x tumour (% of row)

sample   MU3  MU4   MU7   MU8  SC15  SC2b1  SC2b3  SC32  SC33  SC35  SC43  SC46
nearest                                                                        
V1       4.7  4.5   0.6  14.1   2.5    0.8    5.3  35.6  31.2   0.3   0.3   0.1
V2       0.1  0.2   0.5   0.3   3.7    0.1   18.9   3.5   1.0  35.2  27.9   8.7
V3       1.3  0.6  53.0  11.8  31.7    0.0    0.0   0.1   0.0   0.0   0.0   1.5
V4       0.0  0.0   0.0   0.0   0.1   34.4   45.5   0.4   0.3  16.2   0.8   2.4

## Each EMT state: % of its cells at each vertex or interior

group                 V1    V2    V3    V4  interior
cluster_annotation                                  
Epithelial          22.9   0.0  10.2   0.0      66.9
Early-Hybrid        18.4   0.1   6.3   0.2      75.0
Hybrid               0.0  55.6   0.0   1.6      42.8
Late-EMT             0.0   2.8   0.4  34.9      61.9

## Each clusterU: % of its cells at each vertex or interior

group                             V1    V2    V3    V4  interior
clusterU                                                        
Epithelial                      18.8   0.1  11.9   0.0      69.2
Epithelial (Differentiated)     46.2   0.0   0.4   0.0      53.5
Proliferative Early Hybrid EMT  24.6   0.1   0.1   0.0      75.2
Early Hybrid EMT_1              25.3   0.0   5.2   0.0      69.5
Early Hybrid EMT_2               2.1   0.0  16.9   0.6      80.4
Early Hybrid EMT_3              14.6   0.0   6.5   0.1      78.9
Proliferative Hybrid EMT         0.0  67.4   0.0   0.4      32.2
Hybrid EMT                       0.0  20.5   0.5   9.6      69.4
Late EMT                         0.0   1.0   0.1  45.8      53.2

## Each tumour: % of its cells at each vertex or interior

group     V1    V2    V3    V4  interior
sample                                  
MU3     15.3   0.2   1.6   0.0      82.9
MU4     13.3   0.4   0.7   0.1      85.4
MU7      0.6   0.4  20.4   0.0      78.6
MU8      7.4   0.1   2.4   0.0      90.1
SC15     4.4   5.1  21.7   0.1      68.6
SC2b1    2.0   0.3   0.0  70.2      27.5
SC2b3    5.0  13.9   0.0  34.4      46.6
SC32    52.0   4.0   0.1   0.5      43.5
SC33    63.0   1.5   0.0   0.5      35.0
SC35     0.3  30.0   0.0  14.1      55.6
SC43     0.9  58.7   0.0   1.6      38.8
SC46     0.3  12.6   1.1   3.5      82.6

## Tumour check per vertex

        n_specialists top_tumour  top_share_pct  n_tumours_ge5pct
vertex                                                           
V1               2789       SC32           35.6                 4
V2               2190       SC35           35.2                 4
V3               1078        MU7           53.0                 3
V4               2244      SC2b3           45.5                 3

## Mean scores per vertex (specialists) and interior

          KS_EMT  LKPR_Epithelial_andrea  LKPR_EarlyHybridEMT_andrea  LKPR_HybridEMT_andrea  LKPR_LateEMT_andrea  purity
group                                                                                                                   
V1        -0.443                   1.104                       1.475                 -0.103               -0.095   0.499
V2         0.478                  -0.482                       0.224                  0.360                0.300   0.574
V3        -0.272                   0.620                       0.710                 -0.054               -0.081   0.482
V4         0.500                  -0.398                       0.245                  0.373                0.663   0.567
interior  -0.041                   0.429                       0.856                  0.121                0.139   0.162

## Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)

                   label_type                         subtype  k_in_bin  fold_enrichment  q_value
archetype                                                                                        
V1                     sample                            SC33       984           7.1353      0.0
V1                     sample                            SC32      1026           5.3715      0.0
V1                   clusterU  Proliferative Early Hybrid EMT       923           2.7323      0.0
V1                   clusterU     Epithelial (Differentiated)       116           2.1681      0.0
V1         cluster_annotation                      Epithelial       797           1.9223      0.0
V1                   clusterU              Early Hybrid EMT_1       618           1.8387      0.0
V1                   clusterU                      Epithelial       592           1.7184      0.0
V1                      assay                           scRNA      2269           1.6138      0.0
V1         cluster_annotation                    Early-Hybrid      1599           1.6013      0.0
V2                   clusterU        Proliferative Hybrid EMT      1478           7.1055      0.0
V2         cluster_annotation                          Hybrid      2204           6.1084      0.0
V2                     sample                            SC43       624           5.9940      0.0
V2                     sample                            SC35       830           3.2244      0.0
V2                   clusterU                      Hybrid EMT       874           2.4043      0.0
V2                      phase                             G2M      1059           2.2832      0.0
V2                     sample                           SC2b3       521           1.7595      0.0
V2                      phase                               S       760           1.7450      0.0
V2                      assay                           scRNA      2390           1.6999      0.0
V3                     sample                             MU7      1035           3.6950      0.0
V3                   clusterU              Early Hybrid EMT_2       652           2.8458      0.0
V3                     sample                            SC15       445           2.8235      0.0
V3                   clusterU                      Epithelial       807           2.3424      0.0
V3         cluster_annotation                      Epithelial       829           1.9994      0.0
V3                      assay                        multiome      1922           1.9413      0.0
V3                   clusterU              Early Hybrid EMT_1       635           1.8892      0.0
V3                   clusterU              Early Hybrid EMT_3       173           1.5109      0.0
V4                     sample                           SC2b1       818           7.4428      0.0
V4                   clusterU                        Late EMT      1938           4.7393      0.0
V4         cluster_annotation                        Late-EMT      2286           3.6751      0.0
V4                     sample                           SC2b3      1025           3.4615      0.0
V4                     sample                            SC35       466           1.8103      0.0
V4                      assay                           scRNA      2396           1.7042      0.0

## Top 15 genes higher near each vertex (log-expression difference)

            V1_up   V2_up     V3_up   V4_up
0           Krt14  Lgals1    Camk1d  Col3a1
1            Krt5   Csrp2     Pcdh7  Col1a1
2           Krt17  Tm4sf1     Maml2  Col1a2
3           Krt6a     Vim    Frmd4b     Dcn
4          Krtdap     Dcn  AY036118   Postn
5             Sfn  S100a4    Adgrl3     Bgn
6           Fabp5   Timp1     Pard3     Fn1
7           Fxyd3    Pfn1     Magi1  Crabp1
8            Gclm   Actg1      Ext1   Sparc
9          Lgals7   H2afz    Zswim6   Aebp1
10           Rbm3  Crabp1     Mast4    Mmp2
11           Perp    Actb     Aopep  Col5a2
12           Aqp3     Ran     Atxn1  Cxcl12
13           Klk8  Col1a2    Ptpn13     Gsn
14  2200002D01Rik   Stmn2    Hivep2     Id3
