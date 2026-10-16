# Simplex description: gene set all, tag main, k=4

Cells: 19930; inside simplex 35.7%; specialists 42.8%; interior 57.2%
Fit cells inside simplex: 36.6%

## Specialists: vertex x EMT state (counts)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                         267            46           22           32
V2                         159            41           21           40
V3                         368           596          834         1232
V4                         959          1496         1393         1017

## Specialists: vertex x EMT state (% of row)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                        72.8          12.5          6.0          8.7
V2                        60.9          15.7          8.0         15.3
V3                        12.1          19.7         27.5         40.7
V4                        19.7          30.8         28.6         20.9

## Specialists: vertex x clusterU (% of row)

clusterU    K0    K1    K2    K3    K4
nearest                               
V1         3.0  58.6  31.1   1.1   6.3
V2         1.5  18.0  68.2   7.3   5.0
V3         1.2   0.9  31.0  49.5  17.4
V4        61.0   1.9   7.2   1.3  28.4

## Specialists: vertex x phase (% of row)

phase      G1   G2M     S
nearest                  
V1       86.4   4.6   9.0
V2       58.2  14.9  26.8
V3       96.8   1.4   1.8
V4       32.1  28.8  39.1

## Specialists: vertex x assay (% of row)

assay    in_vivo
nearest         
V1         100.0
V2         100.0
V3         100.0
V4         100.0

## Specialists: vertex x tumour (% of row)

sample    4T1     E   EM1   EM2   EM3     M
nearest                                    
V1       12.3  57.2   9.0   8.4  12.8   0.3
V2       48.7  18.0  14.9   6.9   5.7   5.7
V3        2.4   0.0  18.8   1.2   2.0  75.5
V4       43.3   0.2   6.4  31.8   9.9   8.4

## Each EMT state: % of its cells at each vertex or interior

group                V1   V2    V3    V4  interior
cluster_annotation                                
Epithelial          5.4  3.2   7.4  19.2      64.8
Early-Hybrid        0.9  0.8  12.0  30.0      56.3
Late-Hybrid         0.4  0.4  16.7  28.0      54.4
Mesenchymal         0.6  0.8  24.7  20.4      53.4

## Each clusterU: % of its cells at each vertex or interior

group       V1   V2    V3    V4  interior
clusterU                                 
K0         0.2  0.1   0.6  50.2      48.9
K1        14.8  3.2   1.9   6.5      73.6
K2         3.1  4.9  25.7   9.7      56.6
K3         0.1  0.5  40.9   1.8      56.7
K4         0.4  0.2  10.0  26.3      62.9

## Each tumour: % of its cells at each vertex or interior

group     V1   V2    V3    V4  interior
sample                                 
4T1      0.7  2.1   1.2  35.0      60.9
E       26.4  5.9   0.1   1.3      66.2
EM1      1.1  1.3  19.6  10.6      67.3
EM2      0.9  0.5   1.1  45.7      51.8
EM3      2.4  0.8   3.1  24.3      69.5
M        0.0  0.3  47.2   8.5      44.0

## Tumour check per vertex

        n_specialists top_tumour  top_share_pct  n_tumours_ge5pct
vertex                                                           
V1                367          E           57.2                 5
V2                261        4T1           48.7                 6
V3               3030          M           75.5                 2
V4               4865        4T1           43.3                 5

## Mean scores per vertex (specialists) and interior

          KS_EMT  purity
group                   
V1        -0.470   0.475
V2        -0.236   0.488
V3         0.132   0.555
V4         0.033   0.570
interior  -0.029   0.166

## Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)

                   label_type      subtype  k_in_bin  fold_enrichment  q_value
archetype                                                                     
V1                     sample            E       506           6.3728      0.0
V1                   clusterU           K1       650           4.4766      0.0
V1         cluster_annotation   Epithelial      1111           2.2296      0.0
V1                   clusterU           K2       580           1.5921      0.0
V2                     sample            E       396           4.9874      0.0
V2                   clusterU           K1       536           3.6915      0.0
V2                   clusterU           K2       974           2.6736      0.0
V2         cluster_annotation   Epithelial      1158           2.3239      0.0
V3                     sample            M      1835           3.7874      0.0
V3                   clusterU           K3       840           2.2901      0.0
V3         cluster_annotation  Mesenchymal       980           1.9667      0.0
V3                   clusterU           K2       659           1.8089      0.0
V4                      phase            S       959           2.8833      0.0
V4                      phase          G2M       702           2.8759      0.0
V4                     sample          EM2       756           2.2360      0.0
V4                   clusterU           K0      1271           2.1495      0.0

## Top 15 genes higher near each vertex (log-expression difference)

       V1_up   V2_up     V3_up      V4_up
0      Cldn4     Fos    Ctla2a     S100a4
1       Csf3    Apoe       Mgp      Top2a
2       Vmp1    Atf3      Fosb   Hist1h1e
3      Lamc2   Cldn4        Cp      Stmn1
4      Plet1   Dusp1      Tslp       Odc1
5       Apoe     Jun  Trp53i11      Ccnd1
6      Kcnn4    Cd74     Thbs1       Cdk1
7       Ubp1    Ccn1      Sgk1     Tuba1b
8     Col7a1  S100a8    Cited2      Mki67
9      Epcam   Zfp36     Nr4a1       Anln
10    Arglu1    Cd14      Skil  Hist1h2ap
11    Spaca6   H2-Aa      Rbp1       Mcm6
12      Plec  S100a9    Igfbp7       Areg
13      Cdh1  Hspa1b      Ier5       Vars
14  Arhgap42    Klf2     Rdh10      Fabp5
