# Simplex description: gene set all, tag main, k=3

Cells: 19930; inside simplex 61.5%; specialists 71.1%; interior 28.9%
Fit cells inside simplex: 62.9%

## Specialists: vertex x EMT state (counts)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                         195            54           26           49
V2                        2153          3523         3832         3935
V3                         293            55           24           31

## Specialists: vertex x EMT state (% of row)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                        60.2          16.7          8.0         15.1
V2                        16.0          26.2         28.5         29.3
V3                        72.7          13.6          6.0          7.7

## Specialists: vertex x clusterU (% of row)

clusterU    K0    K1    K2    K3    K4
nearest                               
V1         1.2  17.6  68.5   7.1   5.6
V2        34.2   1.4  14.3  18.0  32.1
V3         2.7  58.1  32.5   2.0   4.7

## Specialists: vertex x phase (% of row)

phase      G1   G2M     S
nearest                  
V1       63.9  12.0  24.1
V2       65.0  15.0  20.0
V3       84.6   5.2  10.2

## Specialists: vertex x assay (% of row)

assay    in_vivo
nearest         
V1         100.0
V2         100.0
V3         100.0

## Specialists: vertex x tumour (% of row)

sample    4T1     E   EM1   EM2   EM3     M
nearest                                    
V1       39.8  15.7  19.4   6.2   5.9  13.0
V2       28.8   0.1  11.5  19.8   9.9  30.0
V3       15.1  56.8   8.9   8.7  10.2   0.2

## Each EMT state: % of its cells at each vertex or interior

group                V1    V2   V3  interior
cluster_annotation                          
Epithelial          3.9  43.2  5.9      47.0
Early-Hybrid        1.1  70.7  1.1      27.1
Late-Hybrid         0.5  76.9  0.5      22.1
Mesenchymal         1.0  79.0  0.6      19.4

## Each clusterU: % of its cells at each vertex or interior

group      V1    V2    V3  interior
clusterU                           
K0        0.1  77.7   0.2      22.0
K1        3.9  13.3  16.1      66.7
K2        6.1  52.6   3.6      37.7
K3        0.6  66.1   0.2      33.0
K4        0.3  82.1   0.4      17.2

## Each tumour: % of its cells at each vertex or interior

group    V1    V2    V3  interior
sample                           
4T1     2.1  64.3   1.0      32.5
E       6.4   1.5  28.8      63.2
EM1     2.2  53.1   1.2      43.5
EM2     0.6  78.7   1.0      19.7
EM3     1.0  67.0   2.1      30.0
M       0.9  83.1   0.0      16.0

## Tumour check per vertex

        n_specialists top_tumour  top_share_pct  n_tumours_ge5pct
vertex                                                           
V1                324        4T1           39.8                 6
V2              13443          M           30.0                 5
V3                403          E           56.8                 5

## Mean scores per vertex (specialists) and interior

          KS_EMT  purity
group                   
V1        -0.246   0.509
V2         0.075   0.565
V3        -0.466   0.485
interior  -0.129   0.213

## Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)

                   label_type      subtype  k_in_bin  fold_enrichment  q_value
archetype                                                                     
V1                     sample            E       323           4.0680      0.0
V1                   clusterU           K1       469           3.2300      0.0
V1                   clusterU           K2       958           2.6297      0.0
V1         cluster_annotation   Epithelial      1101           2.2095      0.0
V1                     sample          EM1       463           1.5894      0.0
V2                      phase          G2M       530           2.1712      0.0
V2                     sample            M      1022           2.1094      0.0
V2                      phase            S       701           2.1076      0.0
V2                   clusterU           K4       860           1.6368      0.0
V2         cluster_annotation  Mesenchymal       810           1.6255      0.0
V3                     sample            E       508           6.3980      0.0
V3                   clusterU           K1       652           4.4904      0.0
V3         cluster_annotation   Epithelial      1118           2.2436      0.0
V3                   clusterU           K2       597           1.6388      0.0

## Top 15 genes higher near each vertex (log-expression difference)

      V1_up     V2_up     V3_up
0       Fos      Bcl3     Cldn4
1      Atf3     Ccnd1      Csf3
2       Jun     Stmn1      Vmp1
3     Dusp1     Top2a     Lamc2
4     Cldn4      Odc1      Apoe
5      Apoe     Cebpd     Plet1
6      Ccn1     Galk1      Ubp1
7      Fosb      Pfn1     Kcnn4
8     Zfp36      Has2    Col7a1
9      Klf2     Eif5a    Arglu1
10   S100a8    Lgals1     Epcam
11     Cd74  Hist1h1e    Spaca6
12  Gadd45b     Htra1  Arhgap42
13  Gadd45g    S100a4      Cdh1
14   Hspa1b      Cdk1      Plec
