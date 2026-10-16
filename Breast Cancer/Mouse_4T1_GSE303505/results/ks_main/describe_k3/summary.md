# Simplex description: gene set ks, tag main, k=3

Cells: 19930; inside simplex 85.5%; specialists 47.6%; interior 52.4%
Fit cells inside simplex: 84.1%

## Specialists: vertex x EMT state (counts)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                         210             0            0            0
V2                        1423          2446         2286         1764
V3                         444           294          259          355

## Specialists: vertex x EMT state (% of row)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                       100.0           0.0          0.0          0.0
V2                        18.0          30.9         28.9         22.3
V3                        32.8          21.7         19.2         26.3

## Specialists: vertex x clusterU (% of row)

clusterU    K0     K1    K2    K3    K4
nearest                                
V1         0.0  100.0   0.0   0.0   0.0
V2        61.0    0.0   0.5   0.5  38.0
V3         0.0    0.0  54.1  45.6   0.3

## Specialists: vertex x phase (% of row)

phase      G1   G2M     S
nearest                  
V1       95.2   0.5   4.3
V2       59.0  17.1  23.9
V3       88.3   3.7   8.0

## Specialists: vertex x assay (% of row)

assay    in_vivo
nearest         
V1         100.0
V2         100.0
V3         100.0

## Specialists: vertex x tumour (% of row)

sample    4T1     E   EM1   EM2   EM3     M
nearest                                    
V1       13.8  68.6   3.8  11.0   0.0   2.9
V2       41.4   0.0   7.1  28.0  11.7  11.8
V3       10.4   1.7  31.4   5.2   6.1  45.3

## Each EMT state: % of its cells at each vertex or interior

group                V1    V2   V3  interior
cluster_annotation                          
Epithelial          4.2  28.6  8.9      58.3
Early-Hybrid        0.0  49.1  5.9      45.0
Late-Hybrid         0.0  45.9  5.2      48.9
Mesenchymal         0.0  35.4  7.1      57.5

## Each clusterU: % of its cells at each vertex or interior

group       V1    V2    V3  interior
clusterU                            
K0         0.0  81.7   0.0      18.3
K1        14.5   0.1   0.0      85.5
K2         0.0   1.1  20.1      78.8
K3         0.0   1.0  16.8      82.2
K4         0.0  57.3   0.1      42.7

## Each tumour: % of its cells at each vertex or interior

group     V1    V2    V3  interior
sample                            
4T1      0.5  54.5   2.3      42.7
E       18.1   0.0   2.9      79.0
EM1      0.3  19.4  14.6      65.8
EM2      0.7  65.6   2.1      31.7
EM3      0.0  46.5   4.1      49.4
M        0.1  19.3  12.7      67.9

## Tumour check per vertex

        n_specialists top_tumour  top_share_pct  n_tumours_ge5pct
vertex                                                           
V1                210          E           68.6                 3
V2               7919        4T1           41.4                 5
V3               1352          M           45.3                 5

## Mean scores per vertex (specialists) and interior

          KS_EMT  purity
group                   
V1        -1.451   0.497
V2         0.047   0.552
V3         0.007   0.505
interior  -0.008   0.171

## Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)

                   label_type      subtype  k_in_bin  fold_enrichment  q_value
archetype                                                                     
V1                   clusterU           K1      1448           9.9725      0.0
V1                     sample            E       745           9.3829      0.0
V1         cluster_annotation   Epithelial      1740           3.4919      0.0
V2                     sample          EM3       460           2.3162      0.0
V2                   clusterU           K4      1124           2.1393      0.0
V2                     sample          EM2       661           1.9550      0.0
V2                      phase            S       511           1.5364      0.0
V3                   clusterU           K2       984           2.7011      0.0
V3                   clusterU           K3       940           2.5627      0.0
V3                     sample            M      1140           2.3529      0.0
V3         cluster_annotation  Mesenchymal       753           1.5111      0.0

## Top 15 genes higher near each vertex (log-expression difference)

     V1_up   V2_up     V3_up
0    Cldn4    Areg      Lcn2
1    Epcam  S100a4    Ctla2a
2     Cdh1    Odc1      Ccn1
3    Lamc2    Slpi     Thbs1
4      Jup   Ptgs2      Fosb
5    Cldn7    Csf3        Cp
6     Krt7    Ier3    Igfbp7
7     St14   Procr       Mgp
8    Cldn3   Cxcl1       Fos
9    Krt14    Has2      Klf9
10     Dsp     Sfn  Serpine2
11     Sfn  Errfi1    Pdgfra
12   Plet1   Plet1    S100a8
13   Lypd3     Lif      Ccn2
14  Hspa1b    Npnt     Spry1
