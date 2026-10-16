# Simplex description: gene set ks, tag main, k=4

Cells: 30590; inside simplex 36.0%; specialists 46.5%; interior 53.5%
Fit cells inside simplex: 36.3%

## Specialists: vertex x EMT state (counts)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                         167            21            2            0
V2                         984            37            5            0
V3                           0             8          194         2989
V4                        1292          3906         3451         1170

## Specialists: vertex x EMT state (% of row)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                        87.9          11.1          1.1          0.0
V2                        95.9           3.6          0.5          0.0
V3                         0.0           0.3          6.1         93.7
V4                        13.2          39.8         35.1         11.9

## Specialists: vertex x clusterU (% of row)

clusterU    K0    K1     K2    K3    K4
nearest                                
V1         0.0   1.1    0.0  98.9   0.0
V2         0.0  99.5    0.0   0.0   0.5
V3         0.0   0.0  100.0   0.0   0.0
V4        55.1   0.0    0.0   0.0  44.8

## Specialists: vertex x phase (% of row)

phase      G1   G2M     S
nearest                  
V1       88.9   2.1   8.9
V2       98.5   0.4   1.1
V3       57.3  21.5  21.2
V4       51.9  12.0  36.0

## Specialists: vertex x cell line (% of row)

assay    HuMEC  MCF10A
nearest               
V1         1.1    98.9
V2        94.0     6.0
V3         0.1    99.9
V4       100.0     0.0

## Specialists: vertex x position (% of row)

position  inner  outer
nearest               
V1         86.3   13.7
V2         52.3   47.7
V3         35.0   65.0
V4         56.7   43.3

## Specialists: vertex x arm (% of row)

arm      Mock  TGFB
nearest            
V1       98.9   1.1
V2       52.3  47.7
V3        4.5  95.5
V4       55.9  44.1

## Specialists: vertex x tumour (% of row)

sample   HuMEC_Mock  HuMEC_TGFB  MCF10A_Mock  MCF10A_TGFB
nearest                                                  
V1              1.1         0.0         97.9          1.1
V2             50.7        43.3          1.7          4.4
V3              0.0         0.1          4.5         95.5
V4             55.8        44.1          0.0          0.0

## Each EMT state: % of its cells at each vertex or interior

group                V1    V2    V3    V4  interior
cluster_annotation                                 
Epithelial          2.2  12.9   0.0  16.9      68.1
Early-Hybrid        0.3   0.5   0.1  51.1      48.1
Late-Hybrid         0.0   0.1   2.5  45.1      52.2
Mesenchymal         0.0   0.0  39.1  15.3      45.6

## Each clusterU: % of its cells at each vertex or interior

group      V1    V2    V3    V4  interior
clusterU                                 
K0        0.0   0.0   0.0  53.8      46.2
K1        0.1  26.2   0.0   0.0      73.8
K2        0.0   0.0  63.8   0.0      36.2
K3        5.4   0.0   0.0   0.1      94.5
K4        0.0   0.1   0.0  54.2      45.8

## Each sample: % of its cells at each vertex or interior

group         V1   V2    V3    V4  interior
sample                                     
HuMEC_Mock   0.0  4.6   0.0  48.9      46.5
HuMEC_TGFB   0.0  4.4   0.0  42.7      52.9
MCF10A_Mock  3.7  0.3   2.8   0.1      93.0
MCF10A_TGFB  0.0  1.1  72.0   0.0      26.9

## Each position: % of its cells at each vertex or interior

group      V1   V2    V3    V4  interior
position                                
inner     1.1  3.7   7.6  37.8      49.8
outer     0.2  3.1  13.1  26.8      56.9

## Each arm: % of its cells at each vertex or interior

group   V1   V2    V3    V4  interior
arm                                  
Mock   1.2  3.3   0.9  33.9      60.8
TGFB   0.0  3.4  21.2  30.1      45.3

## Tumour check per vertex

        n_specialists   top_tumour  top_share_pct  n_tumours_ge5pct
vertex                                                             
V1                190  MCF10A_Mock           97.9                 1
V2               1026   HuMEC_Mock           50.7                 2
V3               3191  MCF10A_TGFB           95.5                 1
V4               9819   HuMEC_Mock           55.8                 2

## Mean scores per vertex (specialists) and interior

          KS_EMT  purity
group                   
V1        -0.465   0.479
V2        -0.693   0.502
V3         0.618   0.566
V4        -0.004   0.534
interior  -0.069   0.167

## Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)

                   label_type       subtype  k_in_bin  fold_enrichment  q_value
archetype                                                                      
V1                   clusterU            K3      2226           6.3655      0.0
V1                     sample   MCF10A_Mock      2723           5.4558      0.0
V1                      assay        MCF10A      2790           3.0257      0.0
V1         cluster_annotation    Epithelial      1401           1.8319      0.0
V1                        arm          Mock      2937           1.8126      0.0
V1                   clusterU            K1       652           1.6701      0.0
V2                   clusterU            K1      2843           7.2823      0.0
V2         cluster_annotation    Epithelial      2830           3.7003      0.0
V2                     sample    HuMEC_Mock      1870           1.6679      0.0
V2                      phase            G1      2953           1.5430      0.0
V3                     sample   MCF10A_TGFB      2749           6.4988      0.0
V3                   clusterU            K2      3059           6.1180      0.0
V3         cluster_annotation   Mesenchymal      2928           3.8285      0.0
V3                      assay        MCF10A      3057           3.3153      0.0
V3                      phase           G2M       691           2.0187      0.0
V3                        arm          TGFB      2751           1.9121      0.0
V4                   clusterU            K0      2217           2.2031      0.0
V4         cluster_annotation  Early-Hybrid      1277           1.6699      0.0
V4         cluster_annotation   Late-Hybrid      1238           1.6189      0.0
V4                     sample    HuMEC_TGFB      1614           1.5891      0.0

## Top 15 genes higher near each vertex (log-expression difference)

            V1_up    V2_up    V3_up    V4_up
0           KRT81   SPRR1A    KRT81    KRT15
1             PI3    KRT6B     MT2A   CALML3
2           S100P   SPRR1B     MT1E   PYCARD
3            PSCA    LYPD3    TGFBI    KRT14
4            SAA1    KRT16      FN1   S100A9
5            SBSN    FABP5    CXCL3     GLUL
6            SLPI     SBSN    CXCL1    FXYD3
7   RP11-295G20.2     SLPI    CXCL2     CBR1
8           DUSP1  TACSTD2  ANGPTL4   S100A8
9           UQCRQ    KRT6C   DUSP23     KRT5
10           CSTA   CRABP2     MT1X      CA2
11           SAA2     CES2    POSTN   PDLIM1
12          COX6C    CLIC3    LAMC2  AKR1B10
13          CLDN7      PI3     DKK1      TKT
14          SRRM2  SULT2B1   PMEPA1   FGFBP1
