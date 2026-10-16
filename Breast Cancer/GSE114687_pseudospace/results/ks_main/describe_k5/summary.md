# Simplex description: gene set ks, tag main, k=5

Cells: 30590; inside simplex 16.3%; specialists 28.8%; interior 71.2%
Fit cells inside simplex: 17.0%

## Specialists: vertex x EMT state (counts)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                        1976           135           19            0
V2                           4            19           50           82
V3                         485          1693         1416          392
V4                         281             4            0            0
V5                           0            40          377         1843

## Specialists: vertex x EMT state (% of row)

cluster_annotation  Epithelial  Early-Hybrid  Late-Hybrid  Mesenchymal
nearest                                                               
V1                        92.8           6.3          0.9          0.0
V2                         2.6          12.3         32.3         52.9
V3                        12.2          42.5         35.5          9.8
V4                        98.6           1.4          0.0          0.0
V5                         0.0           1.8         16.7         81.5

## Specialists: vertex x clusterU (% of row)

clusterU    K0    K1    K2    K3    K4
nearest                               
V1         0.1  93.5   0.0   0.0   6.4
V2        76.8   1.9  16.1   5.2   0.0
V3        59.6   0.0   0.0   0.3  40.2
V4         0.0  60.4   0.0  39.6   0.0
V5         0.0   0.0  99.0   1.0   0.0

## Specialists: vertex x phase (% of row)

phase      G1   G2M     S
nearest                  
V1       96.9   0.2   2.9
V2       85.2   2.6  12.3
V3       49.5  15.3  35.2
V4       93.0   2.5   4.6
V5       47.4  24.6  27.9

## Specialists: vertex x cell line (% of row)

assay    HuMEC  MCF10A
nearest               
V1       100.0     0.0
V2        72.9    27.1
V3       100.0     0.0
V4         1.1    98.9
V5         0.0   100.0

## Specialists: vertex x position (% of row)

position  inner  outer
nearest               
V1         64.9   35.1
V2         25.2   74.8
V3         60.8   39.2
V4         94.7    5.3
V5         21.5   78.5

## Specialists: vertex x arm (% of row)

arm      Mock  TGFB
nearest            
V1       67.3  32.7
V2        0.6  99.4
V3       57.4  42.6
V4       94.7   5.3
V5       39.6  60.4

## Specialists: vertex x tumour (% of row)

sample   HuMEC_Mock  HuMEC_TGFB  MCF10A_Mock  MCF10A_TGFB
nearest                                                  
V1             67.3        32.7          0.0          0.0
V2              0.6        72.3          0.0         27.1
V3             57.3        42.6          0.0          0.0
V4              1.1         0.0         93.7          5.3
V5              0.0         0.0         39.6         60.4

## Each EMT state: % of its cells at each vertex or interior

group                 V1   V2    V3   V4    V5  interior
cluster_annotation                                      
Epithelial          25.8  0.1   6.3  3.7   0.0      64.1
Early-Hybrid         1.8  0.2  22.1  0.1   0.5      75.3
Late-Hybrid          0.2  0.7  18.5  0.0   4.9      75.7
Mesenchymal          0.0  1.1   5.1  0.0  24.1      69.7

## Each clusterU: % of its cells at each vertex or interior

group       V1   V2    V3   V4    V5  interior
clusterU                                      
K0         0.0  1.2  23.6  0.0   0.0      75.2
K1        51.0  0.1   0.0  4.4   0.0      44.5
K2         0.0  0.5   0.0  0.0  44.8      54.7
K3         0.0  0.2   0.3  3.2   0.6      95.6
K4         1.7  0.0  19.7  0.0   0.0      78.6

## Each sample: % of its cells at each vertex or interior

group          V1   V2    V3   V4    V5  interior
sample                                           
HuMEC_Mock   12.8  0.0  20.4  0.0   0.0      66.8
HuMEC_TGFB    6.9  1.1  16.7  0.0   0.0      75.3
MCF10A_Mock   0.0  0.0   0.0  5.3  17.9      76.7
MCF10A_TGFB   0.0  1.0   0.0  0.4  32.3      66.4

## Each position: % of its cells at each vertex or interior

group      V1   V2    V3   V4    V5  interior
position                                     
inner     9.4  0.3  16.5  1.8   3.3      68.7
outer     4.7  0.7   9.8  0.1  11.2      73.4

## Each arm: % of its cells at each vertex or interior

group   V1   V2    V3   V4   V5  interior
arm                                      
Mock   8.9  0.0  14.1  1.7  5.5      69.9
TGFB   4.8  1.1  11.8  0.1  9.5      72.7

## Tumour check per vertex

        n_specialists   top_tumour  top_share_pct  n_tumours_ge5pct
vertex                                                             
V1               2130   HuMEC_Mock           67.3                 2
V2                155   HuMEC_TGFB           72.3                 2
V3               3986   HuMEC_Mock           57.3                 2
V4                285  MCF10A_Mock           93.7                 2
V5               2260  MCF10A_TGFB           60.4                 2

## Mean scores per vertex (specialists) and interior

          KS_EMT  purity
group                   
V1        -0.573   0.543
V2         0.287   0.439
V3        -0.009   0.479
V4        -0.791   0.502
V5         0.537   0.510
interior   0.010   0.145

## Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)

                   label_type       subtype  k_in_bin  fold_enrichment  q_value
archetype                                                                      
V1                   clusterU            K1      2637           6.7546      0.0
V1         cluster_annotation    Epithelial      2857           3.7356      0.0
V1                     sample    HuMEC_Mock      2141           1.9096      0.0
V1                      phase            G1      2918           1.5247      0.0
V2         cluster_annotation   Mesenchymal      1743           2.2790      0.0
V2                   clusterU            K0      2180           2.1664      0.0
V2                     sample   MCF10A_TGFB       829           1.9598      0.0
V2                        arm          TGFB      2782           1.9337      0.0
V2                     sample    HuMEC_TGFB      1953           1.9228      0.0
V3                   clusterU            K0      2298           2.2836      0.0
V3         cluster_annotation   Late-Hybrid      1291           1.6882      0.0
V3                     sample    HuMEC_TGFB      1629           1.6038      0.0
V3         cluster_annotation  Early-Hybrid      1174           1.5352      0.0
V4                     sample   MCF10A_Mock      2078           4.1635      0.0
V4                   clusterU            K3      1442           4.1235      0.0
V4                   clusterU            K1      1405           3.5989      0.0
V4         cluster_annotation    Epithelial      2228           2.9132      0.0
V4                      assay        MCF10A      2356           2.5550      0.0
V4                        arm          Mock      2643           1.6312      0.0
V5                   clusterU            K2      3056           6.1120      0.0
V5                     sample   MCF10A_TGFB      2070           4.8936      0.0
V5         cluster_annotation   Mesenchymal      2565           3.3538      0.0
V5                      assay        MCF10A      3058           3.3163      0.0
V5                      phase           G2M       728           2.1268      0.0
V5                     sample   MCF10A_Mock       988           1.9796      0.0

## Top 15 genes higher near each vertex (log-expression difference)

      V1_up     V2_up    V3_up    V4_up   V5_up
0    SPRR1A       FN1    KRT15    KRT81   KRT81
1    SPRR1B     LAMC2   CALML3      PI3    MT2A
2     KRT6B     LAMA3    KRT14     SBSN    MT1E
3     LYPD3     LAMB3     GLUL     SLPI   TGFBI
4     KRT16       VIM   S100A9    S100P  DUSP23
5     FABP5      PLAU   PYCARD   SPRR1A   CXCL3
6      SBSN      MT2A    DDIT4     CSTA   CXCL1
7      SLPI  SERPINE2     KRT5   SPRR1B    MT1X
8      CES2      KRT7     CBR1     PSCA   CXCL2
9     KRT6C    TUBA1A   S100A8   KRTDAP    DKK1
10  TACSTD2   ANGPTL4   PDLIM1    FABP5     FN1
11  SULT2B1     ITGB4      TKT     SAA1  NFKBIA
12   CRABP2     TGFBI   FGFBP1   SPINK6   DUSP1
13   S100A8    FERMT1  AKR1B10  TACSTD2    NNMT
14     CSTA      EMP3    FXYD3   CRABP2    EMP3
