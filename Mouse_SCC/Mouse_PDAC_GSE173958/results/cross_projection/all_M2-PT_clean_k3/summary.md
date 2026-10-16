# Cross-projection into M2-PT (genes all_clean, k=3)

## Host corners (closest 10% of the host's cells)

- V1: top genes H2afz, Hmgb2, Birc5, Cenpa, Cycs, Ran; cycling 84%; KS_EMT -0.01; UMI 38120
- V2: top genes Ttr, 2210407C18Rik, Dmbt1, Sptssb, Lgals2, Gsta4; cycling 0%; KS_EMT -0.54; UMI 32304
- V3: top genes Ero1l, Vegfa, Hilpda, Malat1, Bnip3, Slc2a1; cycling 1%; KS_EMT +0.20; UMI 22074

## % of each library's cells at each host corner (purity gap >= 0.35), interior, and inside

     library     n   V1   V2   V3  interior  inside
       M1-PT  4566  4.7 14.5 14.7      66.0    80.1
       M1-SS  2803  5.1  9.7 19.8      65.4    90.8
      M1-Met  5761  9.7  6.0  9.9      74.4    89.2
    M1-Liver  3938 10.4 11.1  2.9      75.6    90.9
     M1-Lung   216 18.1  8.3  3.7      69.9    83.8
      M1-CTC   660 11.7  1.1  9.2      78.0    86.1
M2-PT (host)  5037 11.3 16.7  7.3      64.8    85.6
      M2-Met 14978 30.4  3.9  2.5      63.2    91.6
    M2-Liver  6061 39.1  5.6  1.2      54.1    90.6
     M2-Lung   246 51.6  2.0  2.0      44.3    81.3
    iv-5cell  7126 84.5  0.0  0.0      15.5    96.8
  iv-100cell  6378 93.4  0.0  0.0       6.6    97.4

## All other in-vivo cells, by EMT state: % at each host corner / interior

group           V1    V2    V3  interior
state                                   
Epithelial    12.4  17.3   1.8      68.4
Early-Hybrid  33.0   0.8   4.4      61.8
Late-Hybrid   21.0   0.0  13.5      65.4
Mesenchymal    9.8   0.0  23.3      66.9
