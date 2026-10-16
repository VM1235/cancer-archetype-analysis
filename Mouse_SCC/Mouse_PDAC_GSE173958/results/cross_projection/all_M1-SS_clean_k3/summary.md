# Cross-projection into M1-SS (genes all_clean, k=3)

## Host corners (closest 10% of the host's cells)

- V1: top genes Lgals1, Lgals2, S100a6, Krt8, Tmsb4x, Ccnd1; cycling 45%; KS_EMT +0.29; UMI 37272
- V2: top genes Atp1b1, Ttr, Fxyd3, C3, Crisp1, Dmbt1; cycling 2%; KS_EMT -0.56; UMI 38082
- V3: top genes Gm26870, Neat1, Id3, Fosb, Trf, Inhba; cycling 2%; KS_EMT +0.54; UMI 26321

## % of each library's cells at each host corner (purity gap >= 0.35), interior, and inside

     library     n   V1   V2   V3  interior  inside
       M1-PT  4566  6.6 21.2  8.0      64.2    78.8
M1-SS (host)  2803 11.6 15.9 16.5      56.0    67.4
      M1-Met  5761 10.8  5.7  7.2      76.2    79.3
    M1-Liver  3938  9.6 12.2  5.1      73.1    90.4
     M1-Lung   216 19.9 14.4  8.3      57.4    79.6
      M1-CTC   660 16.7  3.2  4.4      75.8    90.9
       M2-PT  5037  6.1 22.2  0.4      71.3    81.2
      M2-Met 14978  8.6  9.8  0.0      81.6    84.6
    M2-Liver  6061  6.9  9.0  0.1      83.9    91.7
     M2-Lung   246 15.0  5.3  0.4      79.3    71.1
    iv-5cell  7126  1.5  0.1  0.0      98.4    53.6
  iv-100cell  6378  3.1  0.0  0.0      96.9    31.4

## All other in-vivo cells, by EMT state: % at each host corner / interior

group           V1    V2    V3  interior
state                                   
Epithelial     1.7  27.8   0.1      70.5
Early-Hybrid  10.9   1.7   1.6      85.7
Late-Hybrid   18.5   0.2   7.3      74.1
Mesenchymal   16.1   0.1  13.6      70.3
