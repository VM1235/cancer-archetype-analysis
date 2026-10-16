# Groves closest-5% specialists: all/pooled, k=3

Reference cells: 2200; specialists 15.0%

## Fit cells: % of each time in each class (rows sum to 100)

groves5   V1   V2    V3  generalist    n
time                                    
0d       8.6  1.6   0.5        89.3  429
8h       4.6  6.3   2.5        86.7  607
1d       5.5  7.4  11.5        75.5  564
3d       2.9  2.5   5.7        88.9  477
7d       0.0  8.9   0.8        90.2  123

## Removal-arm guests: % of each time in each class (rows sum to 100)

groves5   V1    V2   V3  generalist    n
time                                    
8h_rm    0.6   7.6  0.4        91.4  709
1d_rm    1.0  18.2  0.1        80.7  726
3d_rm    0.0   5.4  1.0        93.6  297

## Fit cells: % of each cell_line in each class (rows sum to 100)

groves5      V1    V2    V3  generalist    n
cell_line                                   
A549        0.0   0.0   0.0       100.0  550
DU145       0.0  20.0   0.0        80.0  550
MCF7       20.0   0.0   0.0        80.0  550
OVCA420     0.0   0.0  20.0        80.0  550

## Fit cells: % of each phase in each class (rows sum to 100)

groves5   V1    V2   V3  generalist     n
phase                                    
G1       7.5   3.0  8.1        81.4  1130
S        3.8   3.6  1.4        91.2   658
G2M      0.0  12.6  2.4        85.0   412

## Fit cells: % of each mix in each class (rows sum to 100)

groves5    V1    V2   V3  generalist    n
mix                                      
Mix1     14.8   0.8  2.9        81.5  480
Mix2      8.2   0.5  5.6        85.7  426
Mix3a     0.4   8.4  8.1        83.2  273
Mix3b     0.0  15.9  5.6        78.4  232
Mix4a     0.0   4.3  3.8        91.9  397
Mix4b     0.8   6.9  5.6        86.7  392

## Rule comparison (counts, fit cells)

            groves5  weight095  gap035
V1              110         33     549
V2              110         15     553
V3              110          7     533
generalist     1870       2145     565

## % of each time point called specialist, by rule

      groves5  weight095  gap035
time                            
0d       10.7        3.0    72.5
8h       13.3        3.1    73.3
1d       24.5        2.0    77.3
3d       11.1        2.1    72.7
7d        9.8        1.6    78.0
