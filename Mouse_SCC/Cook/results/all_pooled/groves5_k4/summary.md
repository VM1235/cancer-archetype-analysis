# Groves closest-5% specialists: all/pooled, k=4

Reference cells: 2200; specialists 20.0%

## Fit cells: % of each time in each class (rows sum to 100)

groves5    V1   V2   V3   V4  generalist    n
time                                         
0d        0.2  2.6  8.4  1.2        87.6  429
8h        2.8  5.9  4.6  5.6        81.1  607
1d       10.8  7.6  5.7  6.0        69.9  564
3d        6.3  2.3  2.9  7.1        81.3  477
7d        0.8  7.3  0.0  2.4        89.4  123

## Removal-arm guests: % of each time in each class (rows sum to 100)

groves5   V1    V2   V3   V4  generalist    n
time                                         
8h_rm    0.7   5.6  0.6  3.1        90.0  709
1d_rm    0.1  15.6  1.4  0.6        82.4  726
3d_rm    1.3   5.7  0.3  0.3        92.3  297

## Fit cells: % of each cell_line in each class (rows sum to 100)

groves5      V1    V2    V3    V4  generalist    n
cell_line                                         
A549        0.0   0.0   0.0  20.0        80.0  550
DU145       0.0  20.0   0.0   0.0        80.0  550
MCF7        0.0   0.0  20.0   0.0        80.0  550
OVCA420    20.0   0.0   0.0   0.0        80.0  550

## Fit cells: % of each phase in each class (rows sum to 100)

groves5   V1    V2   V3   V4  generalist     n
phase                                         
G1       7.5   2.9  7.7  5.9        75.9  1130
S        1.7   4.1  3.5  5.0        85.7   658
G2M      3.4  12.1  0.0  2.4        82.0   412

## Fit cells: % of each mix in each class (rows sum to 100)

groves5   V1    V2    V3   V4  generalist    n
mix                                           
Mix1     3.1   1.2  14.4  1.5        79.8  480
Mix2     5.9   0.5   8.5  4.5        80.8  426
Mix3a    7.0   7.7   0.7  5.9        78.8  273
Mix3b    4.3  16.4   0.0  1.7        77.6  232
Mix4a    4.8   3.3   0.0  6.3        85.6  397
Mix4b    5.6   7.7   0.8  9.9        76.0  392

## Rule comparison (counts, fit cells)

            groves5  weight095  gap035
V1              110          1     507
V2              110          1     544
V3              110         15     550
V4              110          5     549
generalist     1760       2178      50

## % of each time point called specialist, by rule

      groves5  weight095  gap035
time                            
0d       12.4        1.4    96.5
8h       18.9        1.0    98.0
1d       30.1        0.7    97.7
3d       18.7        1.3    98.5
7d       10.6        0.0    97.6
