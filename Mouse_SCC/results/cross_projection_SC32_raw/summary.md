# Cross-projection into the SC32 triangle (k=3, uncorrected)

Host corners (from 08_within_tumour.py --compare): V1 (Early-Hybrid 50%); V2 (Late-EMT 56%); V3 (Epithelial 79%)

Rows: tumour x state. Columns: % of those cells at each host corner (purity gap >= 0.35),
% interior, and % inside the triangle (all weights >= 0).

     tumour        state    n  V1 (Early-Hybrid 50%)  V2 (Late-EMT 56%)  V3 (Epithelial 79%)  interior  inside
SC32 (host)   Epithelial  589                   31.1                0.0                 15.3      53.7    54.3
SC32 (host) Early-Hybrid  956                   19.2                9.1                  2.5      69.1    75.4
SC32 (host)       Hybrid  139                    0.0               56.1                  0.0      43.9    27.3
SC32 (host)     Late-EMT  226                    0.0               91.2                  0.0       8.8    36.7
       SC15   Epithelial  416                    2.6                0.0                 71.6      25.7    45.4
       SC15 Early-Hybrid  818                    2.8               15.5                  6.7      74.9    62.7
       SC15       Hybrid   98                    0.0               52.0                  0.0      48.0    40.8
       SC15     Late-EMT  244                    0.0               71.3                  0.0      28.7    36.5
       SC46   Epithelial   44                    0.0                0.0                  9.1      90.9    68.2
       SC46 Early-Hybrid  430                    0.0               28.1                  0.5      71.4    74.4
       SC46       Hybrid  330                    0.0               29.7                  0.0      70.3    23.0
       SC46     Late-EMT  715                    0.0               93.6                  0.0       6.4    30.5
       SC33   Epithelial  232                   66.8                0.0                  0.4      32.8    64.7
       SC33 Early-Hybrid 1049                   11.3                3.7                  0.2      84.7    72.9
       SC33       Hybrid   43                    4.7               51.2                  0.0      44.2    18.6
       SC33     Late-EMT   55                    0.0               96.4                  0.0       3.6    60.0
      SC2b1   Epithelial   23                   82.6                0.0                  0.0      17.4   100.0
      SC2b1       Hybrid  326                    0.0               49.1                  0.0      50.9     4.9
      SC2b1     Late-EMT  739                    0.0              100.0                  0.0       0.0    36.9
      SC2b3   Epithelial  187                   28.3                0.0                  9.1      62.6    93.0
      SC2b3 Early-Hybrid  155                    7.1               11.0                  6.5      75.5    91.0
      SC2b3       Hybrid  965                    0.0               49.5                  0.0      50.5    44.7
      SC2b3     Late-EMT 1654                    0.0               93.5                  0.0       6.5    29.0
       SC35 Early-Hybrid   46                    2.2               60.9                  0.0      37.0    37.0
       SC35       Hybrid  981                    0.0               45.5                  0.0      54.5    14.9
       SC35     Late-EMT 1545                    0.0               98.6                  0.0       1.4    25.2
       SC43       Hybrid  546                    0.0               43.2                  0.0      56.8    12.1
       SC43     Late-EMT  478                    0.0               95.4                  0.0       4.6    51.3
        MU3   Epithelial  160                    9.4                0.0                  5.0      85.6    98.8
        MU3 Early-Hybrid  654                    0.0                4.3                  0.6      95.1    89.6
        MU3     Late-EMT   25                    0.0               72.0                  0.0      28.0     0.0
        MU4   Epithelial  393                    2.0                0.0                 28.5      69.5    91.3
        MU4 Early-Hybrid  457                    0.7                0.4                 17.1      81.8    93.2
        MU4       Hybrid   25                    0.0               80.0                  0.0      20.0    76.0
        MU4     Late-EMT   64                    0.0               82.8                  0.0      17.2    10.9
        MU7   Epithelial  782                    0.5                0.0                 70.5      29.0    79.2
        MU7 Early-Hybrid 1774                    0.0                0.1                 14.1      85.8    54.5
        MU7       Hybrid   65                    0.0               12.3                  0.0      87.7    95.4
        MU7     Late-EMT  180                    0.0               18.9                  0.0      81.1    21.1
        MU8   Epithelial 1307                    0.8                0.0                 45.3      53.9    89.1
        MU8 Early-Hybrid 3629                    0.2                0.2                  4.1      95.5    86.1
        MU8       Hybrid   77                    0.0               41.6                  0.0      58.4    84.4
        MU8     Late-EMT  295                    0.0               28.8                  0.0      71.2    10.5

## All scRNA guest cells pooled, by state (%)

group           V1    V2    V3  interior  inside
state                                           
Epithelial    26.7   0.0  35.1      38.3    63.1
Early-Hybrid   6.2  13.5   2.7      77.5    70.1
Hybrid         0.1  45.3   0.0      54.6    23.8
Late-EMT       0.0  95.0   0.0       5.0    31.8
