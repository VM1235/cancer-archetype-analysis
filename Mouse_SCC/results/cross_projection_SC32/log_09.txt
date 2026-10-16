# Cross-projection into the SC32 triangle (k=3, V1-corrected)

Host corners (from 08_within_tumour.py --compare): V1 (Late-EMT 56%); V2 (Epithelial 97%); V3 (Early-Hybrid 99%)

Rows: tumour x state. Columns: % of those cells at each host corner (purity gap >= 0.35),
% interior, and % inside the triangle (all weights >= 0).

     tumour        state    n  V1 (Late-EMT 56%)  V2 (Epithelial 97%)  V3 (Early-Hybrid 99%)  interior  inside
SC32 (host)   Epithelial  589                0.0                 58.1                    0.3      41.6    34.0
SC32 (host) Early-Hybrid  956                8.7                  0.7                   42.4      48.2    63.9
SC32 (host)       Hybrid  139               32.4                  1.4                    0.7      65.5    97.1
SC32 (host)     Late-EMT  226               73.5                  0.0                    0.0      26.5    50.0
       SC15   Epithelial  416                0.0                 15.9                    0.5      83.7    32.0
       SC15 Early-Hybrid  818               27.6                  0.0                   14.1      58.3    80.8
       SC15       Hybrid   98               22.4                  0.0                    0.0      77.6    96.9
       SC15     Late-EMT  244               73.8                  0.0                    0.4      25.8    63.9
       SC46   Epithelial   44                9.1                  0.0                    2.3      88.6    95.5
       SC46 Early-Hybrid  430               21.6                  0.0                    3.7      74.7    84.7
       SC46       Hybrid  330                1.2                  0.0                    2.4      96.4    97.3
       SC46     Late-EMT  715               79.0                  0.0                    0.0      21.0    49.1
       SC33   Epithelial  232                0.4                  2.2                   11.6      85.8    51.7
       SC33 Early-Hybrid 1049                1.2                  0.0                   48.0      50.8    84.9
       SC33       Hybrid   43               16.3                  0.0                    4.7      79.1    93.0
       SC33     Late-EMT   55               67.3                  0.0                    0.0      32.7    38.2
      SC2b1   Epithelial   23                0.0                 34.8                    0.0      65.2    82.6
      SC2b1       Hybrid  326               16.0                  0.0                    0.0      84.0    99.4
      SC2b1     Late-EMT  739               98.4                  0.0                    0.0       1.6    67.8
      SC2b3   Epithelial  187                0.0                 46.5                    0.0      53.5    77.5
      SC2b3 Early-Hybrid  155               11.0                  8.4                    0.6      80.0    92.3
      SC2b3       Hybrid  965               25.1                  0.2                    0.0      74.7    98.4
      SC2b3     Late-EMT 1654               97.8                  0.0                    0.0       2.2    60.1
       SC35 Early-Hybrid   46               65.2                  0.0                    4.3      30.4    58.7
       SC35       Hybrid  981                6.8                  0.0                    0.0      93.2    94.4
       SC35     Late-EMT 1545               92.2                  0.0                    0.0       7.8    49.4
       SC43       Hybrid  546                5.5                  0.0                    0.0      94.5   100.0
       SC43     Late-EMT  478               86.4                  0.0                    0.0      13.6    83.9
        MU3   Epithelial  160                0.0                  0.0                   12.5      87.5    96.2
        MU3 Early-Hybrid  654                2.9                  0.0                   41.1      56.0    73.5
        MU3     Late-EMT   25               88.0                  0.0                    0.0      12.0    12.0
        MU4   Epithelial  393                0.0                 17.8                    0.5      81.7    93.4
        MU4 Early-Hybrid  457                2.8                  0.4                   10.9      85.8    94.5
        MU4       Hybrid   25                4.0                  0.0                    0.0      96.0    88.0
        MU4     Late-EMT   64               85.9                  0.0                    0.0      14.1    25.0
        MU7   Epithelial  782                0.1                  2.9                    0.3      96.7    93.0
        MU7 Early-Hybrid 1774               15.7                  0.0                    6.1      78.2    90.2
        MU7       Hybrid   65                3.1                  0.0                    0.0      96.9    98.5
        MU7     Late-EMT  180               72.8                  0.0                    0.0      27.2    70.6
        MU8   Epithelial 1307                0.0                  9.3                    0.5      90.2    94.5
        MU8 Early-Hybrid 3629                7.9                  0.0                   29.5      62.6    86.7
        MU8       Hybrid   77               36.4                  0.0                    0.0      63.6    89.6
        MU8     Late-EMT  295               91.9                  0.0                    0.0       8.1    73.9

## All scRNA guest cells pooled, by state (%)

group           V1    V2    V3  interior  inside
state                                           
Epithelial     0.5  19.0   3.3      77.2    50.6
Early-Hybrid  15.3   0.5  25.3      58.9    83.6
Hybrid        12.9   0.1   0.3      86.7    97.4
Late-EMT      91.4   0.0   0.0       8.6    58.7
