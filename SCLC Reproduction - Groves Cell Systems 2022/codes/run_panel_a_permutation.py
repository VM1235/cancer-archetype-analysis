#!/usr/bin/env python3
"""Scale Panel A t-ratio tests to 1000 shuffles using saved PCA/archetypes."""

from pathlib import Path
import sys

SCLC = Path(__file__).resolve().parents[1]
ROOT = SCLC.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from src.archetypes import permutation_t_ratio

OUT = SCLC / "results" / "panel_a"
PAPER = {
    4: {"p": 0.059, "t_ratio": 0.247},
    5: {"p": 0.034, "t_ratio": 0.107},
    6: {"p": 0.016, "t_ratio": 0.043},
}
K_PERM = (4, 5, 6)
N_PERM = 1000
DELTA = 0.1
SEED = 0


def main():
    scores12 = np.load(OUT / "pc_scores_12.npy")
    rows = []
    for k in K_PERM:
        archetypes = np.load(OUT / f"archetypes_k{k}.npy")
        print(f"t-ratio permutation test: k={k}, n_perm={N_PERM} ...", flush=True)
        result = permutation_t_ratio(
            scores12,
            k=k,
            n_perm=N_PERM,
            delta=DELTA,
            seed=SEED + k,
            observed_archetypes=archetypes,
            verbose=True,
        )
        np.save(OUT / f"null_t_ratios_k{k}_n{N_PERM}.npy", result["null_t_ratios"])
        paper = PAPER[k]
        rows.append(
            {
                "k": k,
                "t_ratio": result["t_ratio"],
                "p_value": result["p_value"],
                "n_perm": result["n_perm"],
                "n_success": result["n_success"],
                "n_fail": result["n_fail"],
                "paper_p": paper["p"],
                "paper_t_ratio": paper["t_ratio"],
            }
        )
        print(
            f"  k={k}: t={result['t_ratio']:.4f} (paper {paper['t_ratio']:.3f}), "
            f"p={result['p_value']:.4f} from {result['n_success']}/{N_PERM} "
            f"(paper p={paper['p']:.3f})",
            flush=True,
        )
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "t_ratio_1000.csv", index=False)
    print("\nFinal 1000-shuffle table:")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
