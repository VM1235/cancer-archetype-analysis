# Cleaning-rule audit: results (2026-10-17)

All 104 reruns finished with no failures. Full tables: `COMPARE.md`. Hallmark scores: `hallmark_scc_4t1.csv`. Why/how: `HANDOFF.md`.

## What was wrong
The technical-gene filter used mouse symbols, so it removed nothing in the two human datasets (GSE114687, Cook); Cook's cleaning script never applied it. The report's "technical genes dropped" was wrong for these two. The 10,000-UMI cutoff also keeps very different fractions of cells (Cook 89%, PDAC 80%, SCC 71%, 4T1 69%, GSE114687 45%).

## Results
1. **Human gene-filter fix changes nothing.** Cook: 18 of 18 fits keep their pass/fail (A549 all-genes k=3 still p=0.010 and still fails when time-balanced; MCF7 KS k=3 still p=0.030). GSE114687: 18 of 24 cleaned fits pass before and after (all genes 11 of 12, KS 7 of 12); the one visible p change (HuMEC all genes k=3, 0.925 to 0.886) stays a fail. The emailed "18 of 24" stands.
2. **4T1 contamination cutoff (0.2 / 0.3 / 0.5), within clone, k=3.** All genes: 0 of 5 pass at every cutoff. KS genes: EM2 passes at all three (0.035 / 0.015 / 0.050, borderline at 0.5); parental 4T1 passes only at 0.5 (p=0.005, probably contaminating cells; corner genes not checked). "Almost nothing within a clone" holds.
3. **PDAC cleaned KS fits (never run before).** k=3: 0 of 9 pass (all genes: 7 of 9). k=4: 3 of 9 (M2-Met, M1-SS, M2). EMT genes alone do not give a within-library triangle at k=3; the all-genes triangles come from other programmes.
4. **GSE114687 UMI cutoff 5,000 / 10,000 / 15,000 (human fix on), k=3, 6 units.**

| | 5,000 | 10,000 | 15,000 |
|---|---|---|---|
| All genes | 5 of 6 | 5 of 6 | 6 of 6 |
| KS genes | 1 of 6 | 3 of 6 | 4 of 6 |

   All genes: MCF10A Mock, MCF10A TGFB and MCF10A (both arms) pass at every cutoff; HuMEC is fragile (HuMEC TGFB fails at 5k, p=0.025 at 10k; HuMEC both arms fails at 10k only). KS genes: the KS-only triangle appears only at >= 10,000 UMI (fails at 5,000). 15,000 keeps few cells (284 to 1,838), so a pass there is a smaller sample.
5. **Hallmark check of corner names, SCC and 4T1 (new).** 4T1: the cycling corner is top on E2F / G2-M in pooled cleaned and every within-clone fit (94 to 97% dividing cells); clone M's corner scores TGF-b / hypoxia; no corner is named EMT, and the epithelial-leaning corner has no strong Hallmark set. SCC: the proliferating corners are top on E2F / G2-M, the late-EMT corner top on EMT; the "technical" corner has no set above about 0.25 (consistent with its name). Existing checks cover PDAC, Cook, GSE114687.

## What it means for the claims
- "Pooled EMT genes give no triangle beyond sample identity", "corners are programmes, not EMT states" and "within one clone almost nothing" hold under the cleaning variants tried.
- "Triangle inside a sample after cleaning" (GSE114687) holds with all genes at 5k to 15k for MCF10A; HuMEC is cutoff-sensitive.
- "KS EMT genes alone give a triangle" in GSE114687 needs cells >= 10,000 UMI.

## Not done
SCC has no cleaning arm (its issue is the assay); pooled cleaned fits exist only for 4T1 and PDAC; path test not run on SCC, PDAC, 4T1; no UMI sweep for 4T1/PDAC; a planted-triangle power check was not run. Report wording ("technical genes dropped") still needs correcting.
