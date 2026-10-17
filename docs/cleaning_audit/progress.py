#!/usr/bin/env python3
"""Write PROGRESS.md: which queued fits are finished. From Project_1/: python3 docs/cleaning_audit/progress.py"""
import os, re, collections, time
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
jobs = [l.rstrip("\n").split("\t") for l in open("docs/cleaning_audit/jobs.tsv")]
grp = lambda i: re.sub(r"_(k\d.*|t\d+)$", "", i.split("_")[0] + "_" + ("sweep" if re.search(r"_t\d+$", i) else "main")) if False else (
    "Cook (human tech genes)" if i.startswith("cook") else "GSE114687 UMI sweep" if re.search(r"_t\d+$", i) else
    "GSE114687 (human tech genes)" if i.startswith("gse") else "4T1 contamination cutoff 0.2/0.5" if i.startswith("t4") else "PDAC cleaned KS fits")
tot = collections.Counter(); done = collections.Counter(); todo = []
for i, js, _ in jobs:
    g = grp(i); tot[g] += 1
    if os.path.exists(js): done[g] += 1
    else: todo.append(i)
fails = [l for l in open("docs/cleaning_audit/logs/queue.log") if " FAIL " in l] if os.path.exists("docs/cleaning_audit/logs/queue.log") else []
lines = ["# Cleaning-audit progress", f"_Updated {time.strftime('%Y-%m-%d %H:%M')}_", "",
         f"**{sum(done.values())} of {sum(tot.values())} fits finished.**", "", "| Group | Done | Total |", "|---|---|---|"]
lines += [f"| {g} | {done[g]} | {tot[g]} |" for g in tot]
lines += ["", f"Failed jobs: {len(fails)}"] + [f"- {l.strip()}" for l in fails[:20]]
lines += ["", "## To resume (safe to rerun; finished fits are skipped)", "```", "cd /Users/apple/Desktop/thesis_iisc/Project_1", "nohup python3 docs/cleaning_audit/run_queue.py > docs/cleaning_audit/logs/queue_stdout.log 2>&1 &", "```",
          "Check it is running: `ps aux | grep 08_within_unit | grep -v grep | wc -l` (and `tail docs/cleaning_audit/logs/queue.log`).", "",
          "Then: `python3 docs/cleaning_audit/compare.py` writes COMPARE.md (old vs new p-values)."]
open("docs/cleaning_audit/PROGRESS.md", "w").write("\n".join(lines) + "\n")
print("\n".join(lines[:14]))
