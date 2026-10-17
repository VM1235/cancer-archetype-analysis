#!/usr/bin/env python3
"""Restartable job queue. From Project_1/:  python3 docs/cleaning_audit/run_queue.py   (NPAR=4 fits at a time by default)
Reads jobs.tsv in order; a job whose result json already exists is skipped, so rerunning after a stop is safe."""
import os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
subprocess.run([sys.executable, "docs/cleaning_audit/gen_jobs.py"], check=True, stdout=subprocess.DEVNULL)
os.makedirs("docs/cleaning_audit/logs", exist_ok=True)
jobs = [l.rstrip("\n").split("\t") for l in open("docs/cleaning_audit/jobs.tsv")]
def log(msg):
    with open("docs/cleaning_audit/logs/queue.log", "a") as f: f.write(msg + "\n")
log("QUEUE START"); 
with ThreadPoolExecutor(int(os.environ.get("NPAR", 4))) as ex:
    list(ex.map(lambda j: subprocess.run(["bash", "docs/cleaning_audit/run_job.sh", j[0], j[1], j[2]]), jobs))
log("QUEUE END")
