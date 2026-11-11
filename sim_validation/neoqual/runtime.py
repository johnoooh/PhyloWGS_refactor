#!/usr/bin/env python3
"""
runtime.py — Python vs Go runtime for one NeoQual sample -> <work>/out/runtime_<k>.json.

Python: each chain's mcmc_samples.txt has header Iteration/LLH/Time, Time being
the seconds that iteration took; a chain's seconds = sum over its 3500 rows.
The 15 chains ran as separate single-threaded processes, so CPU = sum of chain
seconds, and wall = max chain seconds (assumes the 15 ran concurrently).

Go: wall (Elapsed), CPU (TotalCPU) and cores (AllocCPUS) from sacct for the
Slurm task that produced the chains: <work>/go_time/<k>.job (written by task.sh
for new runs) or the manifest's go_job (reused runs). Go's own per-chain wall
times (chains/summary.json chain_times, contended for 12 cores) are kept too.
A missing or unusable record becomes go_status != "ok", never a number.
"""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

PY_ROWS = 3500
GO_DUR = re.compile(r"(\d+(?:\.\d+)?)(h|ms|us|µs|ns|m|s)")
UNIT = {"h": 3600, "m": 60, "s": 1, "ms": 1e-3, "us": 1e-6, "µs": 1e-6, "ns": 1e-9}


def py_chain_seconds(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        if head != ["Iteration", "LLH", "Time"]:
            raise ValueError(f"{path}: unexpected header {head}")
        t = [float(ln.split("\t")[2]) for ln in f if ln.strip()]
    if len(t) != PY_ROWS:
        raise ValueError(f"{path}: {len(t)} rows, expected {PY_ROWS}")
    return sum(t)


def go_duration(s):
    """Go time.Duration string ('1h2m3.5s', '850ms') -> seconds."""
    parts = GO_DUR.findall(s)
    if not parts or "".join(a + b for a, b in parts) != s:
        raise ValueError(f"bad Go duration {s!r}")
    return sum(float(a) * UNIT[b] for a, b in parts)


def slurm_seconds(s):
    """sacct '[D-][HH:]MM:SS[.mmm]' -> seconds."""
    d, _, hms = s.rpartition("-")
    parts = [float(x) for x in hms.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    return (int(d) if d else 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def sacct_job(job):
    out = subprocess.run(["sacct", "-n", "-P", "-j", job, "-o", "JobID,State,Elapsed,TotalCPU,AllocCPUS"],
                         capture_output=True, text=True, check=True).stdout
    rows = [r.split("|") for r in out.splitlines()]
    top = [r for r in rows if r[0] == job]
    if not top:
        return {"go_status": "sacct_no_record"}
    _, state, el, cpu, n = top[0]
    rec = {"go_job_state": state, "go_wall_s": slurm_seconds(el), "go_cpu_s": slurm_seconds(cpu),
           "go_cores": int(n)}
    if state != "COMPLETED":
        rec["go_status"] = f"job_{state.split()[0].lower()}"
    elif rec["go_cpu_s"] <= 0:
        rec["go_status"] = "cpu_unrecorded"
    else:
        rec["go_status"] = "ok"
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("work")
    ap.add_argument("k", type=int)
    args = ap.parse_args()
    W = Path(args.work)
    e = json.load(open(W / "manifest.json"))[args.k]

    py = [py_chain_seconds(Path(p).parent / "mcmc_samples.txt") for p in e["py_src"]]
    rec = {"k": e["k"], "cls": e["cls"], "n_ssm": e["n_ssm"], "go_src": e["go_src"],
           "py_chain_s": py, "py_cpu_s": sum(py), "py_wall_s": max(py)}

    jf = W / f"go_time/{e['k']}.job"
    job = jf.read_text().strip() if jf.exists() else e.get("go_job")
    rec.update(sacct_job(job) if job else {"go_status": "no_job_record"})
    rec["go_job"] = job
    summ = Path(e["go_dir"]) / "chains/summary.json"
    if summ.exists():
        rec["go_chain_s"] = [go_duration(s) for s in json.load(open(summ))["chain_times"]]
    if rec["go_status"] == "ok":
        rec["wall_ratio_go_py"] = rec["go_wall_s"] / rec["py_wall_s"]
        rec["cpu_ratio_go_py"] = rec["go_cpu_s"] / rec["py_cpu_s"]
    (W / "out").mkdir(exist_ok=True)
    json.dump(rec, open(W / f"out/runtime_{e['k']}.json", "w"), indent=1)
    print({x: rec.get(x) for x in ("go_status", "py_wall_s", "go_wall_s", "wall_ratio_go_py", "cpu_ratio_go_py")},
          file=sys.stderr)


if __name__ == "__main__":
    main()
