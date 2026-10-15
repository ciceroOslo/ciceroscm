"""
Benchmark SSP245 emissions-driven run.

Runs the model N_RUNS times and reports per-phase and total wall-clock time.
Used to track performance before/after optimisation work.
"""

import os
import sys
import time
import statistics

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../", "src"))

from ciceroscm import CICEROSCM

DATA_DIR = os.path.join(os.path.dirname(__file__), "../", "tests", "test-data")
OUTDIR = os.path.join(os.getcwd(), "output_benchmark")
os.makedirs(OUTDIR, exist_ok=True)
N_RUNS = 5

CFG = {
    "gaspam_file": os.path.join(DATA_DIR, "gases_v1RCMIP.txt"),
    "nyend": 2100,
    "concentrations_file": os.path.join(DATA_DIR, "ssp245_conc_RCMIP.txt"),
    "emissions_file": os.path.join(DATA_DIR, "ssp245_em_RCMIP.txt"),
    "nat_ch4_file": os.path.join(DATA_DIR, "natemis_ch4.txt"),
    "nat_n2o_file": os.path.join(DATA_DIR, "natemis_n2o.txt"),
}

init_times = []
run_times = []
total_times = []

print(f"SSP245 emissions-driven benchmark  ({N_RUNS} runs, 1750-2100)\n")

for i in range(N_RUNS):
    t0 = time.perf_counter()
    cscm = CICEROSCM(CFG)
    t1 = time.perf_counter()
    cscm._run({"output_folder": OUTDIR}, make_plot=False)
    t2 = time.perf_counter()

    init_times.append(t1 - t0)
    run_times.append(t2 - t1)
    total_times.append(t2 - t0)
    print(f"  Run {i + 1}: init={t1-t0:.3f}s  _run={t2-t1:.3f}s  total={t2-t0:.3f}s")

print()
print(f"  init  mean={statistics.mean(init_times):.3f}s  std={statistics.stdev(init_times):.3f}s")
print(f"  _run  mean={statistics.mean(run_times):.3f}s  std={statistics.stdev(run_times):.3f}s")
print(f"  total mean={statistics.mean(total_times):.3f}s  std={statistics.stdev(total_times):.3f}s")
