import os
import shutil
import sys
import cProfile, pstats, io
import pandas as pd

# Adding location of source code to system path
# os.path.dirname(__file__) gives the directory of
# current file. Put in updated path if running script from elsewhere
# os.path joins all the folders of a path together in a
# system independent way (i.e. will work equally well on Windows, linux etc)
sys.path.append(os.path.join(os.path.dirname(__file__), "../", "src"))

from ciceroscm import CICEROSCM

data_dir = os.path.join(os.path.dirname(__file__), "../", "tests", "test-data")
# os.getcwd() gets the path of where you are running from
outdir = os.path.join(os.getcwd(), "output_test")

# Starting profiler before starting run
pr = cProfile.Profile()
pr.enable()

cscm = CICEROSCM(
    {
        "gaspam_file": os.path.join(data_dir, "gases_v1RCMIP.txt"),
        "nyend": 2100,
        "concentrations_file": os.path.join(data_dir, "ssp245_conc_RCMIP.txt"),
        "emissions_file": os.path.join(data_dir, "ssp245_em_RCMIP.txt"),
        "nat_ch4_file": os.path.join(data_dir, "natemis_ch4.txt"),
        "nat_n2o_file": os.path.join(data_dir, "natemis_n2o.txt"),
        "carbon_cycle_model": "box",
        "thermal_model": "twolayer",
    }
)


cscm._run({"output_folder": outdir}, make_plot=False)

# Stop counting with profiler
pr.disable()

# Sorting and printing profiler results
s = io.StringIO()
sortby = "cumulative"
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
