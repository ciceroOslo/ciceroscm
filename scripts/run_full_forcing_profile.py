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
#os.getcwd() gets the path of where you are running from
outdir = os.path.join(os.getcwd(), "output_test")

# Starting profiler before starting run
pr = cProfile.Profile()
pr.enable()

cscm = CICEROSCM(
    {
        "gaspamfile": os.path.join(data_dir, "gases_v1RCMIP.txt"),
        "sunvolc": 1,
        "nyend": 2100,
        "forc_file": os.path.join(data_dir, "CO2_1pros.txt"),
    },
)

cscm._run({"output_folder": outdir})

# Stop counting with profiler
pr.disable()

# Sorting and printing profiler results
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr,stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
