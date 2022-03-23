import os
import shutil
import sys
import cProfile, pstats, io
import pandas as pd
#from pstats import SortKey

sys.path.append(os.path.join(os.path.dirname(__file__), "../", "src"))

from ciceroscm import CICEROSCM

data_dir = os.path.join(os.path.dirname(__file__), "../", "tests", "test-data")
outdir = os.path.join(os.getcwd(), "output_test")
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

pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr,stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
