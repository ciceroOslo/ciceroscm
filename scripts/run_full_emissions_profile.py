import os
import shutil
import sys
import cProfile, pstats, io
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../", "src"))

from ciceroscm import CICEROSCM

data_dir = os.path.join(os.path.dirname(__file__), "../", "tests", "test-data")
outdir = os.path.join(os.getcwd(), "output_test") 
pr = cProfile.Profile()
pr.enable()
cscm = CICEROSCM(
    {
        "gaspamfile": os.path.join(data_dir, "gases_v1RCMIP.txt"),
        "nyend": 2100,
        "concentrations_file": os.path.join(data_dir, "ssp245_conc_RCMIP.txt"),
        "emissions_file": os.path.join(data_dir, "ssp245_em_RCMIP.txt"),
        "nat_ch4_file": os.path.join(data_dir, "natemis_ch4.txt"),
        "nat_n2o_file": os.path.join(data_dir, "natemis_n2o.txt"),
    },
)


cscm._run({"output_folder": outdir})

pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr,stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
