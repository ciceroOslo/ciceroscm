import os
import shutil
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../", "src"))

print(sys.path)
from ciceroscm import CICEROSCM

data_dir = os.path.join(os.path.dirname(__file__), "../", "tests", "test-data")

cscm = CICEROSCM(
    {
        "gaspamfile": os.path.join(data_dir, "gases_v1RCMIP.txt"),
        "nystart": 1900,
        "emstart": 1950,
        "nyend": 2050,
        "concentrations_file": os.path.join(data_dir, "ssp245_conc_RCMIP.txt"),
        "emissions_file": os.path.join(data_dir, "ssp245_em_RCMIP.txt"),
        "nat_ch4_file": os.path.join(data_dir, "natemis_ch4.txt"),
        "nat_n2o_file": os.path.join(data_dir, "natemis_n2o.txt"),
    },
)

outdir = os.path.join(os.getcwd(), "output_test")

cscm._run({"output_folder": outdir})
