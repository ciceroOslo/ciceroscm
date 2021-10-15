import pytest
import os
import numpy as np
import pandas as pd
import pandas.testing as pdt

from ciceroscm import CICEROSCM


def test_ciceroscm():
    cscm = CICEROSCM()


def test_cfg(test_data_dir):
    cscm = CICEROSCM()
    cscm._run(
        {"gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),},
        {"forc_file": os.path.join(test_data_dir, "CO2_1pros.txt")},
    )
    cscm._run(
        {"gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"), "sunvolc": 1},
        {"forc_file": os.path.join(test_data_dir, "CO2_1pros.txt")},
    )
