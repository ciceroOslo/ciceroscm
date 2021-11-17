import os

import numpy as np
import pandas as pd

from ciceroscm  import concentrations_emissions_handler

def test_read_components(test_data_dir):
    df_gas = concentrations_emissions_handler.read_components(os.path.join(test_data_dir, "gases_v1RCMIP.txt"))
    print(df_gas.columns)
    assert df_gas.columns.tolist() == [
        "EM_UNIT",
        "CONC_UNIT",
        "BETA",
        "ALPHA",
        "TAU1",
        "TAU2",
        "TAU3",
        "NAT_EM",
    ]
    assert len(df_gas.index) == 46
