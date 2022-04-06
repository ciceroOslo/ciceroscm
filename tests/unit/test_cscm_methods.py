import os

import numpy as np
import pandas as pd

from ciceroscm import ciceroscm


def test_read_volc_sun(test_data_dir):
    cscm = ciceroscm.CICEROSCM(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "forc_file": os.path.join(test_data_dir, "CO2_1pros.txt"),
        }
    )

    # Testing volcano rf reading
    df_volc = cscm.read_data_on_year_row(
        os.path.join(test_data_dir, "meanVOLCmnd_ipcc_NH.txt")
    )
    print(df_volc.columns.tolist)
    assert len(df_volc.index) == 351
    assert df_volc.columns.tolist() == list(range(12))

    # Testing solar rf reading:
    df_sun = cscm.read_data_on_year_row(os.path.join(test_data_dir, "solar_IPCC.txt"))
    assert len(df_sun.index) == 351
    assert df_sun.columns.tolist() == [0]

    # Testing LUCalbedo rf reading
    df_luc = cscm.read_data_on_year_row(
        os.path.join(test_data_dir, "IPCC_LUCalbedo.txt")
    )
    assert len(df_luc.index) == 351
    assert df_luc.columns.tolist() == [0]


def test_read_forc(test_data_dir):
    forc = ciceroscm.read_forc(os.path.join(test_data_dir, "CO2_1pros.txt"))
    assert isinstance(forc, np.ndarray)
    assert len(forc) == 736

    forc = ciceroscm.read_forc(os.path.join(test_data_dir, "AR6_ERF_1750-2019.csv"))
    assert isinstance(forc, pd.DataFrame)
    total = forc["total"].to_numpy()
    assert len(total) == 270
    assert total[1] == 0.2866421077706115
