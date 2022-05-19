import os

from ciceroscm import ciceroscm


def test_read_volc_sun(test_data_dir):
    cscm = ciceroscm.CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
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
