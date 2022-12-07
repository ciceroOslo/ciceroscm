import os

import pandas as pd

from ciceroscm import input_handler
from ciceroscm.parallel._configdistro import _ConfigDistro
from ciceroscm.parallel.calibrator import Calibrator


def test_calibrator(test_data_dir):
    gaspam_data = input_handler.read_components(
        os.path.join(test_data_dir, "gases_v1RCMIP.txt")
    )
    conc_data = input_handler.read_inputfile(
        os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"), True, 1750, 2100
    )
    ih = input_handler.InputHandler({"nyend": 2050, "nystart": 1900, "emstart": 2015})
    em_data = ih.read_emissions(os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"))
    nat_ch4_data = input_handler.read_natural_emissions(
        os.path.join(test_data_dir, "natemis_ch4.txt"), "CH4"
    )
    nat_n2o_data = input_handler.read_natural_emissions(
        os.path.join(test_data_dir, "natemis_n2o.txt"), "N2O"
    )
    scendata = {
        "gaspam_data": gaspam_data,
        "nyend": 2030,
        "nystart": 1900,
        "emstart": 2015,
        "concentrations_data": conc_data,
        "nat_ch4_data": nat_ch4_data,
        "nat_n2o_data": nat_n2o_data,
        "emissions_data": em_data,
        "udir": test_data_dir,
        "scenname": "ssp245",
    }
    calibdata = pd.DataFrame(
        data={
            "Variable Name": [
                "Heat Content|Ocean",
                "Surface Air Ocean Blended Temperature Change",
            ],
            "Yearstart_norm": [1971, 1961],
            "Yearend_norm": [1971, 1990],
            "Yearstart_change": [2018, 2000],
            "Yearend_change": [2018, 2019],
            "Central Value": [320.69251537323, 0.5372],
            "sigma": [17.020342912051203, 0.039028311931729676],
        }
    )
    # calibdata = pd.DataFrame(data = {'Variable Name': ["Heat Content|Ocean"], 'Yearstart_norm': [1971], "Yearend_norm":[1971], "Yearstart_change":[2018], "Yearend_change":[2018], "Central Value": [320.69251537323], "sigma": [17.020342912051203]})
    testconfig = _ConfigDistro(
        distro_array=[],
        ordering=["aerosol_total", "W", "lambda", "beta_f"],
        setvalues={
            "threstemp": 7.0,
            "lm": 40,
            "ldtime": 12,
            "qbmb": 0,
            "qo3": 0.5,
            "qh2o_ch4": 0.091915,
            "rlamdo": 16,
            "akapa": 0.634,
            "cpi": 0.4,
            "beto": 3.5,
            "mixed": 60,
        },
        options={"aerosol_total": [-0.36, -0.97, 0.16, -0.08]},
    )
    assert testconfig.ordering == ["aerosol_total", "W", "lambda", "beta_f"]
    calibrator = Calibrator(calibdata, testconfig, scendata)
    drawn_cfgs = calibrator.get_n_samples(
        1, current_samples=[], kept_configs=[], recurse_num=0
    )
    assert len(drawn_cfgs) >= 1
    testconfig = _ConfigDistro(
        distro_array=[],
        ordering=["W", "lambda"],
        setvalues={
            "threstemp": 7.0,
            "lm": 40,
            "ldtime": 12,
            "qbmb": 0,
            "qo3": 0.5,
            "qh2o_ch4": 0.091915,
            "rlamdo": 16,
            "akapa": 0.634,
            "cpi": 0.4,
            "beto": 3.5,
            "mixed": 60,
        },
    )
    assert testconfig.ordering == ["W", "lambda"]
    calibrator = Calibrator(calibdata, testconfig, scendata)
    drawn_cfgs = calibrator.get_n_samples(
        3, current_samples=[], kept_configs=[], recurse_num=0
    )
    assert len(drawn_cfgs) >= 3
