import os

import numpy as np

from ciceroscm import input_handler
from ciceroscm.parallel.cscmparwrapper import CSCMParWrapper, run_ciceroscm_parallel


def test_cscmparwrapper(test_data_dir):

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
        "nyend": 2050,
        "nystart": 1900,
        "emstart": 2015,
        "concentrations_data": conc_data,
        "nat_ch4_data": nat_ch4_data,
        "nat_n2o_data": nat_n2o_data,
        "emissions_data": em_data,
        "udir": test_data_dir,
        "scenname": "ssp245",
    }
    parwrapper = CSCMParWrapper(scendata)

    assert parwrapper.scen == "ssp245"
    assert parwrapper.model == "ssp245"
    cfgs = [
        {
            "pamset_udm": {
                "rlamdo": 15.1,
                "akapa": 0.657,
                "cpi": 0.208,
                "W": 2.2,
                "beto": 6.9,
                "lambda": 0.606,
                "mixed": 107.0,
            },
            "pamset_emiconc": {
                "qbmb": 0.0,
                "qo3": 0.5,
                "qdirso2": -0.3701,
                "qindso2": -0.4163,
                "qbc": 0.163,
                "qoc": -0.084,
                "qh2o_ch4": 0.171,
            },
            "Index": "13555_old_NR_rounded",
        }
    ]
    output_variables = [
        "Heat Content|Ocean",
        "Surface Air Temperature Change",
        "Effective Radiative Forcing|Anthropogenic",
        "Effective Radiative Forcing|Greenhouse Gases",
        "Emissions|CH4",
        "Atmospheric Concentrations|N2O",
        "Ocean carbon flux",
        "Airborne fraction CO2",
        "Biosphere carbon pool",
    ]
    results = parwrapper.run_over_cfgs(cfgs, output_variables)
    print(results)
    assert set(results["variable"].unique()) == set(
        [
            "Heat Content|Ocean",
            "Surface Air Temperature Change",
            "Effective Radiative Forcing|Anthropogenic",
            "Effective Radiative Forcing|Greenhouse Gases",
            "Emissions|CH4",
            "Atmospheric Concentrations|N2O",
            "Ocean carbon flux",
            "Airborne fraction CO2",
            "Biosphere carbon pool",
        ]
    )
    assert set(results["scenario"].unique()) == set(["ssp245"])
    test_length = results.query(
        'variable=="Surface Air Temperature Change" & scenario=="ssp245" & run_id=="13555_old_NR_rounded"'
    )
    test_not_empty = results[2000]
    assert len(test_length.values[0]) == 158
    assert test_not_empty.values.any()


def test_ciceroscm_run_parallel_many_scenarios(test_data_dir):
    scenarios = []
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
    cfgs = [
        {
            "pamset_carbon": {
                "npp0": 60,
                "solubility_sens": 0,
                "ml_t_half": 10,
                "t_half": 10,
                "t_threshold": 10,
            },
            "Index": "13555_old_NR_rounded",
        }
    ]
    for i in range(20):
        em_data_here = em_data.mul((1 + i / 100))
        new_scen = {
            "gaspam_data": gaspam_data,
            "nyend": 2050,
            "nystart": 1900,
            "emstart": 2015,
            "concentrations_data": conc_data,
            "nat_ch4_data": nat_ch4_data,
            "nat_n2o_data": nat_n2o_data,
            "emissions_data": em_data_here,
            "udir": test_data_dir,
            "scenname": "ssp245-plus-%d-percent" % i,
        }
        scenarios.append(new_scen)
    output_variables = ["Heat Content|Ocean", "Surface Air Temperature Change"]
    results = run_ciceroscm_parallel(scenarios, cfgs, output_variables)
    assert set(results["variable"].unique()) == set(
        ["Heat Content|Ocean", "Surface Air Temperature Change"]
    )
    assert set(results["scenario"].unique()) == set(
        [f"ssp245-plus-{d}-percent" for d in range(20)]
    )
    test_length = results.query(
        'variable=="Surface Air Temperature Change" & scenario=="ssp245-plus-1-percent" & run_id=="13555_old_NR_rounded"'
    )
    test_not_empty = results[2000]

    assert len(test_length.values[0]) == 158

    assert test_not_empty.values.any()


def test_ciceroscm_run_parallel_many_cfgs(test_data_dir):
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
    scenarios = [
        {
            "gaspam_data": gaspam_data,
            "nyend": 2050,
            "nystart": 1900,
            "emstart": 2015,
            "concentrations_data": conc_data,
            "nat_ch4_data": nat_ch4_data,
            "nat_n2o_data": nat_n2o_data,
            "emissions_data": em_data,
            "udir": test_data_dir,
            "scenname": "ssp245",
        }
    ]
    cfgs = [
        {
            "pamset_udm": {
                "rlamdo": 15.1,
                "akapa": 0.657,
                "cpi": 0.208,
                "W": 2.2,
                "beto": 6.9,
                "lambda": 0.606,
                "mixed": 107.0,
            },
            "pamset_emiconc": {
                "qbmb": 0.0,
                "qo3": 0.5,
                "qdirso2": -0.3701,
                "qindso2": -0.4163,
                "qbc": 0.163,
                "qoc": -0.084,
                "qh2o_ch4": 0.171,
            },
            "pamset_carbon": {
                "npp0": 60,
                "solubility_sens": 0,
                "ml_t_half": 10,
                "t_half": 10,
                "t_threshold": 10,
            },
            "Index": "13555_old_NR_rounded",
        },
        {
            "pamset_udm": {
                "rlamdo": 15.08357,
                "akapa": 0.6568376339229769,
                "cpi": 0.2077266,
                "W": 2.205919,
                "beto": 6.89822,
                "lambda": 0.6062529,
                "mixed": 107.2422,
            },
            "pamset_emiconc": {
                "qbmb": 0.0,
                "qo3": 0.5,
                "qdirso2": -0.3562,
                "qindso2": -0.96609,
                "qbc": 0.1566,
                "qoc": -0.0806,
            },
            "pamset_carbon": {
                "npp0": 60,
                "solubility_sens": 0,
                "ml_t_half": 10,
                "t_half": 10,
                "t_threshold": 10,
            },
            "Index": "13555_old_NR_improved",
        },
        {
            "pamset_udm": {
                "rlamdo": 15.08357,
                "akapa": 0.6568376339229769,
                "cpi": 0.2077266,
                "W": 2.205919,
                "beto": 6.89822,
                "lambda": 0.6062529,
                "mixed": 107.2422,
            },
            "pamset_emiconc": {
                "qbmb": 0.0,
                "qo3": 0.5,
                "qdirso2": -0.3701323367808028 / 1.5,
                "qindso2": -0.4162980444986502 * 1.5,
                "qbc": 0.162692563111132,
                "qoc": -0.08377713183167902,
            },
            "pamset_carbon": {
                "npp0": 60,
                "solubility_sens": 0,
                "ml_t_half": 10,
                "t_half": 10,
                "t_threshold": 10,
            },
            "Index": "13555_old_NR",
        },
        {
            "pamset_udm": {
                "rlamdo": 5.269455,
                "akapa": 0.40099950002568496,
                "cpi": 0.2054687,
                "W": 1.95183,
                "beto": 3.278654,
                "lambda": 0.7308369,
                "mixed": 115.1219,
            },
            "pamset_emiconc": {
                "qbmb": 0.0,
                "qo3": 0.3,
                "qdirso2": -0.32211222516087934,
                "qindso2": -0.3622885009449893,
                "qbc": 0.1415852070009452,
                "qoc": -0.07290808089507649,
            },
            "pamset_carbon": {
                "npp0": 60,
                "solubility_sens": 0,
                "ml_t_half": 10,
                "t_half": 10,
                "t_threshold": 10,
            },
            "Index": "10496_old_NR",
        },
        {
            "pamset_udm": {
                "rlamdo": 5.269455,
                "akapa": 0.40099950002568496,
                "cpi": 0.2054687,
                "W": 1.95183,
                "beto": 3.278654,
                "lambda": 0.7308369,
                "mixed": 115.1219,
            },
            "pamset_emiconc": {
                "qbmb": 0.0,
                "qo3": 0.3,
                "qdirso2": -0.32211222516087934,
                "qindso2": -0.3622885009449893,
                "qbc": 0.1415852070009452,
                "qoc": -0.07290808089507649,
            },
            "pamset_carbon": {
                "npp0": 60,
                "solubility_sens": 0,
                "ml_t_half": 10,
                "t_half": 10,
                "t_threshold": 10,
            },
            "Index": "10974_old_NR",
        },
        {
            "pamset_udm": {
                "rlamdo": 23.71469,
                "akapa": 0.5486541129586187,
                "cpi": 0.4059296,
                "W": 2.090927,
                "beto": 4.426507,
                "lambda": 0.6345007,
                "mixed": 97.25478,
            },
            "pamset_emiconc": {
                "qbmb": 0.0,
                "qo3": 0.3,
                "qdirso2": -0.4031202924043847,
                "qindso2": -0.453400507735301,
                "qbc": 0.17719249872571508,
                "qoc": -0.09124374858602939,
            },
            "pamset_carbon": {
                "npp0": 60,
                "solubility_sens": 0,
                "ml_t_half": 10,
                "t_half": 10,
                "t_threshold": 10,
            },
            "Index": "Unknown_old_NR",
        },
        {
            "pamset_udm": {
                "rlamdo": 23.13088,
                "akapa": 0.6588532950589802,
                "cpi": 0.1690946,
                "W": 2.228695,
                "beto": 2.623041,
                "lambda": 0.5402487,
                "mixed": 99.86714,
            },
            "pamset_emiconc": {
                "qbmb": 0.0,
                "qo3": 0.3,
                "qdirso2": -0.882827534329248,
                "qindso2": -0.9929404692583796,
                "qbc": 0.3880489761967636,
                "qoc": -0.1998224726091362,
            },
            "pamset_carbon": {
                "npp0": 60,
                "solubility_sens": 0,
                "ml_t_half": 10,
                "t_half": 10,
                "t_threshold": 10,
            },
            "Index": "28925_old_NR",
        },
    ]

    output_variables = ["Heat Content|Ocean", "Surface Air Temperature Change"]
    results = run_ciceroscm_parallel(scenarios, cfgs, output_variables)
    print(results)
    assert set(results["run_id"].unique()) == set(
        [
            "13555_old_NR_rounded",
            "13555_old_NR_improved",
            "13555_old_NR",
            "10496_old_NR",
            "10974_old_NR",
            "Unknown_old_NR",
            "28925_old_NR",
        ]
    )


def test_ciceroscm_run_parallel_many_forcing(test_data_dir):
    gaspam_data = input_handler.read_components(
        os.path.join(test_data_dir, "gases_v1RCMIP.txt")
    )
    scenarios = [
        {
            "gaspam_data": gaspam_data,
            "nyend": 2100,
            "forc_data": np.loadtxt(os.path.join(test_data_dir, "test_forcing.txt")),
            "udir": test_data_dir,
            "scenname": "forc_data_test_forcing",
        },
        {
            "gaspam_data": gaspam_data,
            "nyend": 2100,
            "forc_data": np.loadtxt(os.path.join(test_data_dir, "zero_forcing.txt")),
            "udir": test_data_dir,
            "scenname": "forc_data_zero_forcing",
        },
        {
            "gaspam_data": gaspam_data,
            "nyend": 2100,
            "forc_data": np.loadtxt(os.path.join(test_data_dir, "test_forcing.txt")),
            "sunvolc": 1,
            "udir": test_data_dir,
            "scenname": "forc_data_test_forcing_sunvolc",
        },
        {
            "gaspam_data": gaspam_data,
            "nyend": 2100,
            "forc_data": np.loadtxt(os.path.join(test_data_dir, "CO2_1pros.txt")),
            "udir": test_data_dir,
            "scenname": "forc_data_CO2_1pros",
        },
        {
            "gaspam_data": gaspam_data,
            "nyend": 2100,
            "forc_file": os.path.join(test_data_dir, "test_forcing_hemisplit.txt"),
            "udir": test_data_dir,
            "scenname": "forc_data_test_forcing_hemisplit",
        },
    ]
    cfgs = [
        {
            "pamset_udm": {
                "rlamdo": 15.1,
                "akapa": 0.657,
                "cpi": 0.208,
                "W": 2.2,
                "beto": 6.9,
                "lambda": 0.606,
                "mixed": 107.0,
            },
            "pamset_emiconc": {
                "qbmb": 0.0,
                "qo3": 0.5,
                "qdirso2": -0.3701,
                "qindso2": -0.4163,
                "qbc": 0.163,
                "qoc": -0.084,
                "qh2o_ch4": 0.171,
            },
            "pamset_carbon": {
                "npp0": 60,
                "solubility_sens": 0,
                "ml_t_half": 10,
                "t_half": 10,
                "t_threshold": 10,
            },
            "Index": "10496_old_NR_rounded",
        }
    ]
    output_variables = ["Heat Content|Ocean", "Surface Air Temperature Change"]
    results = run_ciceroscm_parallel(scenarios, cfgs, output_variables)
    print(results)
    assert set(results["scenario"].unique()) == set(
        [
            "forc_data_test_forcing_hemisplit",
            "forc_data_CO2_1pros",
            "forc_data_test_forcing_sunvolc",
            "forc_data_zero_forcing",
            "forc_data_test_forcing",
        ]
    )
    test_length = results.query(
        'variable=="Surface Air Temperature Change" & scenario=="forc_data_zero_forcing" & run_id=="10496_old_NR_rounded"'
    )
    assert len(test_length.values[0]) == 358


def test_parallel_conc_run(test_data_dir):
    scenarios = [
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "conc_run": True,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "sunvolc": 1,
            "scenname": "ssp245conc",
        }
    ]

    cfgs = [
        {
            "pamset_emiconc": {
                "qh2o_ch4": 0.171,
                "qbmb": 0.03,
                "qo3": 0.4,
                "qdirso2": -0.457 / 57,
                "qindso2": -0.514 / 57,
                "qbc": 0.200 / 7,
                "qoc": -0.103 / 16,
            },
            "pamset_udm": {
                "rlamdo": 16.0,
                "akapa": 0.634,
                "cpi": 0.4,
                "W": 4.0,
                "beto": 3.5,
                "lambda": 0.540,
                "mixed": 60.0,
            },
            "Index": "test_test",
        },
        {
            "pamset_emiconc": {
                "qh2o_ch4": 0.171,
                "qbmb": 0.03,
                "qo3": 0.4,
                "qdirso2": -0.457 / 57,
                "qindso2": -0.514 / 57,
                "qbc": 0.200 / 7,
                "qoc": -0.103 / 16,
            },
            "pamset_udm": {
                "rlamdo": 16.0,
                "akapa": 0.634,
                "cpi": 0.4,
                "W": 4.0,
                "beto": 3.5,
                "lambda": 0.540,
                "mixed": 60.0,
                "ocean_efficacy": 1.2,
            },
            "pamset_carbon": {"mixed_carbon": 120},
            "Index": "test_test_oceff",
        },
    ]

    output_vars = [
        "Surface Air Temperature Change",
        "Effective Radiative Forcing",
        "Emissions|CO2",
        "Heat Content|Ocean",
        "Atmospheric Concentrations|CO2",
        "Effective Radiative Forcing|CO2",
        "Effective Radiative Forcing|N2O",
        "Effective Radiative Forcing|CH4",
        "Effective Radiative Forcing|Aerosols",
        "Effective Radiative Forcing|F-Gases",
        "Biosphere carbon flux",
        "Airborne fraction CO2",
        "Ocean carbon flux",
    ]

    results = run_ciceroscm_parallel(scenarios, cfgs, output_vars)
    print(results.columns)
    assert set(results["variable"].unique()) == set(output_vars)
    print(results["run_id"].unique())
    back_em1 = (
        results.loc[
            (results["variable"] == "Emissions|CO2")
            & (results["run_id"] == "test_test")
        ]
        .iloc[:, 7:]
        .to_numpy(float)
    )
    print(back_em1)
    back_em2 = (
        results.loc[
            (results["variable"] == "Emissions|CO2")
            & (results["run_id"] == "test_test_oceff")
        ]
        .iloc[:, 7:]
        .to_numpy(float)
    )
    assert not np.allclose(back_em1, back_em2)
