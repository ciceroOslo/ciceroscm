import os
import shutil

import numpy as np
import pandas as pd
import pandas.testing as pdt

from ciceroscm import CICEROSCM, input_handler
from ciceroscm.formattingtools import (
    reformat_cscm_results,
    reformat_inputdata_to_cscm_format,
)


def check_output(
    output_dir,
    expected_output_dir,
    update_expected_files=False,
    rtol=1e-2,
    files=["output_em.txt", "output_conc.txt"],
):
    for filename in files:
        file_to_check = os.path.join(output_dir, filename)
        file_expected = os.path.join(expected_output_dir, filename)
        if update_expected_files:
            shutil.copyfile(file_to_check, file_expected)
        else:
            res = pd.read_csv(file_to_check, sep=r"\s+")
            exp = pd.read_csv(file_expected, sep=r"\s+")
            pdt.assert_index_equal(res.index, exp.index)

            pdt.assert_frame_equal(
                res.T,
                exp.T,
                check_like=True,
                rtol=rtol,
            )


def check_output_just_some_lines(
    output_dir,
    expected_output_dir,
    update_expected_files=False,
    rtol=1e-2,
    files=["output_temp.txt", "output_forc.txt"],
    lines=10,
):
    for filename in files:
        file_to_check = os.path.join(output_dir, filename)
        file_expected = os.path.join(expected_output_dir, filename)

        if update_expected_files:
            shutil.copyfile(file_to_check, file_expected)
        else:
            res = pd.read_csv(file_to_check, sep=r"\s+", skiprows=range(lines, 352))
            exp = pd.read_csv(file_expected, sep=r"\s+", skiprows=range(lines, 352))
            pdt.assert_index_equal(res.index, exp.index)

            pdt.assert_frame_equal(
                res.T,
                exp.T,
                check_like=True,
                rtol=rtol,
            )


def test_ciceroscm_run_emi(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            # "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:

    cscm._run(
        {"output_folder": outdir},
        pamset_udm={
            "rlamdo": 15.1,
            "akapa": 0.657,
            "cpi": 0.208,
            "W": 2.2,
            "beto": 6.9,
            "lambda": 0.606,
            "mixed": 107.0,
        },
        pamset_emiconc={
            "qbmb": 0.0,
            "qo3": 0.5,
            "qdirso2": -0.3701 / 57,
            "qindso2": -0.4163 / 57,
            "qbc": 0.163 / 7,
            "qoc": -0.084 / 16,
            "qh2o_ch4": 0.171,
        },
    )

    check_output(outdir, os.path.join(test_data_dir, "ssp245_emis"))
    check_output_just_some_lines(
        outdir,
        os.path.join(test_data_dir, "ssp245_emis"),
        files=["output_forc.txt"],
        lines=14,
    )
    check_output_just_some_lines(
        outdir,
        os.path.join(test_data_dir, "ssp245_emis"),
        files=["output_temp.txt"],
        lines=16,
    )


def test_ciceroscm_short_run(tmpdir, test_data_dir):
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:
    nystart = 1900
    nyend = 2050
    emstart = 1950
    gaspamfile = os.path.join(test_data_dir, "gases_v1RCMIP.txt")
    cscm = CICEROSCM(
        {
            "gaspam_file": gaspamfile,
            "nystart": nystart,
            "emstart": emstart,
            "nyend": nyend,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "idtm": 24,
            "sunvolc": 1,
            "rf_luc_data": pd.read_csv(
                os.path.join(test_data_dir, "land_use_erf_ar6.txt"),
                header=None,
                skiprows=1,
                index_col=0,
            ),
            "rf_sun_data": pd.read_csv(
                os.path.join(test_data_dir, "solar_erf_ar6.txt"),
                header=None,
                skiprows=1,
                index_col=0,
            ),
            "rf_volc_data": pd.read_csv(
                os.path.join(test_data_dir, "volcanic_erf_ar6.txt"),
                header=None,
                skiprows=1,
                index_col=0,
            ),
        },
    )

    cscm._run({"output_folder": outdir})

    file_results = os.path.join(outdir, "output_em.txt")
    exp_index = np.arange(nystart, nyend + 1)
    res = pd.read_csv(file_results, sep=r"\s+")
    np.testing.assert_equal(res.Year.to_numpy(), exp_index)

    cscm._run({"results_as_dict": True, "carbon_cycle_outputs": True})
    expected_keys = [
        "emissions",
        "concentrations",
        "forcing",
        "OHC700",
        "OHCTOT",
        "RIB_glob",
        "RIB_N",
        "RIB_S",
        "dT_glob",
        "dT_NH",
        "dT_SH",
        "dT_glob_air",
        "dT_NH_air",
        "dT_SH_air",
        "dT_glob_sea",
        "dT_NH_sea",
        "dT_SHsea",
        "Total_forcing",
        "Solar_forcing",
        "Volcanic_forcing_NH",
        "Volcanic_forcing_SH",
        "carbon cycle",
    ]
    for key in expected_keys:
        assert key in cscm.results
    for key in cscm.results.keys():
        assert key in expected_keys
    assert len(cscm.results["carbon cycle"]["Airborne fraction CO2"]) == len(
        cscm.ce_handler.years
    )
    assert len(cscm.results["carbon cycle"]["Biosphere carbon flux"]) == len(
        cscm.ce_handler.years
    )
    assert len(cscm.results["carbon cycle"]["Ocean carbon flux"]) == len(
        cscm.ce_handler.years
    )
    carbon_sum = (
        np.cumsum(
            cscm.results["carbon cycle"]["Airborne fraction CO2"]
            * cscm.results["emissions"]["CO2"]
        )
        + cscm.results["carbon cycle"]["Biosphere carbon flux"].values
        + cscm.results["carbon cycle"]["Ocean carbon flux"].values
    )
    print(carbon_sum[-5:])
    print(np.cumsum(cscm.results["emissions"]["CO2"])[-5:])
    # TODO: Why isn't this closer?
    assert np.allclose(
        carbon_sum, np.cumsum(cscm.results["emissions"]["CO2"]), rtol=5e-1
    )

    # Block for testing reformatting of outputs
    sfilewriter = reformat_inputdata_to_cscm_format.COMMONSFILEWRITER(gaspamfile)
    cscm_reader = reformat_cscm_results.CSCMREADER(nystart, nyend)
    print(cscm.results.keys())
    test_variables = [
        "Ocean carbon flux",
        "Airborne fraction CO2",
        "Biosphere carbon pool",
        "Emissions|CH4",
        "Atmospheric Concentrations|CO2",
        "Heat Content|Ocean",
        "Effective Radiative Forcing|Aerosols|Direct Effect",
        "Effective Radiative Forcing|Anthropogenic",
        "Effective Radiative Forcing",
        "Effective Radiative Forcing|F-Gases",
        "Effective Radiative Forcing|Greenhouse Gases",
        "Effective Radiative Forcing|C6F14",
        "Surface Air Ocean Blended Temperature Change",
        "Heat Uptake",
    ]
    unit_list = [
        "Pg C / yr",
        "Unitless",
        "Pg C",
        "TgCH4 / yr",
        "ppm",
        "ZJ",
        "W/m^2",
        "W/m^2",
        "W/m^2",
        "W/m^2",
        "W/m^2",
        "W/m^2",
        "K",
        "W/m^2",
    ]
    for i, variable in enumerate(test_variables):
        format_output = cscm_reader.get_variable_timeseries(
            cscm.results, variable, sfilewriter
        )
        print(variable)
        print(format_output[0])
        assert np.all(format_output[0] == np.arange(nystart, nyend + 1))
        assert len(format_output[0]) == nyend - nystart + 1
        assert format_output[2] == unit_list[i]
    empty_result = cscm_reader.get_variable_timeseries(
        {}, "Ocean carbon flux", sfilewriter
    )
    assert len(empty_result[0]) == 0
    assert len(empty_result[1]) == 0
    assert empty_result[2] == "NoUnit"
    # Now test getting some variables

    # Put this in again, find out what is happening with CF4
    # check_output(
    #    outdir_save,
    #    os.path.join(test_data_dir, "ssp245_emis"),
    #    files=["output_temp.txt"],
    #    rtol=1.0,
    # )

    # check_output(
    #    outdir_save,
    #    os.path.join(test_data_dir, "ssp245_emis"),
    #    files=["output_forc.txt"],
    #    rtol=0.1,
    # )


def test_ciceroscm_run_conc(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "conc_run": True,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "sunvolc": 1,
        },
    )
    outdir = str(tmpdir)
    # outdir_save = os.path.join(os.getcwd(), "output")

    # One year forcing:

    cscm._run(
        {"output_folder": outdir},
        pamset_emiconc={
            "qh2o_ch4": 0.171,
            "qbmb": 0.03,
            "qo3": 0.4,
            "qdirso2": -0.457 / 57,
            "qindso2": -0.514 / 57,
            "qbc": 0.200 / 7,
            "qoc": -0.103 / 16,
        },
        pamset_udm={
            "rlamdo": 16.0,
            "akapa": 0.634,
            "cpi": 0.4,
            "W": 4.0,
            "beto": 3.5,
            "lambda": 0.540,
            "mixed": 60.0,
        },
    )

    check_output(
        outdir,
        os.path.join(test_data_dir, "ssp245_conc"),
        files=["output_conc.txt", "output_em.txt", "output_forc.txt", "output_ohc.txt"],
    )


def test_run_with_data_not_files(tmpdir, test_data_dir):
    ih = input_handler.InputHandler({"nystart": 1750, "nyend": 2100, "emstart": 1850})
    cscm = CICEROSCM(
        {
            "gaspam_data": input_handler.read_components(
                os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt")
            ),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_data": input_handler.read_inputfile(
                os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"), True, 1750, 2100
            ),
            "emissions_data": ih.read_emissions(
                os.path.join(test_data_dir, "ssp245_em_RCMIP.txt")
            ),
            "nat_ch4_data": input_handler.read_natural_emissions(
                os.path.join(test_data_dir, "natemis_ch4.txt"), "CH4"
            ),
            "nat_n2o_data": input_handler.read_natural_emissions(
                os.path.join(test_data_dir, "natemis_n2o.txt"), "N2O"
            ),
        },
    )
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:

    cscm._run(
        {"output_folder": outdir},
        pamset_udm={
            "rlamdo": 15.1,
            "akapa": 0.657,
            "cpi": 0.208,
            "W": 2.2,
            "beto": 6.9,
            "lambda": 0.606,
            "mixed": 107.0,
        },
        pamset_emiconc={
            "qbmb": 0.0,
            "qo3": 0.5,
            "qdirso2": -0.3701 / 57,
            "qindso2": -0.4163 / 57,
            "qbc": 0.163 / 7,
            "qoc": -0.084 / 16,
            "qh2o_ch4": 0.171,
        },
    )

    check_output(outdir, os.path.join(test_data_dir, "ssp245_emis"))
    check_output_just_some_lines(
        outdir,
        os.path.join(test_data_dir, "ssp245_emis"),
        files=["output_forc.txt"],
        lines=14,
    )
    check_output_just_some_lines(
        outdir,
        os.path.join(test_data_dir, "ssp245_emis"),
        files=["output_temp.txt"],
        lines=16,
    )


def test_ciceroscm_just_one(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "conc_run": True,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
        },
    )
    # outdir = str(tmpdir)

    # One year forcing:

    cscm._run(
        {"results_as_dict": True},
        pamset_emiconc={"qh2o_ch4": 0.171, "just_one": "CO2"},
    )
    assert cscm.results["forcing"]["CO2"].equals(
        cscm.results["forcing"]["Total_forcing"]
    )

    cscm2 = CICEROSCM(
        {"nyend": 2100, "forc_data": cscm.results["forcing"]["CO2"].to_numpy()}
    )

    cscm2._run({"results_as_dict": True})

    for key in cscm2.results:
        assert np.array_equal(cscm2.results[key], cscm.results[key])
