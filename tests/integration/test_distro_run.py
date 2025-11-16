import json
import os

import numpy as np
import pytest
from pandas.testing import assert_frame_equal

from ciceroscm import input_handler
from ciceroscm.parallel._configdistro import _ConfigDistro
from ciceroscm.parallel.distributionrun import DistributionRun


def test_distro_run_forcing(test_data_dir):
    testconfig = _ConfigDistro(
        setvalues={"threstemp": 7.0, "lm": 40, "ldtime": 12}, options={"forc": True}
    )
    jsonfile = os.path.join(test_data_dir, "test_config_json.json")

    distrorun1 = DistributionRun(testconfig, numvalues=10)
    gaspam_data = input_handler.read_components(
        os.path.join(test_data_dir, "gases_v1RCMIP.txt")
    )
    scenarios = {
        "gaspam_data": gaspam_data,
        "nyend": 2100,
        "forc_data": np.loadtxt(os.path.join(test_data_dir, "test_forcing.txt")),
        "udir": test_data_dir,
        "scenname": "forc_data_test_forcing",
    }
    output_vars = ["Heat Content|Ocean", "Surface Air Temperature Change"]

    # Testing run from json file without metadata
    distrorun1.write_configs_to_json(json_file_name=jsonfile)
    results1 = distrorun1.run_over_distribution(scenarios, output_vars)
    distrorun2 = DistributionRun(None, jsonfile)
    os.remove(jsonfile)
    results2 = distrorun2.run_over_distribution(scenarios, output_vars)
    assert_frame_equal(results1, results2)

    # Testing run from json file with metadata
    meta_info = {
        "Model_version": "ciceroscm-v1.5.0",
        "carbon_cycle": "Default",
        "thermal_model": "Default",
        "dump_date": "2025-10-08",
    }
    distrorun1.write_configs_to_json(json_file_name=jsonfile, meta_info=meta_info)
    distrorun3 = DistributionRun(None, jsonfile)
    with open(jsonfile, "r", encoding="utf-8") as rfile:
        json_contents = json.load(rfile)

    assert set(json_contents.keys()) == set(["configurations", "meta_info"])
    assert distrorun3.meta_info == meta_info
    os.remove(jsonfile)

    results3 = distrorun3.run_over_distribution(scenarios, output_vars)
    assert_frame_equal(results1, results3)

    testconfig.make_config_lists(10, json_fname=jsonfile, indexer_pre="testnometa")
    distrorun2 = DistributionRun(None, jsonfile)
    assert distrorun2.meta_info is None
    os.remove(jsonfile)

    testconfig.make_config_lists(
        10, json_fname=jsonfile, indexer_pre="testmeta", meta_info=meta_info
    )
    distrorun2 = DistributionRun(None, jsonfile)
    assert distrorun2.meta_info == meta_info
    os.remove(jsonfile)

    del json_contents["configurations"]
    # Write the current json_contents to jsonfile, and try to create a distrorun4 from that file,
    # this should throw a ValueError which we should test for using pytest.raises
    with open(jsonfile, "w", encoding="utf-8") as wfile:
        json.dump(json_contents, wfile)
    with pytest.raises(ValueError, match=f"{jsonfile} not formatted correctly"):
        distrorun1 = DistributionRun(None, jsonfile)
    os.remove(jsonfile)
