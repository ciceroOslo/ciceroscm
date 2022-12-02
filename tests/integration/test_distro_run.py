import os

import numpy as np

from ciceroscm import input_handler
from ciceroscm.parallel._configdistro import _ConfigDistro
from ciceroscm.parallel.distributionrun import DistributionRun


def test_distro_run_forcing(test_data_dir):
    testconfig = _ConfigDistro(
        setvalues={"threstemp": 7.0, "lm": 40, "ldtime": 12}, options={"forc": True}
    )
    jsonfile = os.path.join(test_data_dir, "test_config_json.json")
    distrorun1 = DistributionRun(testconfig, numvalues=40)
    gaspam_data = input_handler.read_components(
        os.path.join(test_data_dir, "gases_v1RCMIP.txt")
    )
    scenarios = {
        "gaspamfile": gaspam_data,
        "nyend": 2100,
        "forc_data": np.loadtxt(os.path.join(test_data_dir, "test_forcing.txt")),
        "udir": test_data_dir,
        "scenname": "forc_data_test_forcing",
    }
    output_vars = ["Heat Content|Ocean", "Surface Air Temperature Change"]
    distrorun1.write_configs_to_json(json_file_name=jsonfile)
    results1 = distrorun1.run_over_distribution(scenarios, output_vars)
    distrorun2 = DistributionRun(None, jsonfile)
    results2 = distrorun2.run_over_distribution(scenarios, output_vars)
    assert np.array_equal(results1.values, results2.values)
    os.remove(jsonfile)
