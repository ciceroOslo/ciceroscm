import os

import numpy as np
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
    distrorun1.write_configs_to_json(json_file_name=jsonfile)
    results1 = distrorun1.run_over_distribution(scenarios, output_vars)
    distrorun2 = DistributionRun(None, jsonfile)
    os.remove(jsonfile)
    results2 = distrorun2.run_over_distribution(scenarios, output_vars)
    assert_frame_equal(results1, results2)
    max_chunk = 5
    distorrun3 = DistributionRun(testconfig, numvalues=10, max_chunk_size=max_chunk)
    results3 = distorrun3.run_over_distribution(scenarios, output_vars)
    compare_index = np.arange(results1.shape[0])
    results3.index = compare_index
    results1.index = compare_index
    print(results1.columns)
    print(results3.columns)
    assert_frame_equal(
        results1.drop(columns=["run_id"]), results3.drop(columns=["run_id"]), atol=100
    )
    print()
    assert all(
        [
            results3["run_id"].values[i]
            == f"_{int(results1['run_id'].values[i])//max_chunk}_{int(results1['run_id'].values[i])%max_chunk}"
            for i in compare_index
        ]
    )
