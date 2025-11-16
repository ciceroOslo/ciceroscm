import numpy as np

from ciceroscm.parallel._configdistro import _ConfigDistro


def test_config_distro_methods():
    config_full = _ConfigDistro(
        distro_dict={
            "rlamdo": [0, 1],
            "akapa": [0, 1],
            "cpi": [0, 1],
            "W": [4, 6],
            "beta_f": [0.110, 0.465],
        },
        setvalues={
            "qbmb": 0,
            "qo3": 0.5,
            "qdirso2": -0.36,
            "qindso2": -0.97,
            "threstemp": 7.0,
            "lm": 40,
            "ldtime": 12,
            "qh2o_ch4": 0.091915,
        },
    )
    assert np.array_equal(
        config_full.prior,
        np.array(
            [
                [4, 6],
                [0, 1],
                [0.110, 0.465],
                [0, 1],
                [0, 1],
            ]
        ),
    )
    assert config_full.ordering == [
        "W",
        "akapa",
        "beta_f",
        "cpi",
        "rlamdo",
    ]

    num_latin = 100
    latin_samples = config_full.get_samples_from_distro_latin(num_latin)
    gaussian_samples = config_full.get_samples_from_distro_gaussian(50)
    assert latin_samples.shape == (100, len(config_full.ordering))
    assert gaussian_samples.shape == (50, len(config_full.ordering))
    for i, interval in enumerate(config_full.prior):
        masked = np.ma.masked_inside(latin_samples[:, i], interval[0], interval[1])
        assert np.ma.allequal(masked, [True] * num_latin)

    config_list = config_full.make_config_lists(
        20, indexer_pre="test", max_chunk_size=7
    )

    assert len(config_list) == 3
    assert len(config_list[0]) == 7
    assert len(config_list[1]) == 7
    assert len(config_list[2]) == 6
    assert config_list[0][4]["Index"] == "test_0_4"
    expected_emiconc = ["qbmb", "qdirso2", "qindso2", "qo3", "qh2o_ch4"]
    expected_udm = [
        "rlamdo",
        "akapa",
        "cpi",
        "W",
        "lm",
        "ldtime",
        "threstemp",
    ]
    # print(config_list[1][0]["pamset_udm"])
    assert all(pam in config_list[1][0]["pamset_udm"] for pam in expected_udm)
    # print(config_list[1][6]["pamset_emiconc"])
    assert all(pam in config_list[1][6]["pamset_emiconc"] for pam in expected_emiconc)
    assert set(["beta_f"]) == set(config_list[2][4]["pamset_carbon"].keys())
