import numpy as np

from ciceroscm.parallel._configdistro import _ConfigDistro


def test_config_distro_methods():
    config_full = _ConfigDistro(
        distro_array=np.array([[0, 1], [0, 1], [0, 1], [4, 6]]),
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
                [0, 1],
                [0, 1],
                [0, 1],
                [4, 6],
                [0, 7],
                [2 / 3.71, 5 / 3.71],
                [25, 125],
                [0.1, 0.2],
                [-0.1, -0.06],
            ]
        ),
    )
    assert config_full.ordering == [
        "rlamdo",
        "akapa",
        "cpi",
        "W",
        "beto",
        "lambda",
        "mixed",
        "qbc",
        "qoc",
    ]

    num_latin = 100
    latin_samples = config_full.get_samples_from_distro_latin(num_latin)
    gaussian_samples = config_full.get_samples_from_distro_gaussian(50)
    assert latin_samples.shape == (100, len(config_full.ordering))
    assert gaussian_samples.shape == (50, len(config_full.ordering))
    for i, interval in enumerate(config_full.prior):
        masked = np.ma.masked_inside(latin_samples[:, i], interval[0], interval[1])
        assert np.ma.allequal(masked, [True] * num_latin)

    config_list = config_full.make_config_list(20, indexer_pre="test")
    assert len(config_list) == 20
    assert config_list[4]["Index"] == "test4"
    expected_emiconc = ["qbc", "qoc", "qbmb", "qdirso2", "qindso2", "qo3", "qh2o_ch4"]
    expected_udm = [
        "rlamdo",
        "akapa",
        "cpi",
        "W",
        "beto",
        "lambda",
        "mixed",
        "lm",
        "ldtime",
        "threstemp",
    ]
    assert all(pam in config_list[7]["pamset_udm"] for pam in expected_udm)
    print(config_list[13]["pamset_emiconc"])
    assert all(pam in config_list[13]["pamset_emiconc"] for pam in expected_emiconc)
