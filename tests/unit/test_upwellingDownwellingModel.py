import pytest
import numpy as np


from ciceroscm import upwelling_diffusion_model


def test_coefic():
    assert upwelling_diffusion_model._coefic(0, 0, 0) == 19652.21
    assert upwelling_diffusion_model._coefic(1, 0, 0) == 19706.96404
    assert upwelling_diffusion_model._coefic(0, 1, 0) == 19798.31704821712
    assert upwelling_diffusion_model._coefic(0, 0, 1) == 19655.449993093498


def test_denso():
    assert upwelling_diffusion_model._denso(0, 0) == 999.842594
    assert upwelling_diffusion_model._denso(0, 1) == 999.9015372849533
    assert upwelling_diffusion_model._denso(1, 0) == 1000.6618454799999


def test_density():
    assert upwelling_diffusion_model._density(0, 0) == 1028.1063314148107
    assert upwelling_diffusion_model._density(1, 0) == 1028.1539628244798
    assert upwelling_diffusion_model._density(0, 1) == 1028.0456085605597


# Rewrite test for class method
def test_band():
    udm = upwelling_diffusion_model.UpwellingDiffusionModel(
        {
            "lambda": 0.540,
            "akapa": 0.341,
            "cpi": 0.556,
            "W": 1.897,
            "rlamdo": 16.618,
            "beto": 3.225,
            "mixed": 107.277,
            "dirso2_forc": -0.457,
            "indso2_forc": -0.514,
            "bc_forc": 0.200,
            "oc_forc": -0.103,
            "threstemp": 7.0,
        }
    )

    testa = np.ones(40)
    testb = np.array(
        [
            1,
            3,
            7,
            5,
            8,
            9,
            1,
            3,
            7,
            5,
            8,
            9,
            1,
            3,
            7,
            5,
            8,
            9,
            1,
            3,
            7,
            5,
            8,
            9,
            1,
            3,
            7,
            5,
            8,
            9,
            1,
            3,
            7,
            5,
            8,
            9,
            7,
            5,
            8,
            9,
        ]
    )
    testc = np.array(
        [
            2,
            3,
            2,
            2,
            3,
            2,
            2,
            3,
            2,
            2,
            3,
            2,
            2,
            3,
            2,
            2,
            3,
            2,
            2,
            3,
            2,
            2,
            3,
            22,
            3,
            2,
            2,
            3,
            2,
            2,
            3,
            2,
            2,
            3,
            2,
            2,
            2,
            3,
            2,
        ]
    )
    testd = np.array(
        [
            8,
            3,
            7,
            11,
            8,
            3,
            7,
            11,
            8,
            3,
            7,
            11,
            8,
            3,
            7,
            11,
            8,
            3,
            7,
            11,
            8,
            3,
            7,
            11,
            8,
            3,
            7,
            11,
            8,
            3,
            7,
            11,
            8,
            3,
            7,
            11,
            8,
            3,
            7,
            11,
        ]
    )

    ans = udm._band(testa, testb, 1.5 * testc, 0.5 * testd)
    assert len(ans) == len(testa)
    assert ans[0] != 0
    assert ans[-1] != 0
    print(ans)
    np.testing.assert_allclose(
        ans,
        [
            7.29695905,
            -3.29695905,
            0.90975958,
            0.14288067,
            1.29194569,
            -1.43965471,
            4.38831556,
            0.18377972,
            0.12452118,
            0.98152402,
            -1.17738042,
            2.65278208,
            -5.73255276,
            2.35992356,
            0.03395157,
            0.30080514,
            1.32067424,
            -1.52582201,
            4.63724129,
            0.12952691,
            0.10537289,
            1.04428763,
            -1.27560368,
            2.81345373,
            -0.56198424,
            0.38856234,
            0.29876574,
            0.34002582,
            0.77802337,
            -0.85473759,
            2.80487166,
            0.34441465,
            0.55396146,
            -0.07404829,
            0.29250667,
            0.41133165,
            0.50183616,
            0.02527175,
            0.19373447,
            0.6416175,
        ],
        rtol=1.0e-5,
    )
