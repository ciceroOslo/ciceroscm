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
    assert ans[0] == 7.296874116241785
    print(ans)
    np.testing.assert_allclose(
        ans,
        [
            7.296874116241785,
            -3.2968741162417845,
            0.9097218294407932,
            0.14294043671874423,
            1.2918586623218289,
            -1.4395132745096393,
            4.387920269421641,
            0.18386433502933267,
            0.1245526056645247,
            0.9814224751063315,
            -1.1772216603987276,
            2.6525224017963307,
            -5.73182665192275,
            2.359768083375473,
            0.03389386706585125,
            0.3009916157211895,
            1.3203826847760671,
            -1.5253451319843838,
            4.635907834361129,
            0.12981243254108488,
            0.10547885955902586,
            1.0439451835152447,
            -1.2750682590450833,
            2.812577975298983,
            -0.5617616217771444,
            0.388707476995147,
            0.29854639693056784,
            0.3404892481636259,
            0.7775571916114006,
            -0.8536489270182765,
            2.801761050517696,
            0.34486397255568446,
            0.5545490106050834,
            -0.0755690155970894,
            0.29406579275119193,
            0.40768089119585127,
            0.5122687288287155,
            0.002146002334379754,
            0.21711139099986348,
            0.5869876232222374,
        ],
        rtol=1.0e-5,
    )
