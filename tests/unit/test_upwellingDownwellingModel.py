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


#Rewrite test for class method    
def test_band():
    testa = np.array([1,1,1])
    testb = np.array([1,3,7])
    testc = np.array([2,3,2])
    testd = np.array([8,3,7])
    
    ans = upwelling_diffusion_model._band(testa,testb,1.5*testc,0.5*testd,len(testa))
    assert(len(ans) == len(testa))
    assert(ans[0] != 0)
    assert(ans[-1] != 0)
    assert((ans==[6.90789474,-2.90789474,0.73684211]).all)
