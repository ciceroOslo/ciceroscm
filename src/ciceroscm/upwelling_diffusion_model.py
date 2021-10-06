import numpy as np


def _coefic(s, t, p):
    """
    Calculate denisty coefficients
    """
    coefic = (
        19652.21
        + 148.4206 * t
        - 2.327105 * t ** 2
        + 1.360477e-2 * t ** 3
        - 5.155288e-5 * t ** 4
        + 3.239908 * p
        + 1.43713e-3 * t * p
        + 1.16092e-4 * t ** 2 * p
        - 5.77905e-7 * t ** 3 * p
        + 8.50935e-5 * p ** 2
        - 6.12293e-6 * t * p ** 2
        + 5.2787e-8 * t ** 2 * p ** 2
        + 54.6746 * s
        - 0.603459 * t * s
        + 1.09987e-2 * t ** 2 * s
        - 6.1670e-5 * t ** 3 * s
        + 7.944e-2 * s ** 1.5
        + 1.6483e-2 * t * s ** 1.5
        - 5.3009e-4 * t ** 2 * s ** 1.5
        + 2.2838e-3 * p * s
        - 1.0981e-5 * t * p * s
        - 1.6078e-6 * t ** 2 * p * s
        + 1.91075e-4 * p * s ** 1.5
        - 9.9348e-7 * p ** 2 * s
        + 2.0816e-8 * t * p ** 2 * s
        + 9.1697e-10 * t ** 2 * p ** 2 * s
    )
    return coefic


def _denso(s, t):
    """
    Calculate density at p0=0
    """
    denso = (
        999.842594
        + 6.793952e-2 * t
        - 9.095290e-3 * t ** 2
        + 1.001685e-4 * t ** 3
        - 1.120083e-6 * t ** 4
        + 6.536332e-9 * t ** 5
        + 8.24493e-1 * s
        - 4.0899e-3 * t * s
        + 7.6438e-5 * t ** 2 * s
        - 8.2467e-7 * t ** 3 * s
        + 5.3875e-9 * t ** 4 * s
        - 5.72466e-3 * s ** 1.5
        + 1.0227e-4 * t * s ** 1.5
        - 1.6546e-6 * t ** 2 * s ** 1.5
        + 4.8314e-4 * s ** 2
    )
    return denso


def _density(p0, t0):
    """
    Calculate water denisity from equation of state
    """
    s = 35.0
    return _denso(s, t0) / (1.0 - p0 / _coefic(s, t0, p0))


def _band(a, b, c, d, lm, ans):
    """
    Calculate band
    """
    alfa = np.zeros(lm - 1)
    ans = np.zeros(lm)
    bbeta = np.zeros(lm)

    alfa[0] = -b[1] / a[1]
    bbeta[0] = d[1] / a[1]

    for i in range(1, lm - 1):
        tem = a[i] * alfa[i - 1] + b[i]
        alfa[i] = -c[i] / tem
        bbeta = (d[i] - a[i] * bbeta[i - 1]) / tem
    tem = a[i] * alfa[i - 1] + b[i]
    ans[lm - 1] = (d[i] - a[i] * bbeta[i - 1]) / tem

    for i in range(1, lm - 1):
        j = lm - 1 - i
        ans[j] = alfa[j] * ans[j + 1] + bbeta[j]

    return ans


class UpwellingDiffusionModel:
    """
    Class to handle energy budget upwelling and downwelling
    """

    def __init__(self, params):
        """
        Intialise
        """
        self.foan = 0.61  # make changable
        self.foas = 0.81  # make changable
        self.day = 86400.0
        self.year = 365.0
