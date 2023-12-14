"""
Module to make parameter distributions
"""

import json
import logging

import numpy as np
from scipy.stats import qmc

LOGGER = logging.getLogger(__name__)

# TODO: sensible standard priors for q-terms
prior_flat = {
    "rlamdo": [5, 25],
    "akapa": [0.06, 0.08],
    "cpi": [0.161, 0.569],
    "W": [0.55, 2.55],
    "beto": [0, 7],
    "lambda": [2 / 3.71, 5 / 3.71],
    "mixed": [25, 125],
    "qbmb": [0, 2],
    "qo3": [0.4, 0.6],
    "qdirso2": [-0.55, -0.2],
    "qindso2": [-1.5, -0.5],
    "qbc": [0.1, 0.2],
    "qoc": [-0.1, -0.06],
    "qh2o_ch4": [0.08, 0.1],
    "aerosol_total": [-3.0, -0.5],
    "beta_f": [0.110, 0.465],
}
"""dict: Containing a default prior for parameters """

prior_flat_array = [
    [5, 25],
    [0.06, 0.08],
    [0.161, 0.569],
    [0.55, 2.55],
    [0, 7],
    [2 / 3.71, 5 / 3.71],
    [25, 125],
    [0, 2],
    [0.4, 0.6],
    [-0.55, -0.2],
    [-1.5, -0.5],
    [0.1, 0.2],
    [-0.1, -0.06],
    [0.08, 0.1],
]
"""array: same as prior_flat, but just as an array ordered according to standard ordering"""

ordering_standard = [
    "rlamdo",
    "akapa",
    "cpi",
    "W",
    "beto",
    "lambda",
    "mixed",
    "threstemp",
    "lm",
    "ldtime",
    "qo3",
    "qdirso2",
    "qindso2",
    "qbc",
    "qoc",
    "qbmb",
    "qh2o_ch4",
]
"""list: Containing a default ordering of parameters """

ordering_standard_forc = [
    "rlamdo",
    "akapa",
    "cpi",
    "W",
    "beto",
    "lambda",
    "mixed",
    "threstemp",
    "lm",
    "ldtime",
]
"""list: Containing a default ordering of parameters for forcing run"""
aerosols = ["qdirso2", "qindso2", "qbc", "qoc"]


class _ConfigDistro:
    """
    Class that holds distribution information for parameters and can be used
    to create sample parameter sets from
    """

    def __init__(
        self,
        distro_array=prior_flat_array,
        setvalues={
            "threstemp": 7.0,
            "lm": 40,
            "ldtime": 12,
            "qbmb": 0,
            "qh2o_ch4": 0.091915,
        },
        ordering="standard",
        options={
            "forc": False,
            "method": "latin",
            "aerosol-total": [-0.461590, -0.519163, 0.202893, -0.104478],
        },
    ):  # pylint:disable=dangerous-default-value
        """
        Intialise _ConfigDistro

        Parameters
        ----------
        distro_array: array
              2D array, with first index over the parameters according
              to the ordering given in the ordering parameter. Second index
              of length 2, and gives the ends of the interval for the
              parameter. If the method is latin, these are assumed to be
              extent of the parameter dimension. If the method is gaussian,
              these are assumed to be points one standard deviation away
              from the center of of the distribution for the parameter
              If the distro_array does not have correct type or dimensions
              the default prior will be used.
        setvalues: dict
              dictionary of parameters that will have a set value, keys
              are parameter names, values are their set values
        ordering: list
              list of parameter ordering. This  should match the ordering
              of values in the distro_array
        options: dict
              for optional key word arguments, such as forc, a bool which
              defines whether parameter distribution is meant for forcing
              runs. method, a string to define the configuration sampling
              method. If method='latin', latin hypercube is used, and this
              is the default method. If gaussian is chosen, samples are
              drawn from gaussian distributions for each of the
              parameters. This method will also be cosen if some other random
              string or object is sent for this keyword argument.
              And aerosol_total, which should be a four element array
              defining the proportions between the aerosol forcings for
              dirso2, indso2, bc and oc in that order
        """
        self.options = options
        if "forc" not in options:
            self.options["forc"] = False
        if "method" not in options:
            self.options["method"] = "latin"
        if "aerosol_total" in ordering and "aerosol_total" not in options:
            options["aerosol_total"] = [-0.461590, -0.519163, 0.202893, -0.104478]
        if "aerosol_total" in ordering:
            self.options["aerosol_total"] = np.array(options["aerosol_total"]) / np.sum(
                options["aerosol_total"]
            )
        if ordering == "standard":
            if self.options["forc"]:
                ordering = [
                    o for o in ordering_standard_forc if o not in set(setvalues)
                ]
            else:
                ordering = [o for o in ordering_standard if o not in set(setvalues)]
                if "aerosol_total" in ordering:
                    for aerosol in aerosols:
                        ordering.remove(aerosol)
        else:
            ordering = [o for o in ordering if o not in set(setvalues)]
        self.ordering = ordering
        self.prior = self._set_prior(distro_array)
        self.setvalues = setvalues
        self._set_pamsets_start()

    def _set_prior(self, distro_array):
        """
        Fill out the prior from a given distro_array

        Parameters
        ----------
        distro_array: array
              2D array, with first index over the parameters that according
              to the ordering given in the ordering parameter. Second index
              of length 2, and gives the ends of the interval for the
              parameter. If the method is latin, these are assumed to be
              extent of the parameter dimension. If the method is gaussian
              These are assumed to be points one standard deviation away
              from the center of of the distribution for the parameter
              If the distro_array does not have correct type or dimensions
              the default prior will be used.

        Returns
        -------
             array
             2D array of intervals for all the parameters that don't
             have setvalues according to the ordering.
        """
        prior = np.zeros([len(self.ordering), 2])
        len_given = 0
        try:
            len_given = len(distro_array[:, 0])
            prior[:len_given, :] = distro_array
        except (ValueError, TypeError):
            LOGGER.warning("distro_array not 2 dimensional, disregarding")
        for i, pam in enumerate(self.ordering[len_given:]):
            prior[len_given + i, :] = prior_flat[pam]
        return prior

    def _set_pamsets_start(self):
        """
        Set starting pamsets from the setvalues
        """
        self.pamset_udm_start = {}
        self.pamset_emiconc_start = {}
        for pam, value in self.setvalues.items():
            if pam in ordering_standard_forc:
                self.pamset_udm_start[pam] = value
            else:
                self.pamset_emiconc_start[pam] = value

    def get_samples_from_distro_latin(self, numvalues):
        """
        Get latin hypercube samples over the prior

        Parameters
        ----------
        numvalues: int
           number of values in returned sample

        Returns
        -------
        np.array
             latin hypercube of numvalues samples in the intervals
             given by the prior
        """
        sampler = qmc.LatinHypercube(d=len(self.ordering))
        samples_unit = sampler.random(n=numvalues)
        samples_scaled = np.array(
            qmc.scale(samples_unit, self.prior[:, 0], self.prior[:, 1])
        )
        return samples_scaled

    def get_samples_from_distro_gaussian(self, numvalues):
        """
        Get samples from gaussian distribution over the prior

        Parameters
        ----------
        numvalues: int
           number of values in returned sample

        Returns
        -------
        np.array
           numvalues samples drawn from gaussian distributions
           in the intervals for the parameters given by the prior
        """
        means = np.mean(self.prior, axis=1)
        std = self.prior[:, 1] - means
        samples = np.array(
            [
                np.random.normal(loc=means[i], scale=std[i], size=numvalues)
                for i in range(len(self.prior))
            ]
        ).transpose()
        return samples

    def make_config_list(self, numvalues, json_fname="No", indexer_pre=""):
        """
        Make configuration list from samples over the distribution

        Parameters
        ----------
        numvalues: int
           number of values in returned configuration set
        json_fname: str
           path to file to output configurations
        indexer_pre: str
           prefix to put on the Index for the samples

        Returns
        -------
        list
           list of numvalue configuration dictionaries with parameter
           values drawn according to the distribution objects prior and
           distribution methods, along with set values for the parameters
           that are defined with setvalues
        """
        config_list = [None] * numvalues
        if self.options["method"] == "latin":
            samples = self.get_samples_from_distro_latin(numvalues)
        else:
            samples = self.get_samples_from_distro_gaussian(numvalues)

        for i in range(numvalues):
            pamset_udm = self.pamset_udm_start
            pamset_emiconc = self.pamset_emiconc_start
            for j, pam in enumerate(self.ordering):
                if pam in ordering_standard_forc:
                    pamset_udm[pam] = samples[i, j]
                elif not self.options["forc"]:
                    if pam == "aerosol_total":
                        for anum, aerosol in enumerate(aerosols):
                            pamset_emiconc[aerosol] = (
                                samples[i, j] * self.options["aerosol_total"][anum]
                            )
                    else:
                        pamset_emiconc[pam] = samples[i, j]

            config_list[i] = {
                "pamset_udm": pamset_udm.copy(),
                "pamset_emiconc": pamset_emiconc.copy(),
                "Index": f"{indexer_pre}{i}",
            }
        if json_fname != "No":
            with open(json_fname, "w", encoding="utf-8") as wfile:
                json.dump(config_list, wfile)
        return config_list
