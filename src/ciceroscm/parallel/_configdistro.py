"""
Module to make parameter distributions
"""

import json
import logging

import numpy as np
from scipy.stats import qmc

from ..carbon_cycle.carbon_cycle_mod import CARBON_CYCLE_MODEL_REQUIRED_PAMSET

LOGGER = logging.getLogger(__name__)

# TODO: sensible standard priors for q-terms
prior_flat = {
    "rlamdo": [5, 25],
    "akapa": [0.06, 0.8],
    "cpi": [0.161, 0.569],
    "W": [0.55, 2.55],
    "beto": [0, 7],
    "lambda": [2 / 3.71, 5 / 3.71],
    "mixed": [25, 125],
    "ocean_efficacy": [0.8, 1.2],
    "qbmb": [0, 2],
    "qo3": [0.4, 0.6],
    "qdirso2": [-0.006, -0.001],
    "qindso2": [-0.03, -0.01],
    "qbc": [0.004, 0.05],
    "qoc": [-0.008, -0.001],
    "qh2o_ch4": [0.08, 0.1],
    "beta_f": [0.110, 0.465],
    "mixed_carbon": [25, 125],
    "ml_w_sigmoid": [2.0, 4.0],
    "ml_fracmax": [0.3, 0.8],
    "ml_t_half": [-0.5, 1.0],
    "npp0": [55, 75],
    "t_half": [0.3, 0.8],
    "w_sigmoid": [5, 10],
    "t_threshold": [3, 7],
    "w_threshold": [5, 10],
    "solubility_sens": [0.01, 0.03],
    "solubility_limit": [0.4, 0.8],
    # Flat carbon cycle decay function parameters
    "rb_coef0": [0.3, 0.7],
    "rb_coef1": [0.1, 0.4],
    "rb_coef2": [0.1, 0.4],
    "rb_tim0": [1.0, 5.0],
    "rb_tim1": [5.0, 20.0],
    "rb_tim2": [40.0, 200.0],
    "rs_coef0": [0.05, 0.15],
    "rs_coef1": [0.5, 0.7],
    "rs_coef2": [0.1, 0.2],
    "rs_coef3": [0.1, 0.2],
    "rs_tim0": [0.5, 3],
    "rs_tim1": [5.0, 20.0],
    "rs_tim2": [40.0, 200.0],
}
"""dict: Containing a default prior for parameters """

ordering_standard_forc = [
    "rlamdo",
    "akapa",
    "cpi",
    "W",
    "beto",
    "lambda",
    "mixed",
    "ocean_efficacy",
    "threstemp",
    "lm",
    "ldtime",
]
"""list: Containing a default ordering of parameters for forcing run"""

set_values_default = {
    "threstemp": 7.0,
    "lm": 40,
    "ldtime": 12,
    "qbmb": 0,
    "qh2o_ch4": 0.091915,
}
"""dict: Containing default set values for parameters"""

options_default = {"forc": False, "method": "latin"}
"""dict: Containing default options for making parameter distributions"""


class _ConfigDistro:
    """
    Class that holds distribution information for parameters and can be used
    to create sample parameter sets from
    """

    def __init__(self, distro_dict=None, setvalues=None, options=None):
        """
        Intialise _ConfigDistro

        Parameters
        ----------
        distro_array: dict
              Listing parameters to make a config distro over. Keys should be
              parameters to make the distribution over and keys should list
              two numerical values that describes the configuration over which
              to make the configurations. If the method is latin, these are assumed to be
              extent of the parameter dimension. If the method is gaussian,
              these are assumed to be points one standard deviation away
              from the center of of the distribution for the parameter
              If the distro_array does not have correct type or dimensions
              the default prior will be used.
        setvalues: dict
              dictionary of parameters that will have a set value, keys
              are parameter names, values are their set values
        options: dict
              for optional key word arguments, such as forc, a bool which
              defines whether parameter distribution is meant for forcing
              runs. method, a string to define the configuration sampling
              method. If method='latin', latin hypercube is used, and this
              is the default method. If gaussian is chosen, samples are
              drawn from gaussian distributions for each of the
              parameters. This method will also be cosen if some other random
              string or object is sent for this keyword argument.
        """
        if options is None:
            options = options_default.copy()
        if setvalues is None:
            setvalues = set_values_default.copy()
        if distro_dict is None:
            distro_dict = prior_flat.copy()
        self.options = options
        if "forc" not in options:
            self.options["forc"] = False
        if "method" not in options:
            self.options["method"] = "latin"
        if self.options["forc"]:
            ordering = [
                o
                for o in sorted(distro_dict.keys())
                if (o in ordering_standard_forc) and (o not in set(setvalues))
            ]
        else:
            ordering = [
                o for o in sorted(distro_dict.keys()) if o not in set(setvalues)
            ]

        self.ordering = ordering
        self.prior = self._set_prior(distro_dict)

        self.setvalues = setvalues
        self._set_pamsets_start()

    def _set_prior(self, distro_dict):
        """
        Fill out the prior from a given distro_array

        Parameters
        ----------
        distro_dict: dict
              Dictionary with parameters over which to make config distro.
              Keys are parameter names and values should be numerical arrays
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
        for i, pam in enumerate(self.ordering):
            prior[i, :] = distro_dict[pam]
        return prior

    def _set_pamsets_start(self):
        """
        Set starting pamsets from the setvalues
        """
        self.pamset_udm_start = {}
        self.pamset_emiconc_start = {}
        self.pamset_carbon_start = {}
        for pam, value in self.setvalues.items():
            if pam in ordering_standard_forc:
                self.pamset_udm_start[pam] = value
            elif pam in CARBON_CYCLE_MODEL_REQUIRED_PAMSET:
                self.pamset_carbon_start[pam] = value
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
        print(self.prior[:, 0])
        print(self.prior[:, 1])
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

    def make_config_lists(
        self,
        numvalues,
        json_fname=None,
        indexer_pre="",
        max_chunk_size=None,
        meta_info=None,
    ):  # pylint:disable=too-many-arguments,too-many-positional-arguments
        """
        Make configuration list or chunked lists from samples over the distribution

        Parameters
        ----------
        numvalues: int
           number of values in returned configuration set
        json_fname: str
           path to file to output configurations
        indexer_pre: str
           prefix to put on the Index for the samples
        max_chunk_size : int
            Maximum number of samples in a chunk, use this
        to split the full distribution into more managable
        chunks of maximum this size
        meta_info: dict, optional
            dictionary with meta optional information to add to the json file

        Returns
        -------
        list
           list of numvalue configuration dictionaries with parameter
           values drawn according to the distribution objects prior and
           distribution methods, along with set values for the parameters
           that are defined with setvalues
        """
        if self.options["method"] == "latin":
            samples = self.get_samples_from_distro_latin(numvalues)
        else:
            samples = self.get_samples_from_distro_gaussian(numvalues)
        if max_chunk_size is None or max_chunk_size > numvalues:
            return self.make_single_config_list(
                samples,
                numvalues,
                json_fname=json_fname,
                indexer_pre=indexer_pre,
                meta_info=meta_info,
            )
        config_nums = np.ceil(numvalues / max_chunk_size).astype(int)
        config_chunk_list = [None] * config_nums
        for config_num in range(config_nums):
            num_this_chunk = np.min(
                (max_chunk_size, numvalues - config_num * max_chunk_size)
            )
            start_num = config_num * max_chunk_size
            json_fname_now = None
            if json_fname is not None:
                json_fname_now = json_fname.replace(".json", f"chunk_{config_num}.json")
            config_chunk_list[config_num] = self.make_single_config_list(
                samples[start_num : start_num + num_this_chunk],
                numvalues=num_this_chunk,
                json_fname=json_fname_now,
                indexer_pre=f"{indexer_pre}_{config_num}_",
                meta_info=meta_info,
            )
        return config_chunk_list

    def make_single_config_list(
        self, samples, numvalues, json_fname=None, indexer_pre="", meta_info=None
    ):  # pylint:disable=too-many-arguments,too-many-positional-arguments
        """
        Make configuration list from samples

        Parameters
        ----------
        samples : list
            list of samples from the distribution to make into list
        numvalues: int
           number of values in returned configuration set
        json_fname: str
           path to file to output configurations
        indexer_pre: str
           prefix to put on the Index for the samples
        meta_info: dict, optional
            dictionary with meta optional information to add to the json file

        Returns
        -------
        list
           list of numvalue configuration dictionaries with parameter
           values drawn according to the distribution objects prior and
           distribution methods, along with set values for the parameters
           that are defined with setvalues
        """
        config_list = [None] * numvalues
        for i in range(numvalues):
            pamset_udm = self.pamset_udm_start.copy()
            pamset_emiconc = self.pamset_emiconc_start.copy()
            pamset_carbon = self.pamset_carbon_start.copy()
            for j, pam in enumerate(self.ordering):
                if pam in ordering_standard_forc:
                    pamset_udm[pam] = samples[i, j]
                elif (
                    not self.options["forc"]
                    and pam in CARBON_CYCLE_MODEL_REQUIRED_PAMSET
                ):
                    pamset_carbon[pam] = samples[i, j]
                elif not self.options["forc"] and pam.startswith(("rb_", "rs_")):
                    pamset_carbon[pam] = samples[i, j]
                elif not self.options["forc"]:
                    pamset_emiconc[pam] = samples[i, j]

            config_list[i] = {
                "pamset_udm": pamset_udm.copy(),
                "pamset_emiconc": pamset_emiconc.copy(),
                "pamset_carbon": pamset_carbon.copy(),
                "Index": f"{indexer_pre}{i}",
            }
        if json_fname is not None:
            if meta_info is None:
                with open(json_fname, "w", encoding="utf-8") as wfile:
                    json.dump(config_list, wfile)
            else:
                with open(json_fname, "w", encoding="utf-8") as wfile:
                    json.dump(
                        {"configurations": config_list, "meta_info": meta_info}, wfile
                    )
        return config_list
