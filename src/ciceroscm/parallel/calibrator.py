"""
Module to perform calibration of the ciceroscm model
"""
import logging

import numpy as np

from .distributionrun import DistributionRun

LOGGER = logging.getLogger(__name__)


class Calibrator:
    """
    Calibrator class used to do calibration
    """

    def __init__(self, calibdata, config_distro, scendata, subsamplesize=50):
        """
        Initialise calibrator
        Set the calibrators data to calibrate to, prior, scenario
        and subsample size

        Parameters
        ----------
        calibdata : pandas.DataFrame
            Data Frame containing datapoints to calibrate to
            For each datapoint the following columns should be filled
            in: 'Variable Name', which should be the name of the variable
            to be fitted, this should be an ScmRun variable. 'Yearstart_norm'
            'Yearend_norm', 'Yearstart_change' and 'Yearend_change'.
            These assume that the variable has an expected change between
            a norm period and a change period, and define the inclusive
            ranges of these extents. If the variable change is just between
            single years, set Yearstart_norm and Yearend_norm equal.
            'Central value' is the observed change in the variable and
            finally 'sigma' is the standard deviation in the observed
            value
        config_distro: :obj: ConfigDistro
            An object of the ConfigDistro object which has defined priors
            an set values, and can produce draws for the calibrator
            to pick from
        scendata: dict
            Dictionary that defines what to run in order to calibrate
        subsamplesize: The calibrator will try to get fitted values by
            Sampling in chunks and doing an accept/reject on the samples
            after running them. This will be done recursively untill the
            maximal recursion depth is reached, or the requested number of
            accepted samples is found. The subsamplesize decides the number
            of configuration in each chunk.
        """
        self.calibdata = calibdata
        self.subsamplesize = subsamplesize
        self.config_distro = config_distro
        self.scendata = scendata
        self.recurse_max = 1

    def assign_log_liklihood(self, res):
        """
        Assign a log likelihood for a result

        From the calibdata determine a log likelihood for
        a results run

        Parameters
        ----------
        res: ScmRun
            A filtered ScmRun for a single configuration run

        Returns
        -------
            float
            The log likelihood for the result matching the
            calibration data
        """
        logli = 0
        for (
            index,  # pylint: disable=unused-variable
            datapoint,
        ) in self.calibdata.iterrows():
            vres = (
                res.filter(
                    variable=datapoint["Variable Name"], year=datapoint["Yearend"]
                ).values[0]
                - res.filter(
                    variable=datapoint["Variable Name"], year=datapoint["Yearstart"]
                ).values[0]
            )
            sigma = datapoint["sigma"]
            vexp = datapoint["Central Value"]
            logli = logli - 0.5 * (vres - vexp) ** 2 / sigma**2  # + np.log(sigma)
        print(np.exp(logli))
        return logli

    def find_distance(self, res):
        """
        Find the probability scaled distance between result
        and calibration data

        Using the distance measure Sum (x-mu)**2/sigma**2
        where the sum is over each of the datapoints in the calibration
        data

        Parameters
        ----------
        res: ScmRun
            A filtered ScmRun for a single configuration run

        Returns
        -------
            float
            The probability distance between the result and the
            calibration data
        """
        distance = 0
        for (
            index,  # pylint: disable=unused-variable
            datapoint,
        ) in self.calibdata.iterrows():

            if datapoint["Yearstart_change"] == datapoint["Yearend_change"]:
                vres = (
                    res.filter(
                        variable=datapoint["Variable Name"],
                        year=datapoint["Yearstart_change"],
                    ).values[0]
                    - res.filter(
                        variable=datapoint["Variable Name"],
                        year=datapoint["Yearstart_norm"],
                    ).values[0]
                )
            else:
                vres = np.mean(
                    res.filter(
                        variable=datapoint["Variable Name"],
                        year=range(
                            datapoint["Yearstart_change"], datapoint["Yearend_change"]
                        ),
                    ).values
                ) - np.mean(
                    res.filter(
                        variable=datapoint["Variable Name"],
                        year=range(
                            datapoint["Yearstart_norm"], datapoint["Yearend_norm"]
                        ),
                    ).values
                )
            sigma = datapoint["sigma"]
            vexp = datapoint["Central Value"]
            print(f"vres: {vres}")
            print(f"vexp: {vexp}")
            print(f"sigma: {sigma}")
            distance = distance + (vres - vexp) ** 2 / sigma**2
        return distance

    def reject_samples(self, results):
        """
        Reject and accept samples probabilistically

        Sample distances to observations are compared
        to the distance between two randomally drawn
        samples from a normal distribution multiplied
        with the datapoint dimensionality. If the distance
        is less than the random draw, the sample is accepted
        otherwise the sample is rejected

        Parameters
        ----------
        results: ScmRun
            ScmRun for a chunk of model configurations

        Returns
        -------
            list
            A list of the configuration indices that passed
            the test
        """
        indices = results.get_unique_meta("run_id")
        print(indices)

        draws = np.square(
            np.random.normal(size=(len(indices)))
            - np.random.normal(size=(len(indices)))
        )
        keep = []
        for i, index in enumerate(indices):
            distance = self.find_distance(results.filter(run_id=index))
            print(f"distance: {distance}")
            print(f"draw: {draws[i]}")
            if draws[i] * len(self.calibdata) > distance:
                keep.append(index)
        return keep

    def get_n_samples(self, numvalues, current_samples, kept_configs, recurse_num):
        """
        Recurse untill numvalues samples are found

        This method creates a new distribution of the object samplesize
        and runs over it if the required number of values to the
        distribution is reached after this run, it returns,
        otherwise it calls it self and increases the recursion depth

        Parameters
        ----------
        numvalues: int
            Number of accepted values requested
        current_samples: list
            Index list of samples kept
        kept_configs: list
            List of configuration dictionaries accepted
        recurse_num: int
            Current recursion depth. Defines when to stop
            recursion to avoid infinite recursion

        Returns
        -------
            list
            list of accepted configuration dictionaries
        """
        distro_run = DistributionRun(
            self.config_distro,
            indexer_pre=f"{recurse_num}_",
            numvalues=self.subsamplesize,
        )
        output_vars = self.calibdata["Variable Name"].values
        results = distro_run.run_over_distribution(self.scendata, output_vars)
        kept = self.reject_samples(results)
        current_samples.extend(kept)
        for item in distro_run.cfgs:
            if item["Index"] in kept:
                kept_configs.append(item)
        if len(current_samples) >= numvalues:
            return kept_configs
        if recurse_num > self.recurse_max:
            LOGGER.warning("Reached recurse limit with only %d", len(current_samples))
            return kept_configs
        return self.get_n_samples(
            numvalues, current_samples, kept_configs, recurse_num + 1
        )
