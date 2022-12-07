"""
Module to simplify doing a run over a configuration distribution
"""

import json
import logging
import os

from .cscmparwrapper import run_ciceroscm_parallel

LOGGER = logging.getLogger(__name__)


class DistributionRun:
    """
    Class that defines the distribution run
    """

    def __init__(
        self, distro_config, json_file_name="no_file", indexer_pre="", numvalues=100
    ):
        """
        Intialise DistributionRun

        Parameters
        ----------
        distro_config: obj
              Instance of the _ConfigDistro class that defines the
              distribution to be run, in this case the configurations
              will be generated from numvalues of this distribution
              If the file json_file_name exists, this will be ignored
        json_file_name: str
              path to json_file with configurations if the file
              exists, configurations will be read from here
        indexer_pre: str
              prefix to put on the Index for the samples
        numvalues: int
              number of configurations to generate, ignored if
              samples_from_json = True
        """
        if os.path.exists(json_file_name):
            with open(json_file_name, "r", encoding="utf-8") as rfile:
                self.cfgs = json.load(rfile)
        else:
            self.cfgs = distro_config.make_config_list(
                numvalues, indexer_pre=indexer_pre
            )

    def run_over_distribution(self, scendata, output_vars, max_workers=4):
        """
        Run over distribution

        scendata: list or dict
            list of dictionaries with the scenarios to run over
            or dictionary of single scenario
        output_vars: list
            list of variables to output
        max_workers: int
            maximum number of parallel node workers to employ
            for the run
        """
        if isinstance(scendata, dict):
            scendata = [scendata]
        results = run_ciceroscm_parallel(scendata, self.cfgs, output_vars, max_workers)
        return results

    def write_configs_to_json(self, json_file_name):
        """
        Write configs to json

        Parameters
        ----------
        json_file_name: str
            path for json file to write configurations in
        """
        with open(json_file_name, "w", encoding="utf-8") as wfile:
            json.dump(self.cfgs, wfile)
