"""
Module for running CICEROSCM in parallel
"""

import logging
from concurrent.futures import ProcessPoolExecutor
from itertools import product

import pandas as pd
import scmdata
from openscm_runner.adapters.ciceroscm_py_adapter.make_scenario_data import (
    SCENARIODATAGETTER,
)
from openscm_runner.adapters.ciceroscm_py_adapter.read_results import CSCMREADER
from openscm_runner.adapters.utils._parallel_process import _parallel_process

from .._utils import check_numeric_pamset
from ..ciceroscm import CICEROSCM

LOGGER = logging.getLogger(__name__)


FRONT_SERIAL = 0
"""int: Number of serial runs to do before starting parallel runs"""

FRONT_PARALLEL = 0
"""int: Number of front parallel runs to do before starting full parallel runs"""


def _execute_run(cfgs, output_variables, scenariodata):
    """
    Run execution method to be used by parallel run

    Parameters
    ----------
    scenariodata : dict
        Dictionary of scenariodicts to run the cscm model for single run

    cfgs : dict
        Dict of configurations with which to run CICEROSCM for single run

    output_vars : list[str]
        Variables to output these follow ScmRun conventions

    Returns
    -------
    :obj:`ScmRun`
        :obj:`ScmRun` instance with all results.

    """
    scenariodata = check_numeric_pamset(
        {"nystart": 1750, "nyend": 2100, "emstart": 1850}, scenariodata
    )
    par_wrapper = CSCMParWrapper(scenariodata)
    try:
        out = par_wrapper.run_over_cfgs(cfgs, output_variables)
    finally:
        LOGGER.info("Finished run")
    return out


def run_ciceroscm_parallel(scendata, cfgs, output_vars, max_workers=4):
    """
    Run CICEROSCM in parallel

    Parameters
    ----------
    scendata : list[dict]
        List of scenariodicts to run the cscm model with
        Scenariodata with which to run

    cfgs : list[dict]
        List of configurations with which to run CICEROSCM
        Should be used for each scendata set

    output_vars : list[str]
        Variables to output, these follow ScmRun conventions

    max_workers : int
        Max number of parallel threads to run. Default value 4.

    Returns
    -------
    :obj:`ScmRun`
        :obj:`ScmRun` instance with all results.
    """
    LOGGER.info("Entered _parallel_ciceroscm")
    if len(scendata) >= max_workers or len(cfgs) * len(scendata) < max_workers:
        runs = [
            {"cfgs": cfgs, "output_variables": output_vars, "scenariodata": scen}
            for scen in scendata
        ]
    else:
        batch_size = len(cfgs) * len(scendata) // max_workers
        runs = [
            {"cfgs": cfgset, "output_variables": output_vars, "scenariodata": scen}
            for scen, cfgset in product(
                scendata,
                [cfgs[i : i + batch_size] for i in range(0, len(cfgs), batch_size)],
            )
        ]
    LOGGER.info("Running in parallel with up to %d workers", max_workers)

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        result = _parallel_process(
            func=_execute_run,
            configuration=runs,
            pool=pool,
            config_are_kwargs=True,
            # no front runs as these defeat the purpose with CICERO-SCM (because
            # it is only parallel on scenarios, not configs)
            front_serial=FRONT_SERIAL,
            front_parallel=FRONT_PARALLEL,
        )

    LOGGER.info("Appending CICERO-SCM results into a single ScmRun")
    result = scmdata.run_append([r for r in result if r is not None])

    return result


class CSCMParWrapper:  # pylint: disable=too-few-public-methods
    """
    CICEROSCM parallel run wrapper
    """

    def __init__(self, scenariodata):
        """
        Intialise CICEROSCM wrapper

        Parameters
        ----------
        scenariodata: dict
              Dictionary of scenariodata to initialise a ciceroscm with
              In addition to nystart and nyend, the scenariodata must contain
              "udir" a path to an utilities directory containing a gaspam file
              "scenname": A scenario name
        """
        self.udir = scenariodata["udir"]
        nystart = scenariodata["nystart"]
        nyend = scenariodata["nyend"]
        self.sdatagetter = SCENARIODATAGETTER(self.udir, nystart, nyend)
        self.resultsreader = CSCMREADER(nystart, nyend)
        self.cscm = CICEROSCM(scenariodata)
        self.scen = scenariodata["scenname"]
        self.model = scenariodata["scenname"]

    def run_over_cfgs(self, cfgs, output_variables):
        """
        Run over each configuration parameter set
        write parameterfiles, run, read results
        and make an ScmRun with results

        Parameters
        ----------
        cfgs : dict
            Dictionary of configurations with which to run CICEROSCM
            This should include dictionaries pamset_udm and a pamset_emiconc
            (though the latter may be empty for a forcing run). In addition
            it should include a string "Index" that identifies this
            configuration set

        output_variables : list[str]
            Variables to output, these follow ScmRun conventions

        Returns
        -------
             :obj:`ScmRun`
             :obj:`ScmRun` instance with all requested output variables
        """
        runs = []
        for pamset in cfgs:
            self.cscm._run(  # pylint: disable=protected-access
                {"results_as_dict": True},
                pamset_udm=pamset["pamset_udm"],
                pamset_emiconc=pamset["pamset_emiconc"],
            )
            for variable in output_variables:
                (
                    years,
                    timeseries,
                    unit,
                ) = self.resultsreader.get_variable_timeseries(
                    self.cscm.results, variable, self.sdatagetter
                )
                if isinstance(years, pd.DataFrame) and years.empty:  # pragma: no cover
                    continue  # pragma: no cover
                runs.append(
                    scmdata.ScmRun(
                        pd.Series(timeseries, index=years),
                        columns={
                            "climate_model": "CICERO-SCM-PY",
                            "model": self.model,
                            "run_id": pamset["Index"],
                            "scenario": self.scen,
                            "region": ["World"],
                            "variable": [variable],
                            "unit": [unit],
                        },
                    )
                )

        return scmdata.run_append(runs)
