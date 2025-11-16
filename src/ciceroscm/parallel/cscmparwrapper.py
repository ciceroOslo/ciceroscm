"""
Module for running CICEROSCM in parallel
"""

import logging
from concurrent.futures import ProcessPoolExecutor
from itertools import product

import pandas as pd

from .._utils import check_numeric_pamset
from ..ciceroscm import CICEROSCM
from ..formattingtools.reformat_cscm_results import CSCMREADER
from ..formattingtools.reformat_inputdata_to_cscm_format import COMMONSFILEWRITER
from ._parallel_process import _parallel_process

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

    result = pd.concat([r for r in result if r is not None])

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
              the data needed for a CICEROSCM instance including
              either "gaspam_file"  a path to a gaspam file or "gaspam_data"
              and "scenname": A scenario name
        """
        nystart = scenariodata["nystart"]
        nyend = scenariodata["nyend"]
        if "gaspam_file" in scenariodata:
            self.sdatagetter = COMMONSFILEWRITER(
                scenariodata["gaspam_file"], nystart, nyend
            )
        else:
            self.sdatagetter = COMMONSFILEWRITER(
                scenariodata["gaspam_data"], nystart, nyend
            )
        self.resultsreader = CSCMREADER(nystart, nyend)
        self.cscm = CICEROSCM(scenariodata)
        self.scen = scenariodata["scenname"]
        self.model = scenariodata["scenname"]

    def run_over_cfgs(self, cfgs, output_variables, carbon_cycle_outputs=True):
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

        carbon_cycle_outputs : bool
            If carbon cycle outputs should be included in results.
            Default is True, but can be muted if user wants

        Returns
        -------
            pd.DataFrame
            A pd.DataFrame instance with all requested output variables
            over the run timeseries
        """
        runs = []
        for pamset in cfgs:
            self.cscm._run(  # pylint: disable=protected-access
                {"results_as_dict": True, "carbon_cycle_outputs": carbon_cycle_outputs},
                pamset_udm=pamset.get("pamset_udm", None),
                pamset_emiconc=pamset.get("pamset_emiconc", None),
                pamset_carbon=pamset.get("pamset_carbon", None),
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
                # TODO: change CICERO-SCM-PY to include version number and/or
                # thermal model and carbon cycle used from the cscm instance
                data = [
                    "CICERO-SCM-PY",
                    self.model,
                    pamset["Index"],
                    self.scen,
                    "World",
                    variable,
                    unit,
                ]
                data.extend(timeseries)
                index = [
                    "climate_model",
                    "model",
                    "run_id",
                    "scenario",
                    "region",
                    "variable",
                    "unit",
                ]
                index.extend(years)
                runs.append(pd.Series(data=data, index=index))
        return pd.concat(runs, axis=1).transpose()
