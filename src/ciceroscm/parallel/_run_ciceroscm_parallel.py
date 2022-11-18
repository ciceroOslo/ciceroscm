"""
Module for running CICEROSCM in parallel
"""
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from itertoools import product
import scmdata

from ....settings import config
from ...utils._parallel_process import _parallel_process

LOGGER = logging.getLogger(__name__)


FRONT_SERIAL = 0
"""int: Number of serial runs to do before starting parallel runs"""

FRONT_PARALLEL = 0
"""int: Number of front parallel runs to do before starting full parallel runs"""


def run_ciceroscm_parallel(scendata, cfgs, output_vars, _execute_run):
    """
    Run CICEROSCM in parallel

    Parameters
    ----------
    scendata : List of scenariodicts to run the cscm model with
        Scenariodata with which to run

    cfgs : list[dict]
        List of configurations with which to run CICEROSCM
        Should be used for each scendata set

    output_vars : list[str]
        Variables to output (may require some fiddling with ``out_x``
        variables in ``cfgs`` to get this right)

    Returns
    -------
    :obj:`ScmRun`
        :obj:`ScmRun` instance with all results.
    """
    LOGGER.info("Entered _parallel_ciceroscm")
    max_workers = int(config.get("CICEROSCM_WORKER_NUMBER", os.cpu_count()))
    if len(scendata) >= max_workers or len(cfgs)*len(scendata) < max_workers:
        runs = [
            {"cfgs": cfgs, "output_variables": output_vars, "scenariodata": scen}
            for scen in scendata
        ]
    else:
        batch_size = len(cfgs)*len(scendata)//max_workers
        runs = [
            {"cfgs": cfgset, "output_variables": output_vars, "scenariodata": scen}
            for scen,cfgset in product(scendata,[cfgs[i:i+batch_size] for i in range(0, len(cfgs), batch_size)])
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

