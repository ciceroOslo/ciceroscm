"""
Code to support running in parallel
"""

import logging
import time
from concurrent.futures import as_completed

from tqdm.autonotebook import tqdm

LOGGER = logging.getLogger(__name__)

# TODO: handle configuration in a consistent manner
_default_tqdm_params = {
    "mininterval": 5,
    "unit": "it",
    "unit_scale": True,
}


def progress(*args, **kwargs):
    """
    Progress bar

    Uses ``tqdm.autonotebook`` to automatically use a native Jupyter widget
    when executing within a Jupyer Notebook.

    Parameters
    ----------
    *args
        Passed to the tqdm

    **kwargs
        Passed to the tqdm

    Returns
    -------
    tqdm.auto_notebook.tqdm
        tqdm instance with consistent configuration
    """
    kwargs = {**_default_tqdm_params, **kwargs}
    return tqdm(*args, **kwargs)


def _run_serial(func, configs, config_are_kwargs, desc):
    LOGGER.debug("Entering _run_serial")

    if config_are_kwargs:
        LOGGER.debug("Treating config as kwargs")
        res = [func(**a) for a in progress(configs, desc=desc)]
    else:
        LOGGER.debug("Treating config as args")
        res = [func(a) for a in progress(configs, desc=desc)]

    LOGGER.debug("Exiting _run_serial")
    return res


def _run_parallel(  # pylint:disable=too-many-arguments,too-many-positional-arguments
    pool, timeout, func, configs, config_are_kwargs, desc, bar_start
):
    LOGGER.debug("Entering _run_parallel")

    if config_are_kwargs:
        LOGGER.debug("Treating config as kwargs")
        futures = [pool.submit(func, **a) for a in configs]
    else:
        LOGGER.debug("Treating config as args")
        futures = [pool.submit(func, a) for a in configs]

    LOGGER.debug("Waiting for jobs to complete")
    for i, future in progress(
        enumerate(as_completed(futures, timeout=timeout)), total=len(futures), desc=desc
    ):
        if future.exception() is not None:
            time.sleep(2)  # let buffer flush out
            print(
                "One of the processes failed, see error below (was something "
                "unable to be pickled?)"
            )
            raise future.exception()

        LOGGER.debug("Job %s completed", i + bar_start)

    res = []

    LOGGER.debug("Collecting results")
    for i, future in enumerate(futures):
        try:
            res.append(future.result())
            LOGGER.debug("Retrived result %s", i + bar_start)
        except Exception as exc:  # pylint:disable=broad-except
            LOGGER.debug("Retrieving result %s failed", i + bar_start)
            res.append(exc)

    LOGGER.debug("Exiting _run_parallel")
    return res


def _parallel_process(  # pylint:disable=too-many-arguments,too-many-positional-arguments
    func,
    configuration,
    pool=None,
    config_are_kwargs=False,
    front_serial=3,
    front_parallel=2,
    timeout=None,
):
    """
    Run a process in parallel with a progress bar.

    Borrowed from: https://github.com/openscm/openscm-runner
    And there adapted from http://danshiebler.com/2016-09-14-parallel-progress-bar/

    Parameters
    ----------
    func : function
        A function to apply to each set of arguments in ``configuration``

    configuration : sequence
        An array of configuration with which to run ``func``.

    pool : :obj:`concurrent.futures.ProcessPoolExecutor`
        Pool in which to execute the jobs. If ``None``, the jobs will be executed
        serially in a single process (useful for debugging and benchmarking).

    config_are_kwargs : bool
        Are the elements of ``configuration`` intended to be used as keyword arguments
        when calling ``func``.

    front_serial : int
        The number of iterations to run serially before kicking off the parallel job.
        Useful for debugging.

    front_parallel : int
        The number of initial iterations to run parallel before kicking off the rest
        of the parallel jobs. Useful for debugging (especially if pickling is
        possible).

    timeout : float
        How long to wait for processes to complete before timing out. If
        ``None``, there is no timeout limit.

    Returns
    -------
    sequence
        Results of calling ``func`` with each configuration in ``configuration``
    """
    front_serial_res = []
    if front_serial > 0:
        LOGGER.debug("Running front serial jobs")
        front_serial_res = _run_serial(
            func=func,
            configs=configuration[:front_serial],
            config_are_kwargs=config_are_kwargs,
            desc="Front serial",
        )

    if pool is None:
        LOGGER.info("No pool provided, running rest of the jobs serially and returning")
        rest = _run_serial(
            func=func,
            configs=configuration[front_serial:],
            config_are_kwargs=config_are_kwargs,
            desc="Serial runs",
        )

        return rest + front_serial_res

    front_parallel_res = []
    if front_parallel > 0:
        LOGGER.debug("Running front parallel jobs")
        front_parallel_res = _run_parallel(
            pool=pool,
            timeout=timeout,
            func=func,
            configs=configuration[front_serial : front_serial + front_parallel],
            config_are_kwargs=config_are_kwargs,
            desc="Front parallel",
            bar_start=front_serial,
        )

    LOGGER.debug("Running rest of parallel jobs")
    rest = _run_parallel(
        pool=pool,
        timeout=timeout,
        func=func,
        configs=configuration[front_serial + front_parallel :],
        config_are_kwargs=config_are_kwargs,
        desc="Parallel runs",
        bar_start=front_serial + front_parallel,
    )

    return front_serial_res + front_parallel_res + rest
