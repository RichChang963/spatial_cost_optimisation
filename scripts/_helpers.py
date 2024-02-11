# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

from pathlib import Path
from typing import Union
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import MultiPolygon
from operator import attrgetter
from itertools import takewhile
import matplotlib

def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.
    Note: Must only be called once from the __main__ section of a script.
    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.
    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """

    import logging

    kwargs = snakemake.config.get("logging", dict())
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath(
            "..", "logs", f"{snakemake.rule}.log"
        )
        logfile = snakemake.log.get(
            "python", snakemake.log[0] if snakemake.log else fallback_path
        )
        kwargs.update(
            {
                "handlers": [
                    # Prefer the 'python' log, otherwise take the first log for each
                    # Snakemake rule
                    logging.FileHandler(logfile),
                    logging.StreamHandler(),
                ]
            }
        )
    logging.basicConfig(**kwargs)


def retrieve_snakemake_keys(snakemake):
    return (
        snakemake.input,
        snakemake.config,
        snakemake.wildcards,
        snakemake.log,
        snakemake.output,
    )


def progress_retrieve(url, file, data=None, disable_progress=False, roundto=1.0):
    """
    Function to download data from a url with a progress bar progress in retrieving data
    Parameters
    ----------
    url : str
        Url to download data from
    file : str
        File where to save the output
    data : dict
        Data for the request (default None), when not none Post method is used
    disable_progress : bool
        When true, no progress bar is shown
    roundto : float
        (default 0) Precision used to report the progress
        e.g. 0.1 stands for 88.1, 10 stands for 90, 80
    """
    import urllib

    from tqdm import tqdm

    pbar = tqdm(total=100, disable=disable_progress)

    def dlProgress(count, blockSize, totalSize, roundto=roundto):
        pbar.n = round(count * blockSize * 100 / totalSize / roundto) * roundto
        pbar.refresh()

    if data is not None:
        data = urllib.parse.urlencode(data).encode()

    urllib.request.urlretrieve(url, file, reporthook=dlProgress, data=data)


def mock_snakemake(rulename, **wildcards):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.
    If a rule has wildcards, you have to specify them in **wildcards.
    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import os

    import snakemake as sm
    from packaging.version import Version, parse
    from pypsa.descriptors import Dict
    from snakemake.script import Snakemake

    script_dir = Path(__file__).parent.resolve()
    assert (
        Path.cwd().resolve() == script_dir
    ), f"mock_snakemake has to be run from the repository scripts directory {script_dir}"
    os.chdir(script_dir.parent)
    for p in sm.SNAKEFILE_CHOICES:
        if os.path.exists(p):
            snakefile = p
            break
    kwargs = dict(rerun_triggers=[]) if parse(sm.__version__) > Version("7.7.0") else {}
    workflow = sm.Workflow(snakefile, overwrite_configfiles=[], **kwargs)
    workflow.include(snakefile)
    workflow.global_resources = {}
    rule = workflow.get_rule(rulename)
    dag = sm.dag.DAG(workflow, rules=[rule])
    wc = Dict(wildcards)
    job = sm.jobs.Job(rule, dag, wc)

    def make_accessable(*ios):
        for io in ios:
            for i in range(len(io)):
                io[i] = os.path.abspath(io[i])

    make_accessable(job.input, job.output, job.log)
    snakemake = Snakemake(
        job.input,
        job.output,
        job.params,
        job.wildcards,
        job.threads,
        job.resources,
        job.log,
        job.dag.workflow.config,
        job.rule.name,
        None,
    )
    # create log and output dir if not existent
    for path in list(snakemake.log) + list(snakemake.output):
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    os.chdir(script_dir)
    return snakemake


def __get_country(target, **keys):
    """
    Function to convert country codes using pycountry
    Parameters
    ----------
    target: str
        Desired type of country code.
        Examples:
            - 'alpha_3' for 3-digit
            - 'alpha_2' for 2-digit
            - 'name' for full country name
    keys: dict
        Specification of the country name and reference system.
        Examples:
            - alpha_3="ZAF" for 3-digit
            - alpha_2="ZA" for 2-digit
            - name="South Africa" for full country name
    Returns
    -------
    country code as requested in keys or np.nan, when country code is not recognized
    Example of usage
    -------
    - Convert 2-digit code to 3-digit codes: get_country('alpha_3', alpha_2="ZA")
    - Convert 3-digit code to 2-digit codes: get_country('alpha_2', alpha_3="ZAF")
    - Convert 2-digit code to full name: get_country('name', alpha_2="ZA")
    """
    import pycountry as pyc

    assert len(keys) == 1
    try:
        return getattr(pyc.countries.get(**keys), target)
    except (KeyError, AttributeError):
        return np.nan


def iso2_to_iso3_country(iso2):
    """
    Convert 2-digit to 3-digit country code:
    Parameters
    ----------
    iso2: str
        2-digit country name
    Returns
    ----------
    iso3: str
        3-digit country name
    """
    if iso2 == "SN-GM":
        return f"{iso2_to_iso3_country('SN')}-{iso2_to_iso3_country('GM')}"
    if iso2 == "XK":  # fix for kosovo
        return "XKO"

    iso3 = __get_country("alpha_3", alpha_2=iso2)
    return iso3


def iso3_to_iso2_country(iso3):
    """
    Convert 3-digit to 2-digit country code:
    Parameters
    ----------
    iso3: str
        3-digit country name
    Returns
    ----------
    iso2: str
        2-digit country name
    """
    # if iso3 == "SEN-GMB":
    #     return f"{three_2_two_digits_country('SN')}-{three_2_two_digits_country('GM')}"

    iso2 = __get_country("alpha_2", alpha_3=iso3)

    replace_list = {"GR": "EL", "SW": "SE", "AU": "AT", "GB": "UK"}

    for i in replace_list.items():
        iso2 = iso2.replace(i[0], i[1])

    return iso2


def simplify_polys(polys, minarea=0.1, tolerance=0.01, filterremote=True) -> gpd.GeoDataFrame:
    """function for simplifying polygons.
    Args:
        polys: the geometry column in original GeoDataFrame
        minarea (float64): the minimal available area
        tolerance (float64): tolerance of simplfying polygons
        filterremote (boolean): activate any remote distance from the polygons
    Returns:
        GeoDataFrame: table of GeoDataframe
    """ 

    if isinstance(polys, MultiPolygon):
        polys = sorted(polys.geoms, key=attrgetter("area"), reverse=True)
        mainpoly = polys[0]
        mainlength = np.sqrt(mainpoly.area / (2.0 * np.pi))
        if mainpoly.area > minarea:
            polys = MultiPolygon(
                [
                    p
                    for p in takewhile(lambda p: p.area > minarea, polys)
                    if not filterremote or (mainpoly.distance(p) < mainlength)
                ]
            )
        else:
            polys = mainpoly
    return polys.simplify(tolerance=tolerance)


def plot_style():
    """function for setting up the style of the graphs.
    """ 

    ############        LINES           ############
    matplotlib.rcParams["lines.linewidth"] = 2
    matplotlib.rcParams["lines.linestyle"] = "--"
    ############        FONT           ############
    matplotlib.rcParams["font.family"] = "Flexo"
    matplotlib.rcParams["font.size"] = 14
    ############        AXES           ############
    matplotlib.rcParams["axes.facecolor"] = "#E3E4EA"
    matplotlib.rcParams["axes.edgecolor"] = "white"
    matplotlib.rcParams["axes.labelcolor"] = "black"
    matplotlib.rcParams["axes.labelweight"] = "bold"
    ############        TICKS           ############
    matplotlib.rcParams["xtick.color"] = "#000000"
    matplotlib.rcParams["ytick.color"] = "#000000"
    ############        LEGEND           ############
    matplotlib.rcParams["legend.facecolor"] = "inherit"
    # matplotlib.rcParams['legend.title_fontweight'] = 'bold'
    ############        FIGURE           ############
    matplotlib.rcParams["figure.facecolor"] = "#E3E4EA"
    matplotlib.rcParams["figure.edgecolor"] = "#E3E4EA"
