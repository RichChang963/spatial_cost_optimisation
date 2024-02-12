# SPDX-FileCopyrightText: : 2024 The spatial_cost_optimisation Author
#
# SPDX-License-Identifier: MIT

import logging
import time
import os
os.environ["USE_PYGEOS"] = "0"
os.environ["OGR_GEOMETRY_ACCEPT_UNCLOSED_RING"] = "NO"

import atlite
import geopandas as gpd
import pandas as pd
from dask.distributed import Client, LocalCluster

logger = logging.getLogger(__name__)

from _helpers import configure_logging


def create_cutout(onshore_shapes:str, offshore_shapes:str, weather_year:int, 
cutout_params: dict, output_path:str):
    """Create cutout of the node boundary to prepare for ERA5 API call.  

    Parameters
    ----------
    onshore_shapes : str
        file path of the onshore shape GeoJSON from build_shape.
    offshore_shapes : str
        file path of the offshore shape GeoJSON from build_shape.
    weather_year : int
        selected weather year (integer) from config.yaml.
    cutout_params : dict
        dictionary of parameters for cutout from config.yaml.
    output_path : str
        output file path to save cutout.
    """
    logger.info(f"Extract coutout features of {str(weather_year)} from ERA5...")
    start = time.time()
    snapshots = (pd.date_range(
        f"{str(weather_year)}-01-01", 
        f"{str(int(weather_year)+1)}-01-01", 
        freq="H", inclusive="left"))
    period = [snapshots[0], snapshots[-1]]
    # set time variable in cutout_params: either slice the time frame if already defined
    # , otherwise return with period variable.
    cutout_params["time"] = slice(*cutout_params.get("time", period))
    # If one of the parameters is there
    if {"x", "y", "bounds"}.isdisjoint(cutout_params):
        # Determine the bounds from bus regions with a buffer of two grid cells
        onshore = gpd.read_file(onshore_shapes)
        offshore = gpd.read_file(offshore_shapes)
        regions = pd.concat([onshore, offshore])
        d = max(cutout_params.get("dx", 0.25), cutout_params.get("dy", 0.25)) * 2
        cutout_params["bounds"] = regions.total_bounds + [-d, -d, d, d]
    elif {"x", "y"}.issubset(cutout_params):
        cutout_params["x"] = slice(*cutout_params["x"])
        cutout_params["y"] = slice(*cutout_params["y"])

    logging.info(f"Preparing cutout with parameters {cutout_params}.")
    
    if snakemake.config["field"] == "energy":
        features = ["wind", "height", "influx", "temperature", "runoff"]
    elif snakemake.config["field"] == "agriculture":
        features = ["wind", "height", "influx", "temperature", "runoff"]
    else:
        features = ["height"]

    cutout = atlite.Cutout(output_path, **cutout_params)
    cutout.prepare(features=features)
        
    end = time.time()
    logger.info(f"Cutout features extracted. Processing Time: {round((end - start), 2)}s")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("02_build_cutout", weather_year=2012)

    configure_logging(snakemake)
    nprocesses = int(snakemake.threads)
    cluster = LocalCluster(n_workers=nprocesses, threads_per_worker=1)
    client = Client(cluster, asynchronous=True)

    cutout_params = snakemake.config["atlite"]["cutouts"]
    create_cutout(snakemake.input.onshore_shapes, 
                  snakemake.input.offshore_shapes, 
                  snakemake.wildcards.weather_year,
                  cutout_params,
                  snakemake.output.cutout)
