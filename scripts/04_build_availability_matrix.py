import logging
import time
import os
os.environ["USE_PYGEOS"] = "0"
os.environ["OGR_GEOMETRY_ACCEPT_UNCLOSED_RING"] = "NO"
import progressbar as pgb
from dask.distributed import Client, LocalCluster
import functools
import collections # fix atlite iterable issue. Import it before importing atlite package
collections.Iterable = collections.abc.Iterable # fix atlite iterable issue.
import atlite
from atlite.gis import ExclusionContainer
import geopandas as gpd
import xarray as xr
import numpy as np

from _helpers import configure_logging


logger = logging.getLogger(__name__)



def calculate_land_availability_matrix(country_shapes, node, natura, gebco, copernicus_map, 
    cutout_name, cutout_params, resolution, output_path, kwargs) -> xr.DataArray:
    """Create a land availability matrix.

    Parameters
    ----------
    country_shapes : str
        file path of shapfile of a whole country
    node : shp
        shapfile of node division of a country under nuts2 level 
    natura: str
        file path of the natura map
    gebco : str
        file path of the underwater depth map
    copernicus_map: str
        file path of copernicus map
    cutout_name : str
        path of the cutout file
    cutout_params: dict
        a dictionary of cutout variables for era5 api
    resolution : int
        resolution of the matrix in meters
    output_path : str
        output path to save DataArray as NetCDF file
    kwargs : dict
        a dictionary including numbers of threads for multi-threading process
    """
    COPERNICUS_CRS = "EPSG:4326"
    GEBCO_CRS = "EPSG:4326"
    start = time.time()
    node_gdf = gpd.read_file(node).set_index("node_name")

    cutout = atlite.Cutout(cutout_name, **cutout_params)

    excluder = ExclusionContainer(crs=area_crs, res=resolution)

    if "natura" in tech_config and tech_config["natura"]:
        excluder.add_raster(natura, nodata=0, allow_no_overlap=True)

    if "copernicus" in tech_config and tech_config["copernicus"]:
        copernicus = tech_config["copernicus"]
        excluder.add_raster(
            copernicus_map,
            codes=copernicus["grid_codes"],
            invert=True,
            crs=COPERNICUS_CRS,
        )
        if "distance" in copernicus and tech_config["copernicus"]["distance"] > 0:
            excluder.add_raster(
                copernicus_map,
                codes=copernicus["distance_grid_codes"],
                buffer=copernicus["distance"],
                crs=COPERNICUS_CRS,
            )

    if "max_depth" in tech_config:
        # lambda not supported for atlite + multiprocessing
        # use named function np.greater with partially frozen argument instead
        # and exclude areas where: -max_depth > grid cell depth
        func = functools.partial(np.greater, -tech_config["max_depth"])
        excluder.add_raster(
            gebco, codes=func, crs=GEBCO_CRS, nodata=-1000
        )

    if "min_shore_distance" in tech_config:
        buffer = tech_config["min_shore_distance"]
        excluder.add_geometry(country_shapes, buffer=buffer)

    if "max_shore_distance" in tech_config:
        buffer = tech_config["max_shore_distance"]
        excluder.add_geometry(country_shapes, buffer=buffer, invert=True)

    if noprogress:
        logger.info("Calculate landuse availabilities...")
        start = time.time()
        availability = cutout.availabilitymatrix(node_gdf, excluder, **kwargs)

        duration = time.time() - start
        logger.info(f"Completed availability calculation ({duration:2.2f}s)")
    else:
        availability = cutout.availabilitymatrix(node_gdf, excluder, **kwargs)

    duration = time.time() - start
    logger.info(f"Completed availability calculation ({duration:2.2f}s)")

    availability.to_netcdf(output_path)



if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        pgb.streams.wrap_stderr()
        snakemake = mock_snakemake("04_build_renewable_profiles", technology="onwind", weather_year=2012)

    configure_logging(snakemake)
    nprocesses = int(snakemake.threads)
    noprogress = not snakemake.config["atlite"].get("show_progress", False)
    cluster = LocalCluster(n_workers=nprocesses, threads_per_worker=1)
    client = Client(cluster, asynchronous=True)
    kwargs = dict(nprocesses=nprocesses, disable_progressbar=noprogress)

    geo_crs = snakemake.config["crs"]["geo_crs"] # EPSG: 4326
    area_crs = snakemake.config["crs"]["area_crs"] # EPSG:3035
    cutout_params = snakemake.config["atlite"]["cutouts"]
    tech_config = snakemake.config["atlite_config"][snakemake.wildcards.technology]
    country_shapes = snakemake.input.country_shapes
    
    if isinstance(tech_config.get("copernicus", {}), list):
        tech_config["copernicus"] = {"grid_codes": tech_config["copernicus"]}

    if snakemake.wildcards.technology == "offwind":
        # country_shapes = snakemake.input.offshore_shapes
        node = snakemake.input.offshore_shapes
    else:
        # country_shapes = snakemake.input.country_shapes
        node = snakemake.input.onshore_shapes      
        
    calculate_land_availability_matrix(country_shapes,
                                        node,
                                        snakemake.input.natura, 
                                        snakemake.input.gebco, 
                                        snakemake.input.copernicus_map,
                                        snakemake.input.cutout,
                                        cutout_params,
                                        snakemake.config["atlite"]["resolution"], 
                                        snakemake.output.land_availability, 
                                        kwargs)
