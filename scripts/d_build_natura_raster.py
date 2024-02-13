# SPDX-FileCopyrightText: : 2024 The spatial_cost_optimisation Author
#
# SPDX-License-Identifier: MIT

import logging
import os
from os import listdir
os.environ["USE_PYGEOS"] = "0"
os.environ["OGR_GEOMETRY_ACCEPT_UNCLOSED_RING"] = "NO"

import atlite
import geopandas as gpd
import pandas as pd
import multiprocessing as mp
import rasterio as rio
from _helpers import configure_logging
from rasterio.features import geometry_mask
from rasterio.warp import transform_bounds

logger = logging.getLogger(__name__)

from _helpers import configure_logging, iso2_to_iso3_country, simplify_polys

CUTOUT_CRS = "EPSG:4326"


def determine_cutout_xXyY(cutout_name:str):
    """Define the coordinate boundary of the cutout area

    Parameters
    ----------
    cutout_name : str
        the name of the cutout

    Returns
    -------
    tuple
        values of boundary of cutout
    """
    logger.info("Stage 1/5: Determine cutout boundaries")
    cutout = atlite.Cutout.from_netcdf(cutout_name)
    assert cutout.crs == CUTOUT_CRS
    x, X, y, Y = cutout.extent
    dx, dy = cutout.dx, cutout.dy
    cutout_xXyY = [
        x - dx / 2.0, 
        X + dx / 2.0, 
        y - dy / 2.0, 
        Y + dy / 2.0, 
    ]

    return cutout_xXyY


def get_transform_and_shape(bounds:list, res:int):
    """transform the graph from shapefile into raster file

    Parameters
    ----------
    bounds : list
        array of cooridnate boundary of cutout
    res : int
        resolution

    Returns
    -------
    tuple[int, int]
        values of transform and new shape data in raster map
    """
    logger.info("Stage 2/5: Get transform and shape")
    left, bottom = [(b // res) * res for b in bounds[:2]]
    right, top = [(b // res + 1) * res for b in bounds[2:]]
    # "latitude, longitude" coordinates order
    shape = int((top - bottom) // res), int((right - left) // res)
    transform = rio.Affine(res, 0, left, 0, -res, top)

    return transform, shape


def uniform_protected_area(country_name:list, input_path:str, natura_crs:int
) -> gpd.GeoDataFrame:
    """Add land use cover
    Parameters
    ----------
    country_name: list
        an array of iso2 code of countries
    input_path : str
        folder path of shapefiles of protected_area data
    natura_crs: int
        coordinate reference systems of natura (3035)
    Returns
    -------
    gpd.GeoDataFrame
        a GeoDataFrame including all protected area
    """
    logger.info("Stage 3/5: Uniform protected area into one shapefile.")
    protected_area_files = (
        [f for f in listdir(input_path) if any(f.endswith(name) for name in "shp")]
    )
    protected_area_paths = []
    for f in protected_area_files:
        protected_area_paths.append(input_path+"/"+f)

    # reading shapefiles with multiprocessing
    pool=mp.pool.ThreadPool(min(mp.cpu_count(), len(protected_area_paths), 4))
    protected_frames=pool.map(gpd.read_file, protected_area_paths, chunksize=1)
    pool.close()

    final_gdf = gpd.GeoDataFrame()

    for i in range(len(protected_frames)):
        for c in country_name:
            country_df = (protected_frames[i].loc[protected_frames[i]["ISO3"] == 
                                    iso2_to_iso3_country(c)])
            country_df = (country_df[["DESIG_ENG", "geometry"]]
                        .rename(columns={"DESIG_ENG": "type"}))
            final_gdf = gpd.GeoDataFrame(pd.concat([final_gdf, country_df]))

    final_gdf["geometry"].map(simplify_polys)
    final_gdf = final_gdf.to_crs(natura_crs)

    return final_gdf


def save_raster_file(shapes:gpd.GeoDataFrame, out_shape:int, transform:int, 
output_path:str):
    """Mask the geometry and save the raster file.

    Parameters
    ----------
    shapes : gpd.GeoDataFrame
        Geojson file of the natura raster.
    out_shape : int
        values from the output of get_transform_and_shape function.
    transform : int
        values from the output of get_transform_and_shape function.
    output_path : str
        file path to save the raster file.
    """
    logger.info("Stage 4/5: Mask geometry")
    raster = ~geometry_mask(shapes.geometry, out_shape, transform)
    raster = raster.astype(rio.uint8)

    logger.info("Stage 5/5: Export as .tiff")
    with rio.open(
        output_path,
        "w",
        driver="GTiff",
        dtype=rio.uint8,
        count=1,
        transform=transform,
        crs=natura_crs,
        compress="lzw",
        width=raster.shape[1],
        height=raster.shape[0],
    ) as dst:
        dst.write(raster, indexes=1)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("d_build_natura_raster", weather_year=2012)
    configure_logging(snakemake)

    natura_crs = snakemake.config["crs"]["area_crs"] # 3035
    country_name = snakemake.config["scenario"]["countries"]
    resolution = snakemake.config["atlite"]["resolution"]

    cutouts = snakemake.input.raw_cutout
    xs, Xs, ys, Ys = zip(
        *(determine_cutout_xXyY(cutout) for cutout in cutouts)
    )
    bounds = transform_bounds(
        CUTOUT_CRS, natura_crs, min(xs), min(ys), max(Xs), max(Ys)
    )
    transform, out_shape = get_transform_and_shape(bounds, res=resolution)
    shapefiles = uniform_protected_area(country_name, 
                                        snakemake.input.shapefiles_land,
                                        natura_crs)
    save_raster_file(shapefiles, 
                    out_shape, 
                    transform, 
                    snakemake.output.natura_raster)