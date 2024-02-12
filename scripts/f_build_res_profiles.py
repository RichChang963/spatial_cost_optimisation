# SPDX-FileCopyrightText: : 2024 The spatial_cost_optimisation Author
#
# SPDX-License-Identifier: MIT

import logging
import time
import os
os.environ["USE_PYGEOS"] = "0"
os.environ["OGR_GEOMETRY_ACCEPT_UNCLOSED_RING"] = "NO"
import progressbar as pgb
from dask.distributed import Client, LocalCluster
import collections # fix atlite iterable issue. Import it before importing atlite package
collections.Iterable = collections.abc.Iterable # fix atlite iterable issue.
import atlite
import geopandas as gpd
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from _helpers import configure_logging, plot_style


logger = logging.getLogger(__name__)
pd.options.mode.chained_assignment = None


tech_mapping = {
    "solar": "PHOT",
    "onwind": "WTON",
    "offwind": "WTOF"
}


def _local_time_convert(data:pd.DataFrame, utc_zone:int, weather_year:str
) -> pd.DataFrame:
    """Convert timezone from UTC+0 to local time zone.

    Parameters
    ----------
    data : pd.DataFrame
        input Dataframe
    utc_zone : int
        utc time zone to convert time
    weather_year : str
        start date of the selected year

    Returns
    -------
    pd.DataFrame
        input Dataframe with local time zone
    """
    data = data.set_index("time")
    data = pd.concat(objs=[data.iloc[24-utc_zone:],data.iloc[:24-utc_zone]])
    data["time"] = (pd.date_range(
        f"{str(weather_year)}-01-01 00:00:00", 
        f"{str(weather_year)}-12-31 23:00:00", 
        freq="H"))
    data = data.set_index("time")

    return data


def calculate_renewables_profile(cutout_name:str, cutout_params:dict, 
avail_matrix:xr.Dataset, client:dict) -> tuple:
    """Calculate hourly profiles, max_potentials of different renewable technologies.

    Parameters
    ----------
    cutout_name : str
        path of the cutout file
    cutout_params: dict
        a dictionary of cutout variables for era5 api
    avail_matrix : netCDF
        netcdf file of availability matrix 
    client : dict
        dictionary of client information for clustering process from dask.distributed

    Returns
    -------
    Tuple[pd.DataFrame, list, atlite.Cutout]
        a DataFrame of potential data
        a DataArray of area data
        atlite.Cutout file of the whole cutout boundary
    """
    capacity_per_sqkm = tech_config["capacity_per_sqkm"]
    evaluation = tech_config["evaluation"]
    resource = tech_config["resource"]

    cutout = atlite.Cutout(cutout_name, **cutout_params) 

    logger.info(f"Calculate renewable profiles (method '{evaluation}')...")
    area = cutout.grid.to_crs(3035).area / 1e6
    area = (xr.DataArray(area.values.reshape(cutout.shape),
            [cutout.coords['y'], cutout.coords['x']]))
    avail_matrix = xr.open_dataarray(avail_matrix)
    potential_matrix = capacity_per_sqkm * avail_matrix * area / 1000 # unit in GW
    resource["dask_kwargs"] = {"scheduler": client}
    avail_matrix.close()

    return potential_matrix, area, cutout


def plot_potential_map(node:str, potential_matrix:pd.DataFrame, technology:str, 
cutout:atlite.Cutout, potential_map_output:str) -> gpd.GeoDataFrame:
    """Plot potential map

    Parameters
    ----------
    node : str
        path of node division of the country
    potential_matrix : pd.DataFrame
        a DataFrame of potential data
    technology : str
        name of technology
    cutout : atlite.Cutout
        cutout file from function extract_cutout_feature
    potential_map_output : str
        output path to save png graph

    Returns
    ----------
    gpd.GeoDataFrame
        a GeoDataFrame of node division
    """
    node_gdf = gpd.read_file(node).set_index("node_name")
    fig, ax = plt.subplots(figsize=(8, 8))
    plot_style()
    potential_matrix.sum("node_name").plot(cmap='Purples')
    node_gdf.plot(ax=ax, edgecolor='k', color='None', lw=0.1)
    cutout.grid.plot(ax=ax, color='None', edgecolor='grey', ls=':', lw=0.1)
    plt.title(f"{technology} energy potentials: total {potential_matrix.sum():2.2f} GW")
    plt.savefig(potential_map_output, dpi=300)

    return node_gdf
    

def generate_hourly_profile(node_gdf:gpd.GeoDataFrame, area:xr.DataArray, 
avail_matrix:xr.Dataset, cutout:atlite.Cutout, weather_year:int, utc_timezone:int,
technology:str, profile_output_path:xr.Dataset) -> pd.DataFrame:
    """Create table of hourly capacity factor

    Parameters
    ----------
    node_gdf : gpd.GeoDataFrame
        a GeoDataFrame of node division
    area: xr.DataArray
        a DataArray of area data
    avail_matrix : netCDF
        netcdf file of availability matrix 
    cutout : atlite.Cutout
        cutout file from function extract_cutout_feature
    weather_year: int
        target year
    utc_timezone: int
        utc time zone
    technology: str
        onshore wind (WTON), offshore wind (WTOF), or solar pv (PHOT)
    profile_output_path: netCDF
        a netCDF of whole renewable profiles

    Returns
    -------
    pd.DataFrame
        a DataFrame of hourly capacity factor
    """
    capacity_per_sqkm = tech_config["capacity_per_sqkm"]
    tech_method = tech_config["tech_method"]
    correction_factor = tech_config["correction_factor"]
    evaluation = tech_config["evaluation"]
    resource = tech_config["resource"]

    avail_matrix = xr.open_dataarray(avail_matrix)
    availability = avail_matrix @ area
    potential = capacity_per_sqkm * avail_matrix.sum("node_name") * area # unit in MW
    func = getattr(cutout, tech_method)
    capacity_factor = correction_factor * func(capacity_factor=True, **resource)
    layout = capacity_factor * area * capacity_per_sqkm
    profile, capacities = func(matrix=avail_matrix.stack(spatial=['y','x']),
                                layout=layout, index=node_gdf.index,
                                per_unit=True, return_capacity=True, **resource)  
    if evaluation == "simple":
        p_nom_max = capacity_per_sqkm * avail_matrix @ area
    elif evaluation == "conservative":
        max_cap_factor = capacity_factor.where(avail_matrix != 0).max(["x", "y"])
        p_nom_max = capacities / max_cap_factor
    else:
        raise AssertionError(
            "Config key `potential` should be one of 'simple' "
            f"(default) or 'conservative', not '{evaluation}'"
        )
    
    avail_matrix.close()
    p_nom_max_df = pd.DataFrame(p_nom_max.to_pandas(), columns = ["p_nom_max_MW"]) # MW
    avail_df = pd.DataFrame(availability.to_pandas(), columns = ["availability_%"])
    print("Total p_nom_max: ", p_nom_max.sum().to_pandas(), "MW") # MW

    complie_ds = xr.merge([
        (correction_factor * profile).rename("profile"), # MW
        capacities.rename("weight"), # MW
        p_nom_max.rename("p_nom_max"), # MW
        potential.rename("potential"), # MW
    ])

    complie_ds.to_netcdf(profile_output_path)

    hourly_cf_df = complie_ds["profile"].to_pandas().squeeze().reset_index()
    hourly_cf_df = _local_time_convert(hourly_cf_df, utc_timezone, weather_year)

    cf_df = pd.DataFrame()
    cf_df = pd.concat([cf_df, hourly_cf_df], axis=1)

    end_hour = 8760 if len(hourly_cf_df) == 8760 else 8784
    cf_df["node"] = np.arange(0, end_hour, 1)
    new_cf_df = cf_df.set_index("node").T 

    tech_name = tech_mapping.get(technology)
    new_cf_df.insert(0, "archetypes", tech_name)
    new_cf_df.index.names = ["node"]
    p_nom_max_df.insert(0, "archetypes", tech_name)
    avail_df.insert(0, "archetypes", tech_name)

    new_nom_max_df = (
        p_nom_max_df.reset_index()
        .rename(columns = {
            "node_name": "node",
            "p_nom_max_MW": "technical potential [MW]"
        })
        .astype({"node": "string"})
        .sort_values(by=["node"])
        .set_index("node")
    )

    new_avail_df = (
        avail_df.reset_index()
        .rename(columns = {"node_name": "node"})
        .astype({"node": "string"})
        .sort_values(by=["node"])
        .set_index("node")
    )

    if technology != "offwind":
        country_df =  (
            node_gdf.reset_index()
            .rename(columns={
                "na": "node", 
                "node_name": "node", 
                "nuts1": "id"}
            )
            .sort_values(by=["node"])
        )
        country_df["country"] = country_df["id"].str[:2]
        country_df = (
            country_df[["id", "country", "node"]]
            .astype({"node": "string"})
            .set_index("node")
        )

        new_cf_df = pd.merge(
            new_cf_df, country_df, 
            left_index=True, right_index=True, how="outer"
        )
        new_nom_max_df = pd.merge(
            new_nom_max_df, country_df, 
            left_index=True, right_index=True, how="outer"
        )
        new_avail_df = pd.merge(
            new_avail_df, country_df, 
            left_index=True, right_index=True, how="outer"
        )
        index_col = ["node", "id", "country", "archetypes"]
    else:
        index_col = ["node", "archetypes"]

    flh_df = pd.DataFrame(new_cf_df.reset_index().set_index(index_col).sum(axis=1))
    flh_df.columns.names = ["Annual_FLH"]
    flh_df = flh_df.reset_index().set_index("node")

    new_cf_df.to_csv(snakemake.output.hourly_profile, encoding="utf-8-sig", sep=",")
    flh_df.to_csv(snakemake.output.annual_flh, encoding="utf-8-sig", sep=",")
    new_nom_max_df.to_csv(snakemake.output.max_potentials, encoding="utf-8-sig", sep=",")
    new_avail_df.to_csv(snakemake.output.availability, encoding="utf-8-sig", sep=",")

    return new_cf_df


def plot_hourly_graph(new_cf_df:pd.DataFrame, technology:str, 
monthly_profile_graph_output:str):
    """Plot hourly capcaity factor graph

    Parameters
    ----------
    new_cf_df : pd.DataFrame
        a DataFrame of hourly capacity factor
    technology : str
        name of technology
    monthly_profile_graph_output : str
        output path to save png graph
    """
    new_cf_df.resample("1M").mean().plot(figsize=(8,4.5)) 
    plot_style()
    plt.ylabel("Capacity Factor")
    plt.title(f"monthly capacity factor of {technology}")
    plt.savefig(monthly_profile_graph_output, dpi=300)



if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        pgb.streams.wrap_stderr()
        snakemake = mock_snakemake("f_build_res_profiles", technology="onwind", weather_year=2012)

    configure_logging(snakemake)
    nprocesses = int(snakemake.threads)
    cluster = LocalCluster(n_workers=nprocesses, threads_per_worker=1)
    client = Client(cluster, asynchronous=True)

    cutout_path = snakemake.input.cutout
    avail_matrix = snakemake.input.land_availability
    potential_map_output = snakemake.output.res_map
    monthly_profile_graph_output = snakemake.output.monthly_cf

    # country_name = snakemake.config["scenario"]["countries"]
    tech_config = snakemake.config["atlite_config"][snakemake.wildcards.technology]
    cutout_params = snakemake.config["atlite"]["cutouts"]
    utc_timezone = snakemake.config["technology_mapping"]["utc"]
    
    if snakemake.wildcards.technology == "offwind":
        country_shapes = snakemake.input.offshore_shapes
        node = snakemake.input.offshore_shapes
    else:
        country_shapes = snakemake.input.country_shapes
        node = snakemake.input.onshore_shapes      

    potential_matrix, area, cutout = calculate_renewables_profile(cutout_path, 
                                cutout_params, 
                                avail_matrix, 
                                client)

    node_gdf = plot_potential_map(node, 
                                  potential_matrix, 
                                  snakemake.wildcards.technology, 
                                  cutout, 
                                  potential_map_output)


    new_cf_df = generate_hourly_profile(node_gdf, 
                                                area, 
                                                avail_matrix, 
                                                cutout, 
                                                snakemake.wildcards.weather_year,
                                                utc_timezone,
                                                snakemake.wildcards.technology,
                                                snakemake.output.profile)

    plot_hourly_graph(new_cf_df, 
                      snakemake.wildcards.technology,
                      monthly_profile_graph_output)