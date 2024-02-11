import logging
import time
import os
os.environ["USE_PYGEOS"] = "0"
os.environ["OGR_GEOMETRY_ACCEPT_UNCLOSED_RING"] = "NO"

from pathlib import Path
import geopandas as gpd
from shapely.validation import make_valid
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

from _helpers import configure_logging, iso2_to_iso3_country, iso3_to_iso2_country, plot_style, simplify_polys


def __load_EEZ(iso3:str, EEZ_gpkg:str, geo_crs:str) -> gpd.GeoDataFrame:
    """Load the database of the Exclusive Economic Zones.
    Parameters
    ----------
    iso3 : str
        3-digits country name converted from 2-digits country name
    EEZ_gpkg : str
        path of input data of eez zone
    geo_crs : str
        CRS system

    Returns
    -------
    gpd.GeoDataFrame
        a GeoDataFrame of eez area within the country
    Raises
    ------
    Exception
        The dataset shall be downloaded independently by the user (see guide) or
        together with pypsa-earth package.
    """
    eez_file = Path(EEZ_gpkg)
    if not eez_file.exists():
        raise Exception(
            f"File EEZ {EEZ_gpkg} not found, please download it from \
https://www.marineregions.org/download_file.php?name=World_EEZ_v11_20191118_gpkg.zip \
            and copy it in {EEZ_gpkg}."
        )

    eez_gdf = gpd.read_file(EEZ_gpkg).to_crs(geo_crs)
    eez_gdf.dropna(axis=0, how="any", subset=["ISO_TER1"], inplace=True)
    eez_gdf = eez_gdf[["ISO_TER1", "geometry"]]
    country_iso3 = [
        iso2_to_iso3_country(x) for x in iso3
    ]
    eez_gdf = eez_gdf[
        [any([x in country_iso3]) for x in eez_gdf["ISO_TER1"]]
    ]
    eez_gdf["ISO_TER1"] = eez_gdf["ISO_TER1"].map(
        lambda x: iso3_to_iso2_country(x)
    )
    eez_gdf.reset_index(drop=True, inplace=True)
    eez_gdf.rename(columns={"ISO_TER1": "node_name"}, inplace=True)

    return eez_gdf


def add_node(nuts:str, output_path_shape:str, output_path_node:str)-> gpd.GeoDataFrame:
    """Create the node shapes out of country.
    Parameters
    ----------
    nuts : str
        input path of division of the country 
    output_path_shape : str
        output path to save Shapefile of a whole country
    output_path_node : str
        output path to save Shapefile of node division
    Returns
    -------
    gpd.GeoDataFrame
        table of GeoDataframe with nodes of country
    """
    logger.info("Create a GeoDataFrame of nodal division in the country...")
    start = time.time()
    gdf = (
        gpd.read_file(nuts)
        .rename(columns={
            "id": "nuts",
            "na": "node_name"
        })
        .replace({
            "Ä": "Ae",
            "Ö": "Oe",
            "Ü": "Ue",
            "ä": "ae",
            "ö": "oe",
            "ü": "ue",
            "ß": "ss"
        }, regex=True)
    )

    country_gdf = gdf.dissolve(by="nuts")[["geometry"]].reset_index()
    country_gdf["node_name"] = country_gdf["nuts"][:2].map({"SW": "SE", "AU": "AT"})
    country_gdf = country_gdf.set_index("node_name")
    country_gdf.to_file(output_path_shape, driver="GeoJSON") 
    print(gdf.columns)
    final_gdf = gdf[["node_name" , "nuts", "geometry"]].set_index(["node_name", "nuts"])
    final_gdf = final_gdf.dissolve(by=["node_name", "nuts"], aggfunc="sum")[["geometry"]]
    final_gdf["geometry"].map(simplify_polys)
    final_gdf.to_file(output_path_node, driver="GeoJSON") 
    end = time.time()
    print(f"processing time: {round((end - start), 2)}s")

    return final_gdf


def add_node_map(gdf: gpd.GeoSeries, eez:gpd.GeoSeries, output_path:str):
    """Create the graph of country map divided by nodes.
    Parameters
    ----------
    gdf : Shapefile
        division of the country under nuts2 level 
    eez : Shapefile
        a GeoDataFrame of offshore shape within the country
    output_path : str
        output path to save png graph
    """
    logger.info("Create graphs and GeoJSONs of nodal division in the country...")
    start = time.time()
    gdf = gdf.reset_index()
    eez = eez.reset_index()

    map_gdf = gdf[["node_name", "geometry"]]

    fig, ax = plt.subplots(figsize=(10, 6))

    plot_style()

    # for nuts2 map
    map_gdf.plot(ax=ax, markersize=10, alpha=0.7, categorical=True, legend=False, color="blue")

    plt.title("Nodal Division")
    plt.savefig(output_path, dpi=300)

    end = time.time()
    logger.info(f"Graphs created. Processing time: {round((end - start), 2)}s")


def add_eez(iso2:list, EEZ_gpkg:str, output_path:str, geo_crs:str="EPSG:4326"
) -> gpd.GeoDataFrame:
    """Add cover of offshore wind shapes
    Parameters
    ----------
    iso2 : list
        lsit of 2-digits country name
    EEZ_gpkg : str
        path of input data of eez zone
    output_path : str
        output path to save Shapefile
    geo_crs : str, optional
        CRS system, by default "EPSG:4326"
        
    Returns
    -------
    gpd.GeoDataFrame
        a GeoDataFrame of offshore shape within the country
    """
    logger.info("Create a GeoDataFrame of offshore shape within the country...")
    start = time.time()

    df_eez = __load_EEZ(iso2, EEZ_gpkg, geo_crs)

    ret_df = df_eez[["node_name", "geometry"]]

    for c_code in iso2:
        selection = ret_df.node_name == c_code
        n_offshore_shapes = selection.sum()

        if n_offshore_shapes > 1:
            geom = ret_df[selection].geometry.unary_union
            ret_df.drop(ret_df[selection].index, inplace=True)
            ret_df = ret_df._append(
                {"node_name": c_code, "geometry": geom}, ignore_index=True
            )
    ret_df = ret_df.set_index("node_name")["geometry"]
    
    ret_df_new = ret_df.map(
        lambda x: x
        if x is None
        else simplify_polys(x, minarea=0.0001, tolerance=0.001)
    )
    ret_df_new = ret_df_new.apply(lambda x: x if x is None else make_valid(x))

    ret_df = ret_df_new.dropna()
    ret_df.to_file(output_path, driver="GeoJSON") 
    end = time.time()
    logger.info(f"GeoDataFrame created. Processing time: {round((end - start), 2)}s")

    return ret_df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("01_build_shape")

    configure_logging(snakemake)
    countries = snakemake.config["scenario"]["countries"]

    node = add_node(snakemake.input.nuts, 
                    snakemake.output.country_shapes,
                    snakemake.output.onshore_shapes)

    eez = add_eez(iso2=countries, 
                  EEZ_gpkg=snakemake.input.eez,
                  output_path=snakemake.output.offshore_shapes)

    add_node_map(node, eez,
                snakemake.output.onshore_map)

   