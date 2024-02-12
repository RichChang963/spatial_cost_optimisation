# SPDX-FileCopyrightText: : 2024 The spatial_cost_optimisation Author
#
# SPDX-License-Identifier: MIT

import logging
import os
os.environ["USE_PYGEOS"] = "0"
os.environ["OGR_GEOMETRY_ACCEPT_UNCLOSED_RING"] = "NO"

import requests
import shutil
from pathlib import Path
import zipfile
import geopandas as gpd
import pandas as pd

logger = logging.getLogger(__name__)

from _helpers import configure_logging, iso2_to_iso3_country


def __get_GADM_filename(iso2_code:str, nuts_level:int=1)-> str:
    """Function to get the GADM filename given the country code.

    Parameters
    ----------
    iso2_code : str
        iso2 of the country
    nuts_level : int, optional
        level of nodal division, by default 1

    Returns
    -------
    str
        filename for GADM source
    """
    special_codes_GADM = {
        "XK": "XKO",  # kosovo
        "CP": "XCL",  # clipperton island
        "SX": "MAF",  # sint maartin
        "TF": "ATF",  # french southern territories
        "AX": "ALA",  # aland
        "IO": "IOT",  # british indian ocean territory
        "CC": "CCK",  # cocos island
        "NF": "NFK",  # norfolk
        "PN": "PCN",  # pitcairn islands
        "JE": "JEY",  # jersey
        "XS": "XSP",  # spratly
        "GG": "GGY",  # guernsey
        "UM": "UMI",  # united states minor outlying islands
        "SJ": "SJM",  # svalbard
        "CX": "CXR",  # Christmas island
    }

    if iso2_code in special_codes_GADM:
        return f"gadm41_{special_codes_GADM[iso2_code]}_{str(nuts_level)}"
    else:
        return f"gadm41_{iso2_to_iso3_country(iso2_code)}_{str(nuts_level)}"


def download_GADM(iso2_code, update=False, nuts_level:int=1) -> tuple:
    """
    Download gpkg file from GADM for a given country code.

    Parameters
    ----------
    iso2_code : str
        Two letter country codes of the downloaded files
    update : bool
        Update = true, forces re-download of files
    nuts_level: int
        level of nodal division, by default 1

    Returns
    -------
    tuple[str, str]
        filename and file path of GADM data
    """
    GADM_filename = __get_GADM_filename(iso2_code, nuts_level=nuts_level)
    GADM_url = f"https://geodata.ucdavis.edu/gadm/gadm4.1/json/{GADM_filename}.json.zip"
    GADM_inputfile_path = os.path.join(
        os.getcwd(),
        "data",
        "geo_input",
        "gadm",
        GADM_filename + ".json.zip",
    )  # Input filepath

    if not os.path.exists(GADM_inputfile_path) or update is True:
        #  create data/osm directory
        os.makedirs(os.path.dirname(GADM_inputfile_path), exist_ok=True)

        try:
            r = requests.get(GADM_url, stream=True, timeout=300)
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
            raise Exception(
                f"GADM server is down at {GADM_url}. Data needed for building shapes can't be extracted.\n\r"
            )
        except Exception as exception:
            raise Exception(
                f"An error happened when trying to load GADM data by {GADM_url}.\n\r"
                + str(exception)
                + "\n\r"
            )
        else:
            with open(GADM_inputfile_path, "wb") as f:
                shutil.copyfileobj(r.raw, f)

    return GADM_inputfile_path, GADM_filename


def clean_geojson(GADM_inputfile_path:str, GADM_filename:str, nuts_level:int=1
) -> gpd.GeoDataFrame:
    """Unzip ZIP file, read, clean, and save Geojson file.

    Parameters
    ----------
    GADM_inputfile_path : str
        filename of GADM data
    GADM_filename : str
        file path of GADM data
    nuts_level : int, optional
        level of nodal division, by default 1

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame of nodal division
    """
    archive = zipfile.ZipFile(GADM_inputfile_path, "r")
    node = gpd.read_file(archive.read(f"{GADM_filename}.json").decode("utf-8"))

    if nuts_level == 1 and "ISO_1" in node.columns.unique():
        filter_col = ["ISO_1", "NAME_1", "COUNTRY", "geometry"]
    else:
        filter_col = [
            f"GID_{str(nuts_level)}", f"NAME_{str(nuts_level)}", "geometry"
        ]
    node = node[filter_col].rename(columns={
        "ISO_1": "nuts", 
        f"GID_{str(nuts_level)}": "nuts", 
        f"NAME_{str(nuts_level)}": "node_name",
        "COUNTRY": "country"
    })[["country", "node_name", "nuts", "geometry"]]

    return node.set_index("nuts")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("a_get_gadm")

    configure_logging(snakemake)
    country = snakemake.config["scenario"]["countries"]
    nuts_level = snakemake.config["scenario"]["gadm_nuts_level"]

    nuts_gdf = gpd.GeoDataFrame()
    for c in country:
        GADM_inputfile_path, GADM_filename = download_GADM(c, nuts_level)

        node = clean_geojson(GADM_inputfile_path, GADM_filename, nuts_level)

        nuts_gdf = pd.concat([nuts_gdf, node], axis=1)
    
    nuts_gdf.to_file(snakemake.output.nuts_file)

    # remove ZIP files
    shutil.rmtree(str(Path(GADM_inputfile_path).parents[0]))