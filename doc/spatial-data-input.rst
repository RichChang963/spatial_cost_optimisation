.. _spatial-data-input: 

##########################################
Spatial Data Input
##########################################


This is documentation for how to create the spatial intput data for country map divided by ndoes (country_map_by_node.png), hourly renewable profiles (Availability.csv) and potential data (Renewables_technical_potential.csv) from atlite package if you do not have your own data or if you do not want to use the default data from `data/global csv templates/ <https://github.com/agoenergy/pypsa-agora/tree/main/data/global%20csv%20templates>`_. 

The first step in Pypsa workflow is to prepare network inputs. Because most of input data will be in geo-format, we have to collect the data by downloading from each website and copy into the `data/geo_input` that we will use as an inputs folders. The structure of the folder inside `data/geo_input` is displayed as follows:


`/node_division.csv`
====================
**Definition**

The node_division.csv consists of several columns (node_name,zone,nuts1,nuts1_eng,other_name). These define the nodes within the region that we will use inside Pypsa network. Please fill in the information according to your target region/country before running the next step.

- ``node_name``: The exact node name that will be used in the model (e.g. VN_SO)
- ``zone``: The general region name of the node (e.g. Southern is the zone name of VN_SO)
- ``nuts1``: The corresponding provinces/states level within the node (e.g. Long An belongs to the VN_SO node). You can refer to `Nomenclature of territorial units for statistics <https://ec.europa.eu/eurostat/web/nuts/background>`_ for more details on the naming rules of different nuts/admin level
- ``nuts1_eng``: The english name of the nuts1 provinces/states in case they are shown in local language
- ``other_name``: Other alias of the province/state


**Instructions on how to fill in the csv**

All values are in string format and hence there is no speical restrictions on filling-in data process. Just remember to cover all provinces/states within the region/country so that the visualisation will be completed.


`/admin_level/`
====================
**Definition**

The data of this folder shall be shapefiles from `GADM <https://gadm.org/download_country.html>`_ and the definition of different nuts/admin levels can be refered from `Nomenclature of territorial units for statistics <https://ec.europa.eu/eurostat/web/nuts/background>`_. There might be some extra specific areas for visualization purpose which will not be used in modelling. 

**Instructions on how to fill in the csv**

1. When downloading the shapefiles from GADM, there will be several files by different nuts levels. For consistency, please rename all the shapefiles in nuts2 level to be [`gadm36_2`] so that the paths can be read properly in snakmake rules. 

2. In the case of extra areas, the shapefiles of these areas shall be also saved in this folder. The path to read the shapefiles shall be also record in `snakemake <https://github.com/agoenergy/pypsa-agora/blob/main/Snakefile>`_ files as following section:



    if config['enable'].get('build_node_n_land_map', True):
        rule build_node_n_land_map:
            input: 
                ...
                extra_nuts= GIDIR + "/admin_level/YOUR_SHAPEFILE_NAME" 

3. Please replace the whole path (`GIDIR + "/admin_level/extra_nuts_vn_humandata.shp"`) with **None** if you do not wish to use it.


`/eez/`
====================
**Definition**

EEZ represents "Exclusive Economic Zone" which is used for offshore wind energy assessment. The data is from `Marineregions.org <https://www.marineregions.org/downloads.php>`_.

**Instructions on how to fill in the csv**

While downloading data from `Marineregions.org <https://www.marineregions.org/downloads.php>`_, please select the [shapefile] in `World EEZ v11` and enter personal email and information (GDPR purpose) to download and save the files inside this folder.


`/gebco/`
====================
**Definition**

`GEBCO <https://www.gebco.net/data_and_products/gridded_bathymetry_data/>`_  data covers underwater depth in meters on a 15 arc-second interval grid. This will be used to exclude the area that are too deep for offshore wind turbines to be built.

**Instructions on how to fill in the csv**

1. Instead of downloading the whole dataset from GEBCO which might take up a very large space of your memory, enter the custormized user panel `GEBCO User Panel <https://download.gebco.net/>`_, enter the corresponding offshore eez boundary of the target region/country including maximum and minimum of longtitude and latitdue.

2. In the SELECT FORMATS section, click the checkbox of [Grid] in [2D netCDF] row. 

3. In the YOUR DATA SELECTION section, click [add to basket] and then complete the download process. Put the `GEBCO_2020_2D.nc` and `GEBCO_2020_2D.tif` in this folder to complete the input process.


`/landuse/`
====================
**Definition**

The landuse data represents the land area which are occupied for other types of applications including residential, industrial, military areas, waterways, railways, etc. The data is from `OpenSteetpMap Data Extracts <http://download.geofabrik.de/>`_. It is derived from `OpenSteetpMap API <https://wiki.openstreetmap.org/wiki/API>`_

**Instructions on how to fill in the csv**

Please download the shapefiles and save them in the folder. For consistency, please change the names in a systematic way. The default names include:

- ``buildings``
  
- ``landuse``
  
- ``railways``
  
- ``roads``
  
- ``waterways``


`/protected_area/`
====================
**Definition**

Protected area mostly covers vegetation or national parks. The source is from `WDPA <https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA>`_.

**Instructions on how to fill in the csv**

1. Once entering the website `WDPA <https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA>`_, please enter the target region/country to find out the location and download the shapefile.

2. For consistency, please take out all the shapefiles within each folder and save it in current folder.


