.. _build-spatial-res-profiles: 

######################################################
Build Spatial Renewable Profile from Atlite Package
######################################################


This is documentation for how to create hourly renewable profiles ``Availability.csv`` and potential data ``Renewables_technical_potential.csv`` from atlite package if you do not have your own data or if you do not want to use the default data from `data/global csv templates/ <https://github.com/agoenergy/pypsa-agora/tree/main/data/global%20csv%20templates>`_. Current version covers onshore wind, offshore wind, and photovoltaic. This step shall be executed after `build_skeleton.py` is created. Make sure that the spatial/geo-input data has been filled before starting this process. You can refer to `spatial-data-input.rst` for more information.

Before starting the process, some settings need to be checked in `config.yaml` and `Snakefile`. The installation of `ERA5` is also required to access the api.


config.yaml
====================
There are three section in this file that needs to be adjusted. 

1. In [snapshots] section, you can write different weather_year in `string` format. The hourly profile will based on this time scale. [utc] requires an `integer` format input to convert the hourly profile into local time zone.

    ``weather_year: "2019"`` 

    ``utc: 7``


2. In [enable] section, you can enable the script by give `true` value. If you turn ``build_availability_profile`` to be `true`, it requires ``build_node_n_land_map`` to be `true` and executed first since the inputs of ``build_availability_profile`` comes from the outputs of ``build_node_n_land_map``.

    ``build_node_n_land_map: false`` # or true
 
    ``build_availability_profile: false`` # or true


3. In ``atlite`` section, the resolution of availability matrix and energy potential is defiend based on square meters. Default value is 100 meters

    ``resolution: 50`` # default: 100


4. In ``atlite_config`` section, different technologies are defined with different characteristics which affect the final energy potentials.  

- ``capacity_per_sqkm``: the allowable density of wind turbine or solar pv panel placement.
- ``turbine``: the model of turbines. The database of turbine model can be referred from `atlite/resources/windturbine/ <https://github.com/PyPSA/atlite/tree/master/atlite/resources/windturbine>`_

- ``evaluation``: either ``simple`` or ``conservative`` in `string` format. According to `Pypsa-eur/configuration <https://pypsa-eur.readthedocs.io/en/latest/configuration.html>`_. ``simple`` adds up the installable potentials of the individual grid cells. If the model comes close to this limit, then the time series may slightly overestimate production since it is assumed the geographical distribution is proportional to capacity factor. ``conservative`` assertains the nodal limit by increasing capacities proportional to the layout until the limit of an individual grid cell is reached.

    [tech]:
        capacity_per_sqkm: 3 # Allowable density of wind turbine placement. # ScholzPhd Tab 4.3.1: 10MW/km^2

        resource:

            turbine: Vestas_V112_3MW

            evaluation: simple # or conservative


Snakefile
====================
The only note in checking the `Snakefile` is in the `extra_nut` input. As mentioned in `spatial-data-input.rst`, if you Please replace the whole path (`GIDIR + "/admin_level/extra_nuts_vn_humandata.shp"`) with **None** if you do not wish to use it.


ERA5 Installation
====================
ERA5 api requires **CDS API key** to access the api. Please refer to the instruction and introduction from `atlite/Creating a Cutout with ERA5 <https://atlite.readthedocs.io/en/latest/examples/create_cutout.html>`_. After getting **CDS API key**, please save the file under the main directory 
`C://Users/your_user_name/.cdsapirc`.


Once all the above pre-settings are completed, the model to evaluate energy potentials and hourly profiles of onshore wind, offshore wind and photovoltaic technologies.

build_node_n_land_map
==========================
This is the first script to be executed. The output of this script will be saved in `data/COUNTRY_ISO2` including:

- ``country_map_by_node.png``
- ``node_onshore.geojson``: a shapefile of onshore area (will be used as a input in the next script)
- ``node_offshore.geojson``: a shapefile of offshore area (will be used as a input in the next script)
- ``land_availability_cutout_onshore.nc``: a NetCDF file of land availability matrix of onshore area (will be used as a input in the next script)
- ``land_availability_cutout_offshore.nc``: a NetCDF file of land availability matrix of offshore area (will be used as a input in the next script)

build_availability_profile
================================
This is the second script to be executed after ``build_node_n_land_map`` is completed. The output of this script will be saved in `data/COUNTRY_ISO2`  including:

- ``solar_potentials.png``
- ``solar_monthly_cf.png``: monthly capacity factor of solar by node
- ``onshore_wind_potentials.png``
- ``onshore_wind_monthly_cf.png``: monthly capacity factor of onshore wind by node
- ``offshore_wind_potentials.png``
- ``offshore_wind_monthly_cf.png``: monthly capacity factor of offshore wind by node

The output of maximum energy potentials and hourly profile (capacity factor) of different technologies will replace the `Renewables_technical_potential.csv` and `Availability.csv` (created by ``build_skeleton``) respectively. 
