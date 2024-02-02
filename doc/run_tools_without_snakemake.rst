.. 
    SPDX-FileCopyrightText: 2022 Agora Energiewende authors
    
    SPDX-License-Identifier: CC-BY-4.0
.. run_tools_without_snakemake: 

######################################################
Run the whole tools without snakemake rules
######################################################


This is documentation for how to create hourly renewable profiles ``Availability.csv`` and potential data ``Renewables_technical_potential.csv`` from atlite package if you do not want to use `snakemake rules` and `config settings`. Make sure that the spatial/geo-input data has been filled before starting this process. You can refer to `spatial-data-input.rst` for more information.
Please note that this module does not record the process in the `log file` and the variables need to adjust manually in each script before running at the terminal.

Manual inputs
====================
The variables need to be adjusted in each script before starting. This will make sure that the output will match your request of node division.
`build_node_n_land_map <https://github.com/agoenergy/Atlite-agora/blob/directly-command-without-snakemake/without_snakemake/build_node_n_land_map.py#L25>`_

`build_availability_profile <https://github.com/agoenergy/Atlite-agora/blob/directly-command-without-snakemake/without_snakemake/build_availability_profile.py#L20`_

Once all the above pre-settings are completed, the model to evaluate energy potentials and hourly profiles of onshore wind, offshore wind and photovoltaic technologies.

Commands to run the scripts
====================
1. build_node_n_land_map

.. code:: bash

    .../Atlite-agora% python main.py build-map-n-matrix


2. build_node_n_land_map

.. code:: bash

    .../Atlite-agora% python main.py build-profile
    

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

