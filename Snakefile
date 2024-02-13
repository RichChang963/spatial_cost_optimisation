from os.path import normpath, exists
from shutil import copyfile, move, rmtree
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

configfile: "config/config_energy.yaml"

ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 4)

CDIR = config["input_dir"] + "geo_output/" + config["scenario"]["output_folder"]
GIDIR = config["input_dir"] + "geo_input"

# main command to enable build_shape process
if config["enable"].get("build_shape", True):
    rule build_shape:
        input: 
            nuts = GIDIR + "/nuts/template_IT.geojson",
            eez = GIDIR + "/eez/eez_v12.gpkg",
        output: 
            country_shapes = CDIR + "/country_shape.geojson",
            onshore_shapes = CDIR + "/node_onshore.geojson",
            offshore_shapes = CDIR + "/node_offshore.geojson",
            onshore_map = CDIR + "/country_map_by_node.png",
        log: "log/b_build_shape.log"
        benchmark: "benchmarks/b_build_shape"
        script:
            "scripts/b_build_shape.py"


# command to enable build_cutout process (without weather_year wildcards)
rule build_cutout:
    input: 
        onshore_shapes = CDIR + "/node_onshore.geojson",
        offshore_shapes = CDIR + "/node_offshore.geojson",
    output: 
        raw_cutout = CDIR + "/netcdf_records/raw_cutout_{weather_year}.nc",
    log: "log/c_build_cutout_{weather_year}.log"
    benchmark: "benchmarks/c_build_cutout_{weather_year}"
    script:
        "scripts/c_build_cutout.py"

# main command to enable build_cutout process with weather_year wildcards
if config["enable"].get("build_year_cutout", True):
    rule build_year_cutout:
        input:
            expand(CDIR + "/netcdf_records/raw_cutout_{weather_year}.nc", **config["technology_mapping"]),


# command to enable build_natura_raster process (without weather_year wildcards)
rule build_natura_raster:
        input:
            shapefiles_land = GIDIR + "/protected_area",
            raw_cutout = CDIR + "/netcdf_records/raw_cutout_{weather_year}.nc",
        output:
            natura_raster =  CDIR + "/natura.tiff",
        log: "log/d_build_natura_raster.log"
        benchmark: "benchmarks/d_build_natura_raster"
        script:
            "scripts/d_build_natura_raster.py"

# main command to enable build_natura_raster with weather_year wildcards
if config["enable"].get("build_all_natura_rasters", True):
    rule build_all_natura_rasters:
        input:
           expand(CDIR + "/netcdf_records/raw_cutout_{weather_year}.nc", **config["technology_mapping"]),


# command to enable build_availability_matrix process (without weather_year & technology wildcards)
rule build_availability_matrix:
    input: 
        natura = lambda w: (
            CDIR + "/natura.tiff"
            if config["atlite_config"][w.technology]["natura"]
            else []),
        gebco=lambda w: (
            GIDIR + "/gebco/GEBCO_2022_TID.nc"
            if "max_depth" in config["atlite_config"][w.technology].keys()
            else []),
        copernicus_map = GIDIR + "/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        country_shapes = CDIR + "/country_shape.geojson",
        onshore_shapes = CDIR + "/node_onshore.geojson",
        offshore_shapes = CDIR + "/node_offshore.geojson",
        raw_cutout = CDIR + "/netcdf_records/raw_cutout_{weather_year}.nc",
    output:
        land_availability = CDIR + "/netcdf_records/land_availability_{technology}_{weather_year}.nc",
    log:
        "log" + "/e_build_availability_matrix_{technology}_{weather_year}.log",
    benchmark:
        "benchmarks" + "/e_build_availability_matrix_{technology}_{weather_year}"
    threads: ATLITE_NPROCESSES
    resources:
        mem_mb=ATLITE_NPROCESSES * 5000,
    script:
        "scripts/e_build_availability_matrix.py"
        
# main command to enable build_availability_matrix with weather_year & technology wildcards
if config["enable"].get("build_all_matrices", True):
    rule build_all_matrices:
        input:
            expand(CDIR + "/netcdf_records/land_availability_{technology}_{weather_year}.nc", **config["technology_mapping"]),


# command to enable build_res_profiles process (without weather_year & technology wildcards)
rule build_res_profiles:
    input: 
        raw_cutout = CDIR + "/netcdf_records/raw_cutout_{weather_year}.nc",
        land_availability = CDIR + "/netcdf_records/land_availability_{technology}_{weather_year}.nc",
        country_shapes = CDIR + "/country_shape.geojson",
        onshore_shapes = CDIR + "/node_onshore.geojson",
        offshore_shapes = CDIR + "/node_offshore.geojson",
    output:
        cutout = CDIR + "/netcdf_records/build_cutout_{weather_year}.nc",
        availability = CDIR + "/Availability_share_{technology}_{weather_year}.csv",
        max_potentials = CDIR + "/Renewables_technical_potentials_{technology}_{weather_year}.csv",
        hourly_profile = CDIR + "/Availability_{technology}_{weather_year}.csv",
        annual_flh = CDIR + "/Annual_FLH_{technology}_{weather_year}.csv",
        profile = CDIR + "/netcdf_records/profile_{technology}_{weather_year}.nc",
        res_map = CDIR + "/{technology}_potentials_{weather_year}.png",
        monthly_cf = CDIR + "/{technology}_monthly_cf_{weather_year}.png",
    log:
        "log" + "/f_build_res_profiles_{technology}_{weather_year}.log",
    benchmark:
        "benchmarks" + "/f_build_res_profiles_{technology}_{weather_year}"
    threads: ATLITE_NPROCESSES
    resources:
        mem_mb=ATLITE_NPROCESSES * 1000,
    script:
        "scripts/f_build_res_profiles.py"

# main command to enable build_res_profiles process with weather_year & technology wildcards
if config["enable"].get("build_all_profiles", True):
    rule build_all_profiles:
        input:
            expand(CDIR + "/netcdf_records/build_cutout_{weather_year}.nc", **config["technology_mapping"]),
            expand(CDIR + "/Availability_share_{technology}_{weather_year}.csv", **config["technology_mapping"]),
            expand(CDIR + "/Renewables_technical_potentials_{technology}_{weather_year}.csv", **config["technology_mapping"]),
            expand(CDIR + "/Availability_{technology}_{weather_year}.csv", **config["technology_mapping"]),
            expand(CDIR + "/Annual_FLH_{technology}_{weather_year}.csv", **config["technology_mapping"]),
            expand(CDIR + "/netcdf_records/profile_{technology}_{weather_year}.nc", **config["technology_mapping"]),
            expand(CDIR + "/{technology}_potentials_{weather_year}.png", **config["technology_mapping"]),
            expand(CDIR + "/{technology}_monthly_cf_{weather_year}.png", **config["technology_mapping"])


# command to pull geojson data from gadm
rule get_gadm:
    output:
        nuts_file = GIDIR + "/nuts/nuts_gadm.geojson",
    log: "log/a_get_gadm.log"
    benchmark: "benchmarks/a_get_gadm"
    script:
        "scripts/a_get_gadm.py"


# command to delete temporary files after completion
rule remove_temporary_file:
# onsuccess:
    run:
        rmtree(".snakemake")
