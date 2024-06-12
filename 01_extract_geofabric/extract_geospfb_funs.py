# functions to extract the MERIT Basins geospatial fabric

# load the needed libraries
import geopandas as gpd # version 0.14.0
import pandas as pd # version 1.4.0
import numpy as np # version 1.22.2
import matplotlib.pyplot as plt # version 3.5.1

from shapely.geometry import Point # version 2.0.1

import hydrant.topology.geom as gm # version 0.1.0-dev1

import subprocess # built-in Python 3.10.2
import os # built-in Python 3.10.2
import glob # built-in Python 3.10.2

##########################

# Read the MERIT Basins geofabric
# This part is static as it reads the entire MERIT geofabric (should only run once)
def read_MERIT_basins(merit_basins_root_path):
    merit_basins_geom_path = os.path.join(merit_basins_root_path, 'pfaf_level_02')
    merit_basins_nca_path = os.path.join(merit_basins_root_path, 'coastal_hillslopes')
    # Reading `MERIT-Basins` Geospatial Fabric Dataset
    ##########################
    ## `MERIT-Basins` Geospatial Layers
    cat_files = [
        'cat_pfaf_71_MERIT_Hydro_v07_Basins_v01_bugfix1.shp',
    ]
    # rivers (river segments)
    riv_files = [
        'riv_pfaf_71_MERIT_Hydro_v07_Basins_v01_bugfix1.shp',
    ]
    # non-contributing catchments (those without any river segments defined for them)
    nca_files = [
        'hillslope_71_clean.shp',
    ]

    # reading in data in an iterative manner
    cat = pd.concat([gpd.read_file(os.path.join(merit_basins_geom_path, f)) for f in cat_files])
    riv = pd.concat([gpd.read_file(os.path.join(merit_basins_geom_path, f)) for f in riv_files])
    nca = pd.concat([gpd.read_file(os.path.join(merit_basins_nca_path, f)) for f in nca_files])
    ##########################
    # reproject `MERIT-Basins` to EPSG `4326`. 
    # specifying epsg:4326 for all the MERIT-Basins layers
    cat.set_crs(epsg=4326, inplace=True)
    nca.set_crs(epsg=4326, inplace=True)
    riv.set_crs(epsg=4326, inplace=True)

    # Show the EPSG of all geospatial layers
    # print(f'`cat` CRS: {cat.crs}')
    # print(f'`riv` CRS: {riv.crs}')
    # print(f'`nca` CRS: {nca.crs}')
    ##########################
    # Preparing `cat`, `riv`, and `nca` objects for the required watershed
    ## Preparing `MERIT-Basins` Layers
    # Hydrant's `geom` module provides the `prepare_cat`

    catchments = gm.prepare_cat(
        cat=cat, # 
        cat_col_id='COMID',
        cst=nca,
        cst_col_mapper={'FID':'COMID'},
        cst_col_id='COMID'
    )
    ##########################
    # prepare the `MERIT-Basins`'s river segments for the next
    # post-processing steps:

    rivers = gm.prepare_riv(
        riv=riv,
        riv_cols={
            'id':'COMID',
            'next_id':'NextDownID',
            'slope':'slope',
            'length':'lengthkm',
            'length_direct':'lengthdir'
        },
        cat=catchments,
        cat_cols={
            'id':'COMID',
            'hillslope':'hillslope',
            'geom':'geometry'
        }
    )
    return catchments, rivers, nca

#####################################
# extract and plots the geofabric upstream of an outlet point
def extract_geofabric(catchments, rivers, outlet_point):
    ## Subsetting Sub-basins and River Segments Upstream of the outlet_point
    # find the sub-basin that this gauge intersects with:
    # catchments[catchments.intersects(outlet_point)]
    # extract catchements and rivers

    extracted_catchments, extracted_rivers = gm.intersect_topology(
        cat=catchments,
        cat_cols={
            'id':'COMID'
        },
        riv=rivers,
        riv_cols={
            'id':'COMID',
            'next_id':'NextDownID'
        },
    outlet_id=catchments[catchments.intersects(outlet_point)]['COMID'].values)
    # plot what have extracted from the larger `MERIT-Basins` geospatial fabric:
    fig, ax = plt.subplots(
    nrows=1,
    ncols=1,
    figsize=(8, 8)
    )

    # sub-basins
    extracted_catchments.plot(ax=ax, color='gray', edgecolor='black', alpha=0.8, zorder=1)
    # river segments
    extracted_rivers.plot(ax=ax, color='blue', alpha=1, zorder=2)
    # gauge location
    ax.scatter(outlet_point.x, outlet_point.y, color='red', alpha=0.8, zorder=3)
    return extracted_catchments, extracted_rivers
#####################################
# save the extracted geofabric to file
def save_extracted_geofabric(extracted_catchments, extracted_rivers, output_path):
    # saving the results into the `output_path` directory

    # first, creating the directory
    try:
        os.makedirs(output_path)
    except FileExistsError:
        pass

    # then, saving the data
    extracted_catchments.to_file(os.path.join(output_path, 'extracted_subbasins.shp'))
    extracted_rivers.to_file(os.path.join(output_path, 'extracted_rivers.shp'))