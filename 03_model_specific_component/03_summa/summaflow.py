# summaflow: a collection of functions to create summa input files
# load the needed packages
# packages are loaded
import xarray as xr
import pint_xarray
import glob
import os
import pandas as pd
from   easymore import Easymore
import numpy as np
import geopandas   as      gpd


############################
# write summa forcing from MESH forcing
def write_summa_forcing(path_to_save, nc_file_source):
    if not os.path.isdir(path_to_save):
        os.makedirs(path_to_save)
    
    # rename the variables:
    ds = xr.open_dataset(nc_file_source)

    # time step of the data in seconds
    ds['data_step'] = (ds['time'][1].values - ds['time'][0].values).astype('timedelta64[s]') #.astype(int)

    # rename the dictionary
    rename_dictionary = {'subbasin': 'hruId',
                         'lat': 'latitude',
                         'lon': 'longitude',
                         'RDRS_v2.1_A_PR0_SFC': 'pptrate',
                         'RDRS_v2.1_P_TT_09944': 'airtemp',
                         'RDRS_v2.1_P_FB_SFC': 'SWRadAtm',
                         'RDRS_v2.1_P_FI_SFC': 'LWRadAtm',
                         'RDRS_v2.1_P_UVC_09944': 'windspd',
                         'RDRS_v2.1_P_HU_09944': 'spechum',
                         'RDRS_v2.1_P_P0_SFC': 'airpres'}
    ds = ds.rename_vars(rename_dictionary)
    ds = ds.rename({'subbasin': 'hru'})
    ds['hru'] = ds['hruId']

    ds = ds.transpose('time', 'hru') # to make sure the varibales are transposed
    
    # save file
    if os.path.isfile(path_to_save+'summa_forcing.nc'):
        os.remove(path_to_save+'summa_forcing.nc')

    ds.to_netcdf(path_to_save+'summa_forcing.nc')
    
    return(ds)
############################
def write_summa_attribute(path_to_save, subbasins_shapefile, rivers_shapefile, gistool_output):
    # read forcing file
    forcing = xr.open_dataset(os.path.join(path_to_save,'summa_forcing.nc'))

    # prepare data by merging the spatial fields into one dataframe
    # read merit geofabric
    # read rivers
    riv = gpd.read_file(rivers_shapefile)
    # reorder river to follow the same order as the hru in the forcing file
    riv = riv.set_index('COMID').loc[forcing['hru']].reset_index()
    
    # read catchments
    cat = gpd.read_file(subbasins_shapefile)
    # reorder river to follow the same order as the hru in the forcing file
    cat = cat.set_index('COMID').loc[forcing['hru']].reset_index()

    # read elevation stats
    elev_stats = pd.read_csv(os.path.join(gistool_output, 'modified_domain_stats_elv.csv'))
    # reorder df to follow the same order as the hru in the forcing file
    elev_stats = elev_stats.set_index('COMID').loc[forcing['hru']].reset_index()
    #rename columns except COMID
    elev_stats = elev_stats.rename(columns=lambda x: x + '_elev' if x != 'COMID' else x)

    # merge riv, cat, and elev_stat in one dataframe on COMID
    geofabric = riv.merge(cat, on='COMID')
    geofabric = geofabric.merge(elev_stats, on='COMID')
    geofabric = geofabric.drop(columns=['geometry_x', 'hillslope_x', 'hillslope_y', 'geometry_y'])

    # read soil and landcover data
    soil_stats = pd.read_csv(os.path.join(gistool_output, 'modified_domain_stats_soil_classes.csv'))
    # reorder df to follow the same order as the hru in the forcing file
    soil_stats = soil_stats.set_index('COMID').loc[forcing['hru']].reset_index()
    # rename majority column
    soil_stats = soil_stats.rename(columns={'majority': 'soil_majority'})
    
    landuse_stats = pd.read_csv(os.path.join(gistool_output, 'modified_domain_stats_NA_NALCMS_landcover_2020_30m.csv'))
    # reorder df to follow the same order as the hru in the forcing file
    landuse_stats = landuse_stats.set_index('COMID').loc[forcing['hru']].reset_index()
    # rename majority column
    landuse_stats = landuse_stats.rename(columns={'majority': 'landuse_majority'})
    
    soil_landuse_stats = landuse_stats.merge(soil_stats, on='COMID')

    ######---------------------------
    # link landuse and soil to vegtbl and soiltbl
    # link NALCMS to USGS ids in the VEGPARM.TBL
    # Manually matched NALCMS to USGS land cover types
    matched_landuse = [
        (0, 'Unknown', None, None),
        (1, 'Temperate or sub-polar needleleaf forest', 14, 'Evergreen Needleleaf Forest'),
        (2, 'Sub-polar taiga needleleaf forest', 12, 'Deciduous Needleleaf Forest'),
        (3, 'Tropical or sub-tropical broadleaf evergreen forest', 13, 'Evergreen Broadleaf Forest'),
        (4, 'Tropical or sub-tropical broadleaf deciduous forest', 11, 'Deciduous Broadleaf Forest'),
        (5, 'Temperate or sub-polar broadleaf deciduous forest', 11, 'Deciduous Broadleaf Forest'),
        (6, 'Mixed forest', 15, 'Mixed Forest'),
        (7, 'Tropical or sub-tropical shrubland', 8, 'Shrubland'),
        (8, 'Temperate or sub-polar shrubland', 8, 'Shrubland'),
        (9, 'Tropical or sub-tropical grassland', 10, 'Savanna'),
        (10, 'Temperate or sub-polar grassland', 7, 'Grassland'),
        (11, 'Sub-polar or polar shrubland-lichen-moss', 8, 'Shrubland'),
        (12, 'Sub-polar or polar grassland-lichen-moss', 7, 'Grassland'),
        (13, 'Sub-polar or polar barren-lichen-moss', 19, 'Barren or Sparsely Vegetated'),
        (14, 'Wetland', 17, 'Herbaceous Wetland'),
        (15, 'Cropland', 2, 'Dryland Cropland and Pasture'),
        (16, 'Barren lands', 19, 'Barren or Sparsely Vegetated'),
        (17, 'Urban', 1, 'Urban and Built-Up Land'),
        (18, 'Water', 16, 'Water Bodies'),
        (19, 'Snow and Ice', 24, 'Snow or Ice')
    ]
    
    # Create DataFrame
    matched_lanuse = pd.DataFrame(matched_landuse, columns=['NALCMS_ID', 'NALCMS_Description', 'USGS_ID', 'USGS_Description'])
    
    # Create a dictionary for mapping NALCMS IDs to USGS IDs
    nalcms_to_usgs = dict(zip(matched_lanuse['NALCMS_ID'], matched_lanuse['USGS_ID']))

    # link NALCMS to USGS ids in the VEGPARM.TBL
    # Manually matched soilgrid to ROSETTA types
    matched_soils = [
        (1, 'clay', 1, 'CLAY'),
        (2, 'silty clay', 10, 'SILTY CLAY'),
        (3, 'sandy clay', 6, 'SANDY CLAY'),
        (4, 'clay loam', 2, 'CLAY LOAM'),
        (5, 'silty clay loam', 11, 'SILTY CLAY LOAM'),
        (6, 'sandy clay loam', 7, 'SANDY CLAY LOAM'),
        (7, 'loam', 3, 'LOAM'),
        (8, 'silty loam', 12, 'SILT LOAM'),
        (9, 'sandy loam', 8, 'SANDY LOAM'),
        (10, 'silt', 9, 'SILT'),
        (11, 'loamy sand', 4, 'LOAMY SAND'),
        (12, 'sand', 5, 'SAND')
    ]
    
    # Create DataFrame
    matched_soils = pd.DataFrame(matched_soils, columns=['Soilgrid_ID', 'Soilgrid_Description', 'ROSETTA_ID', 'ROSETTA_Description'])
    # Create a dictionary for mapping Soilgrid IDs to ROSETTA IDs
    soilgrid_to_rosetta = dict(zip(matched_soils['Soilgrid_ID'], matched_soils['ROSETTA_ID']))
    ######---------------------------

    # Create a new xarray dataset
    attr = xr.Dataset()
    
    # prepare for the summa attr file
    attr ['hruId']          = xr.DataArray(geofabric['COMID'].values, dims=('hru'), 
                                           attrs={'long_name': 'Index of hydrological response units (HRU)', 'units': '-'})
    
    attr ['gruId']          = xr.DataArray(geofabric['COMID'].values, dims=('gru'),
                                          attrs={'long_name': 'Index of group of response unit (GRU)', 'units': '-'})
    
    attr ['hru2gruId']      = xr.DataArray(geofabric['COMID'].values, dims=('hru'),
                                          attrs={'long_name': 'Index of GRU to which the HRU belongs', 'units': '-'})
    
    attr ['downHRUindex']   = xr.DataArray(np.zeros(len(geofabric['COMID'])), dims=('hru'),
                                          attrs={'long_name': 'Index of downslope HRU (0 = basin outlet)', 'units': '-'}).astype(int)
    
    attr ['elevation']      = xr.DataArray(geofabric['mean_elev'].values, dims=('hru'),
                                          attrs={'long_name': 'Elevation of HRU\'s centriod point', 'units': 'm'})
    
    attr ['HRUarea']        = xr.DataArray(geofabric['unitarea'].values, dims=('hru'),
                                          attrs={'long_name': 'Area of each HRU', 'units': 'm^2'})
    
    attr ['tan_slope']      = xr.DataArray(geofabric['slope'].values, dims=('hru'),
                                          attrs={'long_name': 'Average tangent slope of HRU', 'units': 'm/m'})
    
    attr ['contourLength']  = xr.DataArray(geofabric['lengthkm'].values*1000, dims=('hru'),
                                            attrs={'long_name': 'ContourLength of HRU', 'units': 'm'})
    
    attr ['slopeTypeIndex'] = xr.DataArray(np.ones(len(geofabric['COMID'])), dims=('hru'),
                                           attrs={'long_name': 'Index defining slope', 'units': '-'}).astype(int)
    
    attr ['soilTypeIndex']  = xr.DataArray(list(map(soilgrid_to_rosetta.get, soil_landuse_stats['soil_majority'].values)), dims=('hru'),
                                            attrs={'long_name': 'Index defining soil type - ROSETTA', 'units': '-'}).astype(int)
    
    attr ['vegTypeIndex']   = xr.DataArray(list(map(nalcms_to_usgs.get, soil_landuse_stats['landuse_majority'].values)), dims=('hru'),
                                           attrs={'long_name': 'Index defining vegetation type - USGS', 'units': '-'}).astype(int)
    
    attr ['mHeight']        = xr.DataArray(np.ones(len(geofabric['COMID']))*40, dims=('hru'),
                                          attrs={'long_name': 'Measurement height above bare ground', 'units': 'm'})
    
    # write attribute file
    if os.path.isfile(path_to_save+'summa_zLocalAttributes.nc'):
        os.remove(path_to_save+'summa_zLocalAttributes.nc')
    
    attr.to_netcdf(path_to_save+'summa_zLocalAttributes.nc')
    # return the attribute file
    return(attr)
#################################
def write_summa_paramtrial(attr, path_to_save):
    
    # Define dimensions
    hru_size = len(attr['hruId'].values)
    
    # Create a new xarray dataset
    ds = xr.Dataset()
    
    # Add dimensions to the dataset
    ds['hru'] = xr.DataArray(attr['hruId'].values, dims=('hru'), attrs={'units': '-'})
    
    # Add variables to the dataset
    ds['hruId'] = xr.DataArray(attr['hruId'].values, dims=('hru'), attrs={'units': '-', 'long_name': 'Index of hydrological response unit (HRU)'})
    ds
    if os.path.isfile(path_to_save+'summa_zParamTrial.nc'):
        os.remove(path_to_save+'summa_zParamTrial.nc')
    
    ds.to_netcdf(path_to_save+'summa_zParamTrial.nc')
    return(ds)

#################################
def write_summa_initial_conditions(attr, soil_mLayerDepth, path_to_save):
    
    # Define dimensions
    hru_size = len(attr['hruId'].values)
    scalarv_size = 1
    midSoil_size = midToto_size = len(soil_mLayerDepth)
    ifcToto_size = midSoil_size + 1

    # Create a new xarray dataset
    ds = xr.Dataset()
    
    # Add dimensions to the dataset
    ds['hru'] = xr.DataArray(attr['hruId'].values, dims=('hru'), attrs={'units': '-'})
    ds['midSoil'] = xr.DataArray(range(midSoil_size), dims=('midSoil'))
    ds['midToto'] = xr.DataArray(range(midToto_size), dims=('midToto'))
    ds['ifcToto'] = xr.DataArray(range(ifcToto_size), dims=('ifcToto'))
    ds['scalarv'] = xr.DataArray(range(scalarv_size), dims=('scalarv'))
    
    # Add variables to the dataset
    ds['hruId'] = xr.DataArray(attr['hruId'].values, dims=('hru'), attrs={'units': '-', 'long_name': 'Index of hydrological response unit (HRU)'})
    ds['dt_init'] = xr.DataArray([[3600.0] * hru_size], dims=('scalarv', 'hru'))
    ds['nSoil'] = xr.DataArray([[midSoil_size] * hru_size], dims=('scalarv', 'hru'))
    ds['nSnow'] = xr.DataArray([[0] * hru_size], dims=('scalarv', 'hru'))
    ds['scalarCanopyIce'] = xr.DataArray([[0] * hru_size], dims=('scalarv', 'hru'))
    ds['scalarCanopyLiq'] = xr.DataArray([[0] * hru_size], dims=('scalarv', 'hru'))
    ds['scalarSnowDepth'] = xr.DataArray([[0] * hru_size], dims=('scalarv', 'hru'))
    ds['scalarSWE'] = xr.DataArray([[0] * hru_size], dims=('scalarv', 'hru'))
    ds['scalarSfcMeltPond'] = xr.DataArray([[0] * hru_size], dims=('scalarv', 'hru'))
    ds['scalarAquiferStorage'] = xr.DataArray([[1.0] * hru_size], dims=('scalarv', 'hru'))
    ds['scalarSnowAlbedo'] = xr.DataArray([[0] * hru_size], dims=('scalarv', 'hru'))
    ds['scalarCanairTemp'] = xr.DataArray([[283.16] * hru_size], dims=('scalarv', 'hru'))
    ds['scalarCanopyTemp'] = xr.DataArray([[283.16] * hru_size], dims=('scalarv', 'hru'))
    ds['mLayerTemp'] = xr.DataArray([[283.16] * hru_size ] * midToto_size , dims=('midToto', 'hru'))
    ds['mLayerVolFracIce'] = xr.DataArray([[0] * hru_size] * midToto_size, dims=('midToto', 'hru'))
    ds['mLayerVolFracLiq'] = xr.DataArray([[0.2] * hru_size] * midToto_size, dims=('midToto', 'hru'))
    ds['mLayerMatricHead'] = xr.DataArray([[-1] * hru_size] * midToto_size, dims=('midSoil', 'hru'))
    ds['iLayerHeight'] = xr.DataArray(np.transpose([[0]+list(np.cumsum(soil_mLayerDepth))] * hru_size) , dims=('ifcToto', 'hru'))
    ds['mLayerDepth'] = xr.DataArray(np.transpose([soil_mLayerDepth] * hru_size), dims=('midToto',  'hru'))
    
    
    if os.path.isfile(path_to_save+'summa_zInitialCond.nc'):
        os.remove(path_to_save+'summa_zInitialCond.nc')
    
    ds.to_netcdf(path_to_save+'summa_zInitialCond.nc')
    
    return(ds)

##################################
def write_summa_filemanager(path_to_save, forcing):


    start_date = pd.to_datetime(forcing['time'][0].values).strftime('%Y-%m-%d %H:%M')
    end_date = pd.to_datetime(forcing['time'][-1].values).strftime('%Y-%m-%d %H:%M')

    template_string="""controlVersion       'summa_FILE_MANAGER_V3.0.0' !  fman_ver 
simStartTime         '{start_date}' ! 
simEndTime           '{end_date}' ! 
tmZoneInfo           'localTime' ! 
settingsPath         './' !  setting_path
forcingPath          './' !  input_path
outputPath           './results/' !  output_path
decisionsFile        './summa_zDecisions_celia1990.txt' !  decision
outputControlFile    './Model_Output.txt' !  OUTPUT_CONTROL
globalHruParamFile   './summa_zLocalParamInfo.txt' !  local_par
globalGruParamFile   './summa_zBasinParamInfo.txt' !  basin_par
attributeFile        './summa_zLocalAttributes.nc' !  local_attr
trialParamFile       './summa_zParamTrial.nc' !  para_trial
forcingListFile      './summa_zForcingFileList.txt' !  forcing_list
initConditionFile    './summa_zInitialCond.nc' !  initial_cond
outFilePrefix        'myTest' !  output_prefix
vegTableFile         'VEGPARM.TBL' ! 
soilTableFile        'SOILPARM.TBL' ! 
generalTableFile     'GENPARM.TBL' ! 
noahmpTableFile      'MPTABLE.TBL' ! 
! history written by summaflow
"""

    # Replace the placeholder with the actual start date using an f-string
    fileManager = template_string.format(start_date=start_date, end_date = end_date)


    # Optionally, write the formatted string to a file
    output_file = path_to_save+'summa_fileManager.txt'
    with open(output_file, 'w') as file:
        file.write(fileManager)

