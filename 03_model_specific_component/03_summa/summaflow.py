# summaflow: a collection of functions to create SUMMA input files
# load the needed packages
# packages are loaded
import xarray as xr
import pint_xarray
import glob
import netCDF4 as nc4
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
    ds['data_step'] = (ds['time'][1].values - ds['time'][0].values).astype('timedelta64[s]').astype(int)

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
    if os.path.isfile(path_to_save+'SUMMA_forcing.nc'):
        os.remove(path_to_save+'SUMMA_forcing.nc')

    ds.to_netcdf(path_to_save+'SUMMA_forcing.nc')
    
    # replace T in the time unit with space

    ncid = nc4.Dataset(path_to_save + 'SUMMA_forcing.nc', 'r+')

    # Access the 'units' attribute and replace 'T' with a space
    units_attribute = ncid['time'].units
    units_attribute = units_attribute.replace('T', ' ')

    # Update the 'units' attribute in the netCDF file
    ncid['time'].setncattr('units', units_attribute)

    # Close the netCDF file
    ncid.close()
    
############################
def write_summa_attribute():
    