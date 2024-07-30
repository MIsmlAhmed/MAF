# Load needed packages
import xarray as xr
import pint_xarray
import glob
import netCDF4 as nc4
import os
import cdo
import pandas as pd
from   easymore import Easymore
import numpy       as      np
import geopandas   as      gpd
import sys
from   itertools   import  product
import datetime



# sort geodata from upstream to downstream
def sort_geodata(geodata):
    # find the subbasin order
    geodata['n_ds_subbasins'] = 0
    for index, row in geodata.iterrows():
        n_ds_subbasins = 0
        nextID = row['maindown']
        nextID_index = np.where(geodata['subid']==nextID)[0]

        while (len(nextID_index>0) or nextID > 0) :
            n_ds_subbasins += 1
            nextID = geodata['maindown'][nextID_index[0]]
            nextID_index = np.where(geodata['subid']==nextID)[0]
            
        geodata.loc[index, 'n_ds_subbasins']= n_ds_subbasins

    # Sort the DataFrame by 'n_ds_subbasins' in descending order
    geodata = geodata.sort_values(by='n_ds_subbasins', ascending=False, ignore_index=True)
    geodata = geodata.drop(columns=['n_ds_subbasins'])
    return geodata

# write HYPE forcing from MESH nc file

def write_hype_forcing(easymore_output, path_to_save):
    if not os.path.isdir(path_to_save):
        os.makedirs(path_to_save)

    def convert_hourly_to_daily (input_file_name,
                                 variable_in,
                                 variable_out,
                                 variable_out_long_name = None,
                                 var_unit_conversion = None,
                                 var_time = 'time',
                                 var_id = 'id',
                                 time_diff = 0,
                                 stat = 'max', 
                                 output_file_name_nc = None,
                                 output_file_name_txt = None,
                                 Fill_value = -9999.0): # 'max', 'min', 'mean'

        # read the input houtly nc file
        ds = xr.open_dataset(input_file_name)
        # set id as integer
        ds.coords[var_id] = ds.coords[var_id].astype(int)
        # ds = ds.rename({'subbasin': 'id'})

        # drop all the other variables except the mentioned varibale, time and id
        variables_to_keep = [variable_in, var_time]
        if not var_id is None:
            variables_to_keep.append(var_id)

        # Drop all variables except the specified ones
        ds = ds.drop([v for v in ds.variables if v not in variables_to_keep])

        # roll the time based on hour of difference to have more accurate
        if time_diff !=0:
            ds[var_time] = ds[var_time].roll(time=time_diff)
            # Remove the first or last roll_steps time steps
            if time_diff < 0:
                ds = ds.isel(time=slice( None, time_diff))
            elif time_diff > 0:
                ds = ds.isel(time=slice( time_diff, None))

        # to create the xarray dataframe with daily time
        if stat == 'max':
            ds_daily = ds.resample(time='D').max()
        elif stat == 'min':
            ds_daily = ds.resample(time='D').min()
        elif stat == 'mean':
            ds_daily = ds.resample(time='D').mean()
        elif stat == 'sum':
            ds_daily = ds.resample(time='D').sum()
        else:
            sys.exit('input stat should be max, min, mean or sum')

        # conversion of units based on provided conversion unit
        ds_daily[variable_in] = ds_daily[variable_in].pint.quantify(var_unit_conversion['in_unit'])
        ds_daily[variable_in] = ds_daily[variable_in].pint.to(var_unit_conversion['out_unit'])
        ds_daily = ds_daily.pint.dequantify()

        # drop the vairiable in
        ds_daily = ds_daily.rename({variable_in: variable_out})

        # add long name
        if not variable_out_long_name is None:
            ds_daily[variable_out].attrs['long_name'] = variable_out_long_name

        # transpose the variable
        ds_daily[variable_out] = ds_daily[variable_out].transpose()

        # this section is written to avoid issues with netcdf and HYPE!
        # I could not find what is the issue, however, when the data is 
        # transferred to df, tranfer back to xarray and saved, the issue
        # with HYPE is resolved. this need closer look. Also HYPE netcdf
        # is in its initial stage of developement and can have issue as
        # well
        df = ds_daily[variable_out].to_dataframe()
        df = df.unstack()
        df = df.T
        df = df.droplevel(level=0, axis=0)
        df.columns.name = None
        df.index.name = var_time
        if not output_file_name_txt is None:
            df.to_csv(output_file_name_txt,\
                      sep='\t', na_rep='', index_label='time', float_format='%.3f')
        esmr = Easymore()
        ds_daily = esmr.dataframe_to_netcdf_xr(df,
                                         data_frame_DateTime_column = var_time,
                                         variable_name = variable_out,
                                         variable_dim_name = 'id',
                                         unit_of_variable = var_unit_conversion['out_unit'],
                                         variable_long_name = variable_out_long_name,
                                         Fill_value = Fill_value)

        # save the file if path is provided
        if not output_file_name_nc is None:
            if os.path.isfile(output_file_name_nc):
                os.remove(output_file_name_nc)
            ds_daily.to_netcdf(output_file_name_nc,\
                               encoding = {variable_out:{'_FillValue':Fill_value}})

        # return
        return ds_daily
    ############
    print('Merging easymore outputs to one NetCDF file \n')
    # Replace with your file path pattern
    easymore_nc_files = sorted(glob.glob(easymore_output+'/*.nc'))
    # split the files in batches as cdo cannot mergetime long list of file names
    batch_size = 15
    # avoid splitting files if their number is too small
    if(len(easymore_nc_files) < batch_size):
        batch_size = len(easymore_nc_files)
    files_split = np.array_split(easymore_nc_files,batch_size)
    cdo_obj = cdo.Cdo()  # CDO object
    intermediate_files = []

    # split files into intermediate files
    # Combine in batches
    for i in range(batch_size):
        print(f'Processing easymore outputs: batch no {i+1} out of {batch_size} \n')
        batch_files = files_split[i].tolist()
        batch_output = f'forcing_batch_{i}.nc'
        cdo_obj.mergetime(input=batch_files, output=batch_output)
        intermediate_files.append(batch_output)

    # Combine intermediate results into one netcdf file
    cdo_obj.mergetime(input=intermediate_files, output='RDRS_forcing.nc')

    # Clean up intermediate files if needed
    for f in intermediate_files:
        os.remove(f)

    # open the forcing file
    forcing = xr.open_dataset('RDRS_forcing.nc')
    # convert calendar to 'standard'
    forcing = forcing.convert_calendar('standard')
    # write to netcdf
    forcing.to_netcdf('RDRS_forcing.nc')
    forcing.close()
    ############
    
    ds1= convert_hourly_to_daily('RDRS_forcing.nc',
                                'RDRS_v2.1_P_TT_09944',
                                'TMAXobs',
                                var_unit_conversion = {'in_unit':'degreeC','out_unit':'degreeC'},
                                var_time = 'time',
                                var_id = 'COMID',
                                time_diff = -7,
                                stat = 'max',
                                # output_file_name_nc = path_to_save+'TMAXobs.nc',
                                output_file_name_txt = path_to_save+'TMAXobs.txt')

    ds2= convert_hourly_to_daily('RDRS_forcing.nc',
                                'RDRS_v2.1_P_TT_09944',
                                'TMINobs',
                                var_unit_conversion = {'in_unit':'degreeC','out_unit':'degreeC'},
                                var_time = 'time',
                                var_id = 'COMID',
                                time_diff = -7,
                                stat = 'min',
                                # output_file_name_nc = path_to_save+'TMINobs.nc',
                                output_file_name_txt = path_to_save+'TMINobs.txt')

    ds3= convert_hourly_to_daily('RDRS_forcing.nc',
                                'RDRS_v2.1_P_TT_09944',
                                'Tobs',
                                var_unit_conversion = {'in_unit':'degreeC','out_unit':'degreeC'},
                                var_time = 'time',
                                var_id = 'COMID',
                                time_diff = -7,
                                stat = 'mean',
                                # output_file_name_nc = path_to_save+'Tobs.nc',
                                output_file_name_txt = path_to_save+'Tobs.txt')

    ds4= convert_hourly_to_daily('RDRS_forcing.nc',
                                'RDRS_v2.1_A_PR0_SFC',
                                'Pobs',
                                var_unit_conversion = {'in_unit':'m h**-1',\
                                                       'out_unit':'mm day**-1'},
                                var_time = 'time',
                                var_id = 'COMID',
                                time_diff = -7,
                                stat = 'mean',
                                # output_file_name_nc = path_to_save+'Pobs.nc',
                                output_file_name_txt = path_to_save+'Pobs.txt')
    
    # remove the merged netcdf file
    os.remove('RDRS_forcing.nc')

################################################################
# write GeoData and GeoClass files
def write_hype_geo_files(gistool_output, subbasins_shapefile, rivers_shapefile, frac_threshold, path_to_save):
    
    if not os.path.isdir(path_to_save):
        os.makedirs(path_to_save)
        
    # load the information from the gistool for soil and land cover and find the number of geoclass
    soil_type = pd.read_csv(gistool_output+'modified_domain_stats_soil_classes.csv')
    landcover_type = pd.read_csv(gistool_output+'modified_domain_stats_NA_NALCMS_landcover_2020_30m.csv')
    elevation_mean = pd.read_csv(gistool_output+'modified_domain_stats_elv.csv')

    soil_type = soil_type.sort_values(by='COMID').reset_index(drop=True)
    landcover_type = landcover_type.sort_values(by='COMID').reset_index(drop=True)
    elevation_mean = elevation_mean.sort_values(by='COMID').reset_index(drop=True)
    
    # find the combination of the majority soil and land cover
    combinations_set_all = set()
    for index, row in landcover_type.iterrows():
        # get the fraction for land cover for each row
        fractions = [col for col in landcover_type.columns if col.startswith('frac') and row[col] > frac_threshold]
        # remove frac_ from the list
        fractions = [col.split('_')[1] for col in fractions]
        fractions = [int(name) for name in fractions]

        # get the majority soil type for each row
        majority_soil = [soil_type['majority'].iloc[index].item()]

        # Combine as combination of soil and land cover and keep as a set
        combinations = list(product(fractions, majority_soil))
        combinations_set = set(combinations)
        combinations_set_all.update(combinations_set)

    # print(combinations_set_all)
    # print(len(combinations_set_all))

    data_list = [{'landcover': item[0], 'soil': item[1]} for item in combinations_set_all]

    # Create a pandas DataFrame from the list of dictionaries
    combination = pd.DataFrame(data_list)

    combination ['SLC'] = 0
    combination ['SLC'] = np.arange(len(combination))+1
    # renumber landcover and soil sequentially
    
    #######################
    landcover_type_prepared = landcover_type.copy()

    for i in range(1, len(combination)+1):
        column_name = f'SLC_{i}'
        landcover_type_prepared[column_name] = 0.00

    landcover_type_prepared['soil'] = soil_type['majority']

    def get_non_zero_columns(row):
        return [col for col in row.index if col.startswith('frac_') and row[col] > frac_threshold]

    # Apply the function to each row
    landcover_type_prepared['non_zero_columns'] = landcover_type_prepared.apply(get_non_zero_columns, axis=1)


    for index, row in landcover_type_prepared.iterrows():
        # get the soil type
        soil_type_value = soil_type['majority'].iloc[index]

        for i in row['non_zero_columns']:

            # remove frac from column name 
            land_cover_value = i.replace("frac_", "")

            # get the SLC value
            result = combination[(combination['landcover'] == int(land_cover_value)) & (combination['soil'] == int(soil_type_value))]['SLC']
            column_name = 'SLC_'+str(result.values[0])
            landcover_type_prepared.loc[index, column_name] = landcover_type_prepared[i].iloc[index]
#######################
    riv = gpd.read_file(rivers_shapefile)
    riv.sort_values(by='COMID').reset_index(drop=True)
    riv['lengthm'] = 0.00
    riv['lengthm'] = riv['lengthkm'] * 1000

    cat = gpd.read_file(subbasins_shapefile)
    cat.sort_values(by='COMID').reset_index(drop=True)
    cat['area'] = 0.00
    cat['area'] = cat['unitarea'] * 1000000 # km2 to m2
    cat['latitude'] = cat.centroid.y
    cat['longitude'] = cat.centroid.x
    
    # add information to the geodata dataframe
    landcover_type_prepared['NextDownID'] = riv['NextDownID']
    landcover_type_prepared['area'] = cat['area']
    landcover_type_prepared['latitude'] = cat['latitude']
    landcover_type_prepared['longitude'] = cat['longitude']
    landcover_type_prepared['elev_mean'] = elevation_mean['mean']
    landcover_type_prepared['slope_mean'] = riv['slope']
    landcover_type_prepared['rivlen'] = riv['lengthm']
    landcover_type_prepared['uparea'] = riv['uparea']
    
    column_name_mapping = {
    'COMID': 'subid',
    'NextDownID': 'maindown',
    'area': 'area',
    'latitude': 'latitude',
    'longitude': 'longitude',
    'elev_mean': 'elev_mean',
    'slope_mean': 'slope_mean',
    'rivlen': 'rivlen'
    }

    # Rename the columns based on the dictionary
    landcover_type_prepared = landcover_type_prepared.rename(columns=column_name_mapping)
    # landcover_type_prepared

    #
    slc_columns = [col for col in landcover_type_prepared.columns if col.startswith('SLC_')]

    # Sort the columns as per your requirements
    column_order = ['subid', 'maindown', 'area', 'latitude', 'longitude', 'elev_mean', 'slope_mean', 'rivlen'] + slc_columns #+ ['uparea']

    landcover_type_prepared = landcover_type_prepared[column_order]
    # landcover_type_prepared
    #######################
    # sort geodata file from upstream to downstream
    landcover_type_prepared = sort_geodata(landcover_type_prepared)
    
    # normalize fracs
    # Identify columns starting with 'SLC_'
    slc_columns = [col for col in landcover_type_prepared.columns if col.startswith('SLC_')]

    # Normalize SLC values so that they sum to 1 for each row
    landcover_type_prepared[slc_columns] = landcover_type_prepared[slc_columns].div(landcover_type_prepared[slc_columns].sum(axis=1), axis=0)

    # landcover_type_prepared = landcover_type_prepared.drop(columns=['uparea'])
    landcover_type_prepared.to_csv(path_to_save+'GeoData.txt', sep='\t', index=False)
    # landcover_type_prepared
    #######################
    
    # write geoclass file

    combination = combination.rename(columns={'landcover': 'LULC'})
    combination = combination.rename(columns={'soil': 'SOIL TYPE'})
    combination = combination[['SLC','LULC','SOIL TYPE']]
    combination['Main crop cropid'] = 0
    combination['Second crop cropid'] = 0
    combination['Crop rotation group'] = 0
    combination['Vegetation type'] = 1
    combination['Special class code'] = 0
    combination['Tile depth'] = 0
    combination['Stream depth'] = 2.296
    combination['Number of soil layers'] = 3
    combination['Soil layer depth 1'] = 0.091
    combination['Soil layer depth 2'] = 0.493
    combination['Soil layer depth 3'] = 2.296

    # combination
    #######################
    # Add commented lines
    commented_lines = [
    """! MODIS landcover													
! Add legend (raster value) and discription													
!	original legend (raster_value)	description											
!   1: 'Temperate or sub-polar needleleaf forest',
!   2: 'Sub-polar taiga needleleaf forest',
!   3: 'Tropical or sub-tropical broadleaf evergreen forest',
!   4: 'Tropical or sub-tropical broadleaf deciduous forest',
!   5: 'Temperate or sub-polar broadleaf deciduous forest',
!   6: 'Mixed forest',
!   7: 'Tropical or sub-tropical shrubland',
!   8: 'Temperate or sub-polar shrubland',
!   9: 'Tropical or sub-tropical grassland',
!   10: 'Temperate or sub-polar grassland',
!   11: 'Sub-polar or polar shrubland-lichen-moss',
!   12: 'Sub-polar or polar grassland-lichen-moss',
!   13: 'Sub-polar or polar barren-lichen-moss',
!   14: 'Wetland',
!   15: 'Cropland',
!   16: 'Barren lands',
!   17: 'Urban',
!   18: 'Water',
!   19: 'Snow and Ice',											
!													
!													
!													
! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------													
!	SoilGrid V1												
!		original legend (raster_value)	description										
!	 C 	    1	 clay										
!	 SIC 	2	 silty clay										
!	 SC 	3	 sandy clay										
!	 CL 	4	 clay loam										
!	 SICL 	5	 silty clay loam										
!	 SCL 	6	 sandy clay loam										
!	 L   	7	 loam										
!	 SIL 	8	 silty loam										
!	 SL 	9	 sandy loam										
!	 SI 	10	 silt										
!	 LS 	11	 loamy sand										
!	 S  	12	 sand										
!													
!													
! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	"""
    ]

    # Open the file in write mode
    with open(path_to_save+'GeoClass.txt', 'w') as file:
        # Write the commented lines
        for line in commented_lines:
            file.write(line + '\n')

    # re-number landcover and soil   
    for i in ['LULC', 'SOIL TYPE']:
        # Identify unique values and create a mapping from old values to new sequential values
        unique_values = pd.unique(combination[i])
        value_mapping = {value: idx + 1 for idx, value in enumerate(unique_values)}

        # Apply the new numbering to the LULC column
        combination[i] = combination[i].map(value_mapping)

        # Create a DataFrame for the mapping and write to geoclass
        mapping_df = pd.DataFrame(list(value_mapping.items()), columns=['Old Value', 'New Value'])

        with open(path_to_save+'GeoClass.txt', 'a') as f:
            f.write('! changes (reclassification) to '+i+'\n')
            for _, row in mapping_df.iterrows():
                # Write each row to the file with "!" at the beginning
                f.write(f"! {row['Old Value']} -> {row['New Value']}\n")


    # writing the `GeoClass.txt` file
    with open(path_to_save+'GeoClass.txt', 'a') as file:
            file.write("""! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
!          SLC	LULC	SOIL TYPE	Main crop cropid	Second crop cropid	Crop rotation group	Vegetation type	Special class code	Tile depth	Stream depth	Number of soil layers	Soil layer depth 1	Soil layer depth 2	Soil layer depth 3 \n""")
            combination.to_csv(file, sep='\t', index=False, header=False)
            

################################################################

# write par.txt file
def write_hype_par_file(path_to_save):

    output_file = path_to_save+'par.txt'

    if os.path.isfile(output_file):
        os.remove(output_file)

    # write out first text section
    header= {
        'section_head' :"""!!	=======================================================================================================																
!! Parameter file for: Case Name																	
!!	------------------------																
!!																	
!!	=======================================================================================================																"""
    }

    General_meteo_param = {
        'section_head' :"""!! METEOROLOGICAL INPUT - general parameters related to temperature and precipitation corrections																	
!!	-----																
!! All of these will be kept as 0 because we are not correcting the temperature or the precpitation																""",
        'tcobselev'    :{'value': 0, 'comment': "!! parameter for temperature correction due to observation elevation deviation from subbasin elevation (deg C)"},
        'tcalt'        :{'value': 0, 'comment': "!! parameter for temperature’s elevation dependence"},
        'tcelevadd'    :{'value': 0, 'comment': "!! parameter for temperature’s elevation dependence"},
        'pcaddg'       :{'value': 0, 'comment': "!! correction parameter for precipitation"},
        'pcusnow'      :{'value': 0, 'comment': "!! undercatch correction for snowfall"},
        'pcurain'      :{'value': 0, 'comment': "!! undercatch correction for rainfall"}
    }

    Snow_param = {
        'section_head' :"""!!	=======================================================================================================																
!!	"SNOWMELT, DISTRIBUTION, DENSITY, HEAT; sublimation is sorted under Evapotranspiration"																
!!	-----																
!!	"General snow accumulation and melt related parameters (baseline values from SHYPE, unless noted otherwise)"																
!! Snow distribution submodel: 0, snow melt submodel: 2, snow density submodel: 0, snow heat submodel: 1 (including snkika below)""",
        'ttpi'         :{'value': 1.7083, 'comment': "!! half of temperature interval with mixed snow and rainfall"},
        'sdnsnew'      :{'value': 0.13,   'comment': "!! density of new-fallen snow (kg/dm3)"},
        'snowdensdt'   :{'value': 0.0016, 'comment': "!! increase of snow density per day"},
        'fsceff'       :{'value': 0.99,   'comment': "!! efficiency of snow cover to influence snow melt and snow evaporation"},
        'cmrefr'       :{'value': 0.05,   'comment': "!! refreeze efficiency compared to the degree-day snow melt factor Used for second snow melt model"},
        'whcsnow'      :{'value': 0.08,   'comment': "!! water holding capacity of snow"}
    }

    # for 19 land covers

    Snow_land_submodel_1 = {
        'section_head' :"""!!	-----																
!!	SNOWMELT Landuse dependent parameters													
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19
!! snowmelt submodel:2, snow heat submodel: 1""",
        'ttmp'         :{'value'  : [-9.7740,-2.4419,2.5478,-3.8724,2.9143,-7.2759,-6.1012,-6.5266,\
                                     -1.8872,-1.2143,-9.9603,-5.4364,-9.774,-9.774,-9.7740,-9.7740,\
                                     -9.7740,-9.7740,-9.7740],
                         'comment': "!! threshold temperature for snow melt snow density and evapotranspiration"},
        'cmlt'         :{'value'  : [9.7021,6.0035,1.1786,9.3525,1.7176,5.8523,4.1957,8.6383,\
                                     8.0090,5.4865,1.1010,5.5150,9.7021,9.7021,9.7021,9.7021,\
                                     9.7021,9.7021,9.7021],
                         'comment': "!! melting parameter for snow"},
        'cmrad'        :{'value'  : [0.249065876,0.249065876,1.5,0.176534861,0.685361445,0.174564317,\
                                     0.174564317,0.174564317,0.174564317,0.685361445,0.501842737,\
                                     0.011482887,0.249065876,0.249065876,0.249065876,0.249065876,\
                                     0.249065876,0.249065876,0.249065876],
                         'comment': "!! coefficient for radiation snow melt, parameter for second snowmelt model"},
        'snalbmin'     :{'value'  : [0.524781764,0.524781764,0.45,0.250044137,0.243925437,0.251664609,\
                                    0.251664609,0.251664609,0.251664609,0.243925437,0.409460604,0.22856541,\
                                    0.524781764,0.524781764,0.524781764,0.524781764,0.524781764,\
                                    0.524781764,0.524781764],
                         'comment': "!! parameter for second snowmelt model"},
        'snalbmax'     :{'value'  : [0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,\
                                    0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,\
                                    0.85,0.85,0.85],
                         'comment': "!! parameter for second snowmelt model"},
        'snalbkexp'    :{'value'  : [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,\
                                    0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1],
                         'comment': "!! parameter for second snowmelt model"},
        'snkika'       :{'value'  : [50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50,50],
                         'comment': "!! snow heat model, relation between snow thermal conductivity and surface heat exchange coefficient"}
    }

    # for 19 land covers

    Snow_land_submodel_2 = {
        'section_head' :"""!!	-----																
!!	SNOWCOVER parameters (general and landuse) - baseline from Rossby RCA model (Samuelsson et al 2006;Lindström et al)																
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19
!! used in SNOWMELT submodel:2 """,
        'fscdistmax'   :{'value'  : [0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,\
                                     0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8],
                         'comment': "!! maximum snow distribution factor"},
        'fscdist0'     :{'value'  : [0.571998656,0.571998656,0.6,0.672227979,0.718591213,\
                                     0.672161579,0.672161579,0.672161579,0.672161579,\
                                     0.718591213,0.302164137,0.663832068,0.663832068,\
                                     0.663832068,0.663832068,0.663832068,0.663832068,\
                                     0.663832068,0.663832068],
                         'comment': "!! minimum snow distribution factor"},
        'fscdist1'     :{'value'  : [0.001,0.001,0.001,0.001,0.001,0.001,0.001,\
                                     0.001,0.001,0.001,0.001,0.001,0.001,0.001,\
                                     0.001,0.001,0.001,0.001,0.001],
                         'comment': "!! std coefficient for snow distribution factor parameter for second snowmelt model"},
        'fscmax'       :{'value'  : 1.00     , 'comment': "!! maximum fractional snow cover area"},
        'fscmin'       :{'value'  : 0.01     , 'comment': "!! minimum fractional snow cover area"},
        'fsclim'       :{'value'  : 0.001    , 'comment': "!! limit of fractional snow cover area for onset of snowmax"},
        'fsck1'        :{'value'  : 0.2      , 'comment': "!! Snowmass threshold to start decreasing the internal snowmax variable towards the end of the melt season"},
        'fsckexp'      :{'value'  : 0.000001 , 'comment': "!! Coefficient in exponential decrease of the internal Snowmax variable"}
    }

    # Glacier

    Glacier  = {
        'section_head' :"""!!	=======================================================================================================																
!!	"GLACIER - parameters for volume-area scaling, accumulation and melt (sublimation sorted under Evapotranspiration)"																
!!	-----																
!!	Glacier volume-area scaling	
!! the parameters used calculate the area of the glaciers""",
        'glacvexp'     :{'value'  : 1.395    , 'comment': "!! exponent of glacier area-volume relationship for glacier of type zero"},
        'glacvcoef'    :{'value'  : 0.17157  , 'comment': "!! coefficient of glacier area-volume relationship for glacier of type zero"},
        'glacvexp1'    :{'value'  : 1.25     , 'comment': "!! exponent of glacier area-volume relationship for glacier of type one)"},
        'glacvcoef1'   :{'value'  : 2.88364  , 'comment': "!! coefficient of glacier area-volume relationship for glacier of type one"},
        'glac2arlim'   :{'value'  : 25000000 , 'comment': "!! area limit for determine glacier type which is used only if glacier type is given in GlacierData.txt"},
        'glacannmb'    :{'value'  : 0        , 'comment': "!! annual mass balance for correction of initial glacier volume"}
    }

    # Glacier melt

    Glacier_melt = {
        'section_head' :"""!!	-----																
!!	Glacier melt parameters 																
!!	----																
!! considered with snowevaporation submodel: 1, snowmelt submodel 2""",
        'glacttmp'     :{'value'  : 0          , 'comment': "!! threshold temperature for glacier melt"},
        'glaccmlt'     :{'value'  : 1.58595482 , 'comment': "!! melting parameter for glacier"},
        'glaccmrad'    :{'value'  : 0.19090136 , 'comment': "!! coefficient for radiation glacier melt parameter for second snowmelt model"},
        'glaccmrefr'   :{'value'  : 0.06259448 , 'comment': "!! refreeze efficiency compared to the degree-day glacier melt factor parameter for second snow meltmodel"},
        'glacalb'      :{'value'  : 0.35       , 'comment': "!! albedo for glacier ice"},
        'fepotglac'    :{'value'  : 0          , 'comment': "!! fraction of snow-free potential evapotranspiration for first snowevaporation model"}
    } 

    # evaporation

    Evap = {
        'section_head' :"""!!	=======================================================================================================																
!!	EVAPOTRANSPIRATION PARAMETERS																
!!	-----																
!!	General evapotranspiration parameters																
!! used for petmodel""",
        'lp'           :{'value'  : 0.546 , 'comment': "!! Threshold for water content reduction of transpiration as fraction of field capacity"},
        'epotdist'     :{'value'  : 0.546 , 'comment': "!! Coefficient in exponential function for potential evapotranspiration's depth dependency"},
        'krs'          :{'value'  : 0.546 , 'comment': "!! parameter for estimating shortwave radiation used in the third petmodel"},
        'jhtadd'       :{'value'  : 0.546 , 'comment': "!! parameter for second petmodel"},
        'jhtscale'     :{'value'  : 0.546 , 'comment': "!! parameter for second petmodel"}
    }


    # evapoation land use

    Evap_land = {
        'section_head' :"""!!	-----																
!!																	
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19""",
        'kc3'          :{'value'  :[1.017511845,1.017511845,1.201224208,1.334493399,1.265059352,\
                                    1.020708799,1.020708799,1.020708799,1.020708799,1.265059352,\
                                    1.342448354,1.024959087,1.024959087,1.024959087,1.024959087,\
                                    1.024959087,1.024959087,1.024959087,1.024959087],
                         'comment': "!! crop coefficient for third petmodel"},
        'alb'          :{'value'  :[0.476534277,0.476534277,0.7,0.45542863,0.669192433,0.799822092,\
                                    0.799822092,0.799822092,0.799822092,0.669192433,0.400103867,\
                                    0.479658425,0.479658425,0.479658425,0.479658425,0.479658425,\
                                    0.479658425,0.479658425,0.479658425],
                         'comment': "!! albedo for petmodels"},
        'ttrig'        :{'value'  :[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                         'comment': "!! temperature threshold for soil temperature control on soil evapotranspiration"},
        'treda'        :{'value'  :[0.84,0.84,0.84,0.84,0.95,0.95,0.95,\
                                    0.95,0.95,0.7,0.9,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8],
                         'comment': "!! soil temperature control on soil evapotranspiration"},
        'tredb'        :{'value'  :[0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,\
                                    0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4],
                         'comment': "!! soil temperature control on soil evapotranspiration"},
        'cevp'         :{'value'  :[0.22,0.22,1.6,1.9,0.17,0.17,0.17,0.17,\
                                    0.17,0.1,0.21,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07],
                         'comment': "!! evapotranspiration parameter"},
        'fepotsnow'    :{'value'  :[0.912879467,0.912879467,0.18,0.533387661,\
                                    0.460848987,0.12002416,0.12002416,0.12002416,\
                                    0.12002416,0.460848987,0.206956849,0.197802201,\
                                    0.197802201,0.197802201,0.197802201,0.197802201,\
                                    0.197802201,0.197802201,0.197802201],
                         'comment': "!! fraction of snow-free potential evapotranspiration, used for calculation of snow evaporation"}
    }

    # Frozen soil infiltration for soil
    Forzen_soil_infil_soil = {
        'section_head' :"""!!======================================================																	
!! Frozen soil infiltration parameters																	
!! General and specific frozen soil 
!! SOIL:	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10	S11	S12					
!! for frozen soil submodel: 2""",
        'deepmem'      :{'value': 1000,  'comment': "!! deep soil temperature memory"},
        'bfroznsoil'   :{'value': [2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1],
                         'comment': "!! ??"},
        'logsatmp'     :{'value': [1.15,1.15,1.88,1.59,1.88,2.17,2.46,2.46,2.46,2.46,2.46,2.46],
                         'comment': "!! coefficient in unfrozen soil water content function"},
        'bcosby'       :{'value': [4.74,4.74,5.33,7.22,5.33,3.44,1.55,1.55,1.55,1.55,1.55,1.55],
                         'comment': "!! coefficient in unfrozen soil water content function"}
    }


    # Frozen soil infiltration for land use
    Forzen_soil_infil_LU = {
        'section_head' :"""!! ------------
!! Frozen soil infiltration parameters															
!!Land use parameters																	
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19""",
        'surfmem'      :{'value': [17.8,17.8,17.8,17.8,5.15,5.15,5.15,\
                                   5.15,5.15,5.15,5.15,5.15,5.15,5.15,\
                                   5.15,5.15,5.15,5.15,5.15],
                         'comment': "!! upper soil layer soil temperature memory"},
        'depthrel'     :{'value': [1.1152,1.1152,1.1152,1.1152,2.47,2.47,\
                                   2.47,2.47,2.47,2.47,2.47,2.47,2.47,2.47,\
                                   2.47,2.47,2.47,2.47,2.47],
                         'comment': "!! depth relation for soil temperature memory"},
        'frost'        :{'value': [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],
                         'comment': "!! frost depth parameter"}
    }



    # soil class parameters
    Soil_class_param = {
        'section_head' :"""!!	=======================================================================================================																
!!	"SOIL/LAND HYDRAULIC RESPONSE PARAMETERS - recession coef., water retention, infiltration, macropore, surface runoff; etc."																
!!	-----																
!!	Soil-class parameters																
!!	S1	S2	S3	S4	S5	S6	S7	S8	S9	S10	S11	S12					
!! surfacerunoff submodel: 0, soilleakage submodel: 0
!! recession coefficient for surface runoff should be set to one for lake and riverclasses with floodplains""",
        'rrcs1'        :{'value': [0.3201,0.2982,0.2663,0.451,0.1637,0.1637,0.1637,0.1637,0.1637,0.1637,0.1637,0.1637],
                         'comment': "!! recession coefficient for uppermost soil layer"},
        'rrcs2'        :{'value': [0.1612,0.0858,0.1422,0.0112,0.1914,0.1914,0.1914,0.1914,0.1914,0.1914,0.1914,0.1914],
                         'comment': "!! recession coefficient for lowest soil layer"},
        'trrcs'        :{'value': [0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15],
                         'comment': "!! recession coefficient for tile drains"},
        'mperc1'       :{'value': [63.9842,119.5863,93.9854,111.8318,20.1177,20.1177,20.1177,20.1177,20.1177,20.1177,20.1177,20.1177],
                         'comment': "!! maximum percolation capacity from soil layer one to soil layer two"},
        'mperc2'       :{'value': [97.6492,12.5429,20.0276,79.481,12.0754,12.0754,12.0754,12.0754,12.0754,12.0754,12.0754,12.0754],
                         'comment': "!! maximum percolation capacity from soil layer two to soil layer three"},
        'sfrost'       :{'value': [1,1,1,1,1,1,1,1,1,1,1,1],
                         'comment': "!! frost depth parameter"},
        'srrate'       :{'value': [0.4975,0.4489,0.1874,0.4499,0.4956,0.4956,0.4956,0.4956,0.4956,0.4956,0.4956,0.4956],
                         'comment': "!! fraction for surface runoff"},
        'wcwp1'        :{'value': [0.2732,0.214,0.1479,0.4233,0.4941,0.4941,0.4941,0.4941,0.4941,0.4941,0.4941,0.4941],
                         'comment': "!! wilting point as a fraction for uppermost soil layer"},
        'wcwp2'        :{'value': [0.2293,0.256,0.0984,0.4303,0.1308,0.1308,0.1308,0.1308,0.1308,0.1308,0.1308,0.1308],
                         'comment': "!! wilting point as a fraction for second soil layer"},
        'wcwp3'        :{'value': [0.3257,0.058,0.3639,0.0433,0.3632,0.3632,0.3632,0.3632,0.3632,0.3632,0.3632,0.3632],
                         'comment': "!! wilting point as a fraction for lowest soil layer"},
        'wcfc1'        :{'value': [0.4344,0.1758,0.3526,0.4812,0.4894,0.4894,0.4894,0.4894,0.4894,0.4894,0.4894,0.4894],
                         'comment': "!! fraction of soil available for evapotranspiration but not for runoff for uppermost soil layer"},
        'wcfc2'        :{'value': [0.1392,0.1966,0.3818,0.1163,0.2385,0.2385,0.2385,0.2385,0.2385,0.2385,0.2385,0.2385],
                         'comment': "!! fraction of soil available for evapotranspiration but not for runoff for second soil layer"},
        'wcfc3'        :{'value': [0.2307,0.2075,0.4055,0.077,0.343,0.343,0.343,0.343,0.343,0.343,0.343,0.343],
                         'comment': "!! fraction of soil available for evapotranspiration but not for runoff for lowest soil layer"},
        'wcep1'        :{'value': [0.8729,0.4168,0.8743,0.5142,0.0117,0.0117,0.0117,0.0117,0.0117,0.0117,0.0117,0.0117],
                         'comment': "!! effective porosity as a fraction for uppermost soil layer"},
        'wcep2'        :{'value': [0.1177,0.2773,0.0329,0.8547,0.087,0.087,0.087,0.087,0.087,0.087,0.087,0.087],
                         'comment': "!! effective porosity as a fraction for second soil layer"},
        'wcep3'        :{'value': [0.3064,0.8004,0.5832,0.474,0.8299,0.8299,0.8299,0.8299,0.8299,0.8299,0.8299,0.8299],
                         'comment': "!! effective porosity as a fraction for lowest soil layer"},
        'rrcs3'        :{'value': 0.1612,
                         'comment': "!! recession coefficient for slope dependence (upper layer) ???"}
    }


    # land use recession coeffcient

    Land_recession = {
        'section_head' :"""!! ------------
!! Frozen soil infiltration parameters															
!!Land use parameters																	
!!LUSE:	LU1	LU2	LU3	LU4	LU5	LU6	LU7	LU8	LU9	LU10	LU11	LU12	LU13	LU14	LU15	LU16	LU17	LU18	LU19""",
        'srrcs'        :{'value': [0.1259,0.0701,0.187,0.1977,0.0951,0.1208,0.1594,0.0694,0.1136,0.0575,1,0.1213,0.1213,0.1213,0.1213,0.1213,0.1213,0.1213,0.1213],
                         'comment': "!! (landuse) recession coefficient for surface runoff should be set to one for lake and riverclasses with floodplains"}
    }


    # Regional or deep Groundwater outflow

    Deep_groundwater = {
        'section_head' :"""!!	-----																
!!	Regional or deep groundwater outflow																
!! for deepground submodel: 0""",
        'rcgrw'        :{'value': 0.0,
                         'comment': "!! recession coefficient for regional groundwater outflow from soil layers (set to zero because we are not considering deepground)"}
    }


    # river routing

    River_routing ={
        'section_head' :"""!!																	
!!	=======================================================================================================																
!!	RIVER ROUTING																
!!	-----			
!! for riverflow submodel:1?""",
        'damp'         :{'value': 0.1614,
                         'comment': "!! fraction of delay in the watercourse which also causes damping"},
        'rivvel'       :{'value': 9.9267,
                         'comment': "!! celerity of flood in watercourse"},
        'qmean'        :{'value': 200,
                         'comment': "!! initial value for calculation of mean flow "}
    }



    # add to the par.txt
    def write_dictionary(output_file, dictionary, num = None):
        with open(output_file, 'a') as file:
            for key, value in dictionary.items():
                if key != 'section_head':
                    if isinstance(value['value'], list):
                        if not num is None: # check the length
                            if num != len(value['value']):
                                sys.exit('the list is not with the same length as soil\
    or land cover for key: '+ key)
                        values_line = '\t'.join(map(str, value['value']))
                        file.write(f"{key}\t{values_line}\t{value['comment']}\n")
                    else:
                        file.write(f"{key}\t{value['value']}\t{value['comment']}\n")
                else:
                    file.write(dictionary['section_head']+"\n")

    write_dictionary(output_file, header)
    write_dictionary(output_file, General_meteo_param)
    write_dictionary(output_file, Snow_param)
    write_dictionary(output_file, Snow_land_submodel_1, num = 19)
    write_dictionary(output_file, Snow_land_submodel_2, num = 19)
    write_dictionary(output_file, Glacier)
    write_dictionary(output_file, Glacier_melt)
    write_dictionary(output_file, Evap)
    write_dictionary(output_file, Evap_land, num = 19)
    write_dictionary(output_file, Forzen_soil_infil_soil, num = 12)
    write_dictionary(output_file, Forzen_soil_infil_LU, num = 19)
    write_dictionary(output_file, Soil_class_param, num = 12)
    write_dictionary(output_file, Land_recession, num = 19)
    write_dictionary(output_file, Deep_groundwater)
    write_dictionary(output_file, River_routing)


    file_path = output_file
    with open(file_path, 'r') as file:
        # Read the contents of the file
        file_contents = file.read()
        # Print the contents
        # print(file_contents)
        
################################################################
# write info and filedir files
def write_hype_info_filedir_files(path_to_save, spinup_days):
    # write filedir file
    output_file = path_to_save+'filedir.txt'

    if os.path.isfile(output_file):
        os.remove(output_file)

    with open(output_file, 'w') as file:
            file.write('./')
    # create results directory
    
    if not os.path.isdir(path_to_save+'/results'):
        os.makedirs(path_to_save+'/results')
        
    ###########

    # Output par to a .txt file

    output_file = path_to_save+'info.txt'
    if os.path.isfile(output_file):
        os.remove(output_file)


    # define start time, end time, based on input forcing
    # spinup period is a user defined inputs

    Pobs = pd.read_csv(path_to_save+'Pobs.txt', sep='\t', parse_dates=['time'])
    Pobs['time'] = Pobs['time'].dt.date
    start_date = Pobs['time'].iloc[0]
    end_date = Pobs['time'].iloc[-1]
    spinup_date = start_date + pd.Timedelta(days=spinup_days)


    # write out first text section
    s1= [
    """!! ----------------------------------------------------------------------------							
!!							
!! HYPE - Model Agnostic Framework
!!							
!! -----------------------------------------------------------------------------							
!! Check Indata during first runs (deactivate after first runs) """
    ]

    # write s1 in output file
    with open(output_file, 'w') as file:
        # Write the commented lines
        for line in s1:
            file.write(line + '\n')

    # create first dataframe
    df1_row=['indatacheckonoff','indatachecklevel']
    df1_val=[2,2]
    df1=pd.DataFrame(df1_val, index=df1_row, columns=None)

    # append df1
    with open(output_file, 'a') as file:
        # Write the DataFrame to the file
        df1.to_csv(file, sep='\t', index=True, header=False)#, line_terminator='\n')

    # write out first text section
    s2= [
    """!!
!! -----------------------------------------------------------------------------							
!!						
!! Simulation settings:							
!!							
!! -----------------	 """
    ]

    # write s2 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s2:
            file.write(line + '\n')

    # create df2
    df2_row=['bdate','cdate','edate','resultdir','instate', 'warning']
    df2_val=[start_date,spinup_date,end_date,'./results/', 'n','y']
    df2=pd.DataFrame(df2_val, index=df2_row, columns=None)

    # append df2
    with open(output_file, 'a') as file:
        # Write the DataFrame to the file
        df2.to_csv(file, sep='\t', index=True, header=False)#, line_terminator='\n')

    # write out s3
    s3= [
    """!! outstatedate """
    ]

    # write s3 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s3:
            file.write(line + '\n')

    # create df3
    df3_row=['readdaily','submodel','calibration','readobsid','soilstretch']
    df3_val=['y','n','n','n','n']
    df3=pd.DataFrame(df3_val, index=df3_row, columns=None)

    # append df3
    with open(output_file, 'a') as file:
        # Write the DataFrame to the file
        df3.to_csv(file, sep='\t', index=True, header=False)#, line_terminator='\n')

    # write out s4
    s4= [
    """!! Soilstretch enable the use of soilcorr parameters (strech soildepths in layer 2 and 3)
steplength	1d							
!! -----------------------------------------------------------------------------							
!!							
!! Enable/disable optional input files
!!							
!! -----------------					"""
    ]

    # write s4 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s4:
            file.write(line + '\n')

    # create df4
    df4_row=['readsfobs','readswobs','readuobs','readrhobs','readtminobs','readtmaxobs','soiliniwet','usestop84']
    df4_val=['n','n','n','n','y','y','n','n']
    df4=pd.DataFrame(df4_val, index=df4_row, columns=None)

    # create the corresponding comments
    c1 = [
        "!! For observed snowfall fractions in SFobs.txt",
        "!! For observed shortwave radiation in SWobs.txt",
        "!! For observed wind speeds in Uobs.txt",
        "!! For observed relative humidity in RHobs.txt",
        "!! For observed min air temperature in TMINobs.txt",
        "!! For observed max air temperature in TMAXobs.txt",
        "!! initiates soil water to porosity instead of field capacity which is default (N). Set Y to use porosity.",
        "!! initiates soil water to porosity instead of field capacity which is default (N). Set Y to use porosity.",
    ]

    # append c1 and df4
    with open(output_file, 'a') as file:
        # Iterate over DataFrame rows
        for i, (index, row) in enumerate(df4.iterrows()):
            # Check if there is a comment line for the current row
            if i < len(c1):
                # Write the row name, values, and comment on the same line
                line = str(index) + '\t' + '\t'.join(str(val) for val in row.values) + '\t' + c1[i] + '\n'
            else:
                # Write the row name and values without comment on the same line
                line = str(index) + '\t' + '\t'.join(str(val) for val in row.values) + '\n'

            # Write the line to the file
            file.write(line)

    # write out s5
    s5= [
    """!! -----------------------------------------------------------------------------							
!!							
!! Define model options (optional)
!!							
!! -----------------							
!!snowfallmodel:								
!!                  0 threshold temperature model							
!!                  1 inputdata (SFobs.txt)							
!!snowmeltmodel:							
!!                  0,1 temperature index             (with/without snowcover scaling)							
!!                  2   temperature + radiation index (with/without snowcover scaling)							
!!							
!!  snowevapmodel   0 off							
!!                  1 on							
!!                   							
!!  petmodel:  (potential evapotranspiration) (is shown in geodata for WWH)							
!!                  0 original HYPE temperature model (with Xobs epot replacement)							
!!                  1 original HYPE temperature model (without Xobs epot replacement)							
!!                  2 Modified Jensen-Haise 							
!!                  3 Modified Hargreaves-Samani							
!!                  4 Priestly-Taylor							
!!                  5 FAo Penman-Monteith							
!!							
!! lakeriverice:							
!!                  0 off							
!!                  1 on, old (simple) air-water heat exchange              (requires T2 water temperature model)							
!!                  2 on, new heatbalance model for air-water heat exchange (requires T2 water temperature model)							
!!							
!! substance T2     switching on the new water temperature trace model							
!!							
!! deepground       0   off    Deep groundwater (Aquifer) model options							
!!                  1,2 on
!! Glacierini	0 off 1 on	(1 used for statefile preparation)	
!! Floodplain		0, 1, 2, 3 (3 used for WWH)					
!! -----------------							"""
    ]

    # write s5 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s5:
            file.write(line + '\n')

    # create df5
    df5_row=['modeloption snowfallmodel','modeloption snowdensity','modeloption snowfalldist',
             'modeloption snowheat','modeloption snowmeltmodel','modeloption snowevaporation','modeloption lakeriverice',
             'modeloption deepground','modeloption glacierini','modeloption frozensoil',
             'modeloption infiltration','modeloption surfacerunoff','modeloption petmodel',
             'modeloption riverflowmodel','modeloption soilleakage']
    df5_val=[0,0,0,1,2,1,0,0,1,2,3,0,1,0,0]
    df5=pd.DataFrame(df5_val, index=df5_row, columns=None)

    # append df5
    with open(output_file, 'a') as file:
        # Write the DataFrame to the file
        df5.to_csv(file, sep='\t', index=True, header=False)#, line_terminator='\n')

    # write out s6
    s6= [
    """!! ------------------------------------------------------------------------------------							
!!							
!! Define outputs
!!							
!! -----------------							
!! meanperiod 1=daymean, 2=weekmean, 3=monthmean, 4=yearmean, 5=total period mean							
!! output variables: see http://www.smhi.net/hype/wiki/doku.php?id=start:hype_file_reference:info.txt:variables 
!! -----------------							
!! BASIN outputs 
!! The present basins are some large rivers distributed over different continents
!! -----------------							"""
    ]

    # write s6 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s6:
            file.write(line + '\n')

    df6={'!! basinoutput variable': 'cout\trout\tctmp\tsnow\tsdep\tsoim\tsom2\tsml1\tsml2\tsml3\tsmrz\tsm13\tstsw\tsrff\tsmfd\tsrfd\tsmfp\tsrfp\tsmdf\tgwat\tsfst\tstmp\tstm1\tstm2\tstm3\tcfsc\tsmax\tcilv\tclbv\tcgwl\tcloc\tclof\tclrv\tcmrv\tqerr\tcobc\tcmri\tclri\tcmrb\tclrb\tcmrs\tclrs\tclic\tglcv\tglca\tlrdp\tmrdp\tcgmb\tcgma\tC106\tC108\tC111\tC114\tC206\tC208\tC211\tC214\tcoT1\tcoT2\tcprc\tcpSF\tcpRF\tevap\tepot\ticpe\tevsn\tlevp\tevpt\tpsim\tcrun\tcro1\tcro2\tcro3\tcrod\tcros\tros1\tros2\tacdf\taqin\taqut\tspeq\tgmlt\tloff\tlrfa\tmrfa\tlred\tmred\tinfi\tsnwc\tsnht\tsnte\tsnts\tdtmp\tcmrp\tcrgl\tcrnt\tcmrr\tcrpt\tcrex\tpsub\tesub\tisub\tcmrq',
         '!! basinoutput meanperiod':1,
         '!! basinoutput decimals':3,
         '!! basinoutput subbasin':'subid1\tsubid2\tsubid3\tsubid4',
         '!! printwaterbal':'N' }

    df6

    # append df6
    with open(output_file, 'a') as file:
        for key,value in df6.items():
            a=str(key)+'\t'+str(value)+'\n'
            file.write(a)


    # write out s7
    s7= [
    """!! -----------------							
!! TIME outputs 
!! -----------------	"""
    ]

    # write s7 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s7:
            file.write(line + '\n')

    df7={'timeoutput variable': 'cout\tevap\tsnow', 
         'timeoutput meanperiod':1,
         'timeoutput decimals':3}

    df7

    with open(output_file, 'a') as file:
        for key,value in df7.items():
            a=str(key)+'\t'+str(value)+'\n'
            file.write(a)


    # write out s8
    s8= [
    """!! -----------------							
!! MAP outputs
!! -----------------							
!! mapoutput variable	cout cprc ctmp
!! mapoutput decimals	3						
!! mapoutput meanperiod	5						
!! ------------------------------------------------------------------------------------							
!!							
!! Select criteria for model evaluation and automatic calibration
!!							
!! -----------------							
!! General settings
!! -----------------			"""
    ]

    # write s8 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s8:
            file.write(line + '\n')

    df8={'!! crit meanperiod': 1, 
         '!! crit datalimit':30,
         '!! crit subbasin':'subid1\tsubid2\tsubid3\tsubid4'}

    df8

    with open(output_file, 'a') as file:
        for key,value in df8.items():
            a=str(key)+'\t'+str(value)+'\n'
            file.write(a)


    # write out s9
    s9= [
    """!! -----------------			
!! Criterion-specific settings
!! -----------------				"""
    ]

    # write s9 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s9:
            file.write(line + '\n')

    # create df9 for basin outputs
    d1=['MKG']
    d2=['cout']
    d3=['rout']
    d4=[1]
    df9_row=['crit 1 criterion','crit 1 cvariable','crit 1 rvariable','crit 1 weight']
    df9=pd.DataFrame([d1,d2,d3,d4], index=df9_row)

    # append df9
    with open(output_file, 'a') as file:
        # Write the DataFrame to the file
        df9.to_csv(file, sep='\t', index=True, header=False)#, line_terminator='\n')