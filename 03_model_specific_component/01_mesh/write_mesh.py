# %%
import meshflow
import os

# %%
print('Running MESH flow')

# %%
# main work path
work_path = '/scratch/mia725/SCRB'

config = {
    'riv': os.path.join(work_path, 'MERIT_geofabric', 'extracted_rivers.shp'),
    'cat': os.path.join(work_path, 'MERIT_geofabric', 'extracted_subbasins.shp'),
    'landcover': os.path.join(work_path, 'gistool-outputs', 'modified_domain_stats_NA_NALCMS_landcover_2020_30m.csv'),
    'forcing_files': os.path.join(work_path, 'easymore-outputs'),
    'forcing_vars': [ # does the variable list, match those of the "agnostic" step?
        "RDRS_v2.1_P_P0_SFC",
        "RDRS_v2.1_P_HU_09944",
        "RDRS_v2.1_P_TT_09944",
        "RDRS_v2.1_P_UVC_09944",
        "RDRS_v2.1_A_PR0_SFC",
        "RDRS_v2.1_P_FB_SFC",
        "RDRS_v2.1_P_FI_SFC",
    ],
    'forcing_units': { # Here, enter RDRS's original variable units
        'RDRS_v2.1_P_P0_SFC': 'millibar',
        'RDRS_v2.1_P_HU_09944': 'kg/kg',
        'RDRS_v2.1_P_TT_09944': 'celsius',
        'RDRS_v2.1_P_UVC_09944': 'knot',
        'RDRS_v2.1_A_PR0_SFC': 'm/hr',
        'RDRS_v2.1_P_FB_SFC': 'W/m^2',
        'RDRS_v2.1_P_FI_SFC': 'W/m^2',
    },
    'forcing_to_units': { # And here, the units that MESH needs to read
         'RDRS_v2.1_P_UVC_09944': 'm/s',
         'RDRS_v2.1_P_FI_SFC': 'W/m^2',
         'RDRS_v2.1_P_FB_SFC': 'W/m^2',
         'RDRS_v2.1_A_PR0_SFC': 'mm/s',
         'RDRS_v2.1_P_P0_SFC': 'pascal',
         'RDRS_v2.1_P_TT_09944': 'kelvin',
         'RDRS_v2.1_P_HU_09944': 'kg/kg',
    },
    'main_id': 'COMID', # what is the main ID of each river segment? Column name in the `cat` Shapefile
    'ds_main_id': 'NextDownID', # what is the downstream segment ID for each river segment? ditto.
    'landcover_classes': { # these are the classes defined for NALCMS-Landsat 2015 dataset. Is this accurate?
        0: 'Unknown',
        1: 'Temperate or sub-polar needleleaf forest',
        2: 'Sub-polar taiga needleleaf forest',
        3: 'Tropical or sub-tropical broadleaf evergreen forest',
        4: 'Tropical or sub-tropical broadleaf deciduous forest',
        5: 'Temperate or sub-polar broadleaf deciduous forest',
        6: 'Mixed forest',
        7: 'Tropical or sub-tropical shrubland',
        8: 'Temperate or sub-polar shrubland',
        9: 'Tropical or sub-tropical grassland',
        10: 'Temperate or sub-polar grassland',
        11: 'Sub-polar or polar shrubland-lichen-moss',
        12: 'Sub-polar or polar grassland-lichen-moss',
        13: 'Sub-polar or polar barren-lichen-moss',
        14: 'Wetland',
        15: 'Cropland',
        16: 'Barren lands',
        17: 'Urban',
        18: 'Water',
        19: 'Snow and Ice',
    },
    'ddb_vars': { # the stuff that MESH needs: slope, river length, etc... Let me know if there is any issues here!
        'slope': 'ChnlSlope',
        'lengthkm': 'ChnlLength',
        'Rank': 'Rank',
        'Next': 'Next',
        'landcover': 'GRU',
        'unitarea': 'GridArea',
        'landcover_names': 'LandUse',
    },
    'ddb_units': {
        'ChnlSlope': 'm/m',
        'ChnlLength': 'km', # is it in km or m? Please check the units of the Shapefile you created!
        'Rank': 'dimensionless',
        'Next': 'dimensionless',
        'GRU': 'dimensionless',
        'GridArea': 'km^2', # what was the unit of the GridArea, or Shape_Area in the `catchments` Shapefile?
        'LandUse': 'dimensionless',
    },
    'ddb_to_units': {
        'ChnlSlope': 'm/m',
        'ChnlLength': 'm', # This is what MESH needs, no need to change.
        'Rank': 'dimensionless',
        'Next': 'dimensionless',
        'GRU': 'dimensionless',
        'GridArea': 'm^2', # This is what MESH needs, no need to change.
        'LandUse': 'dimensionless',
    },
    'ddb_min_values': {
        'ChnlSlope': 1e-10, # in case there are 0s in the `rivers` Shapefile, we need minimums for certain variables
        'ChnlLength': 1e-3,
        'GridArea': 1e-3,
    },
    'gru_dim': 'NGRU', # change to `NGRU` for 'MESH>=r1860', keep for 'MESH<=1860', for example for r1813.
    'hru_dim': 'subbasin',
    'outlet_value': -9999,
}

# %%
exp1 = meshflow.MESHWorkflow(**config)

# %%
exp1.run()

# %%
print('Finished processing MESH inputs')

# %%
# exp1.forcing

# %%
# exp1.ddb

# %%
if not os.path.isdir(work_path+'/MESH/'):
    os.makedirs(work_path+'/MESH/')
exp1.save(work_path+'/MESH/')

# %%
# copy MESH static files (using the terminal)
# !cp -r setting_files/* {work_path}/MESH/
os.system("cp -r setting_files/* "+ work_path +'/MESH/')
# !mkdir -p {work_path}/MESH/results/
os.makedirs(work_path+'/MESH/results/')

# %% [markdown]
# ____


