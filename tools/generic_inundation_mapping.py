import os
import sys
import glob
import argparse
import shutil
#from inundation import inundate
import ray
import math
import numpy as np
import pandas as pd
from numba import njit, typed, types
#from concurrent.futures import ThreadPoolExecutor,as_completed
#from subprocess import run
import netCDF4
from os.path import splitext
import rasterio
from rasterio.merge import merge
from rasterio.mask import mask
from rasterio.io import DatasetReader,DatasetWriter
from rasterio.features import shapes
import rioxarray as rxr
import fiona
from shapely.geometry import shape
from collections import OrderedDict
import argparse
#from warnings import warn
from osgeo import gdal
#from gdal import BuildVRT
import geopandas as gpd
import urllib.request

import xarray as xr
import fsspec
import numpy as np
import pandas as pd
import fiona
import s3fs
import datetime
from dask.distributed import Client

ray.init()

'''
aoi_bb (ex [-80.6,-79.84,32.42,33.32] in epsg 4356) -or- list of hucs -or- [rem, catch, and hydro_tables] 
    -and-
forcast file (csv, netcdf, or )
    -and-
output_name
config: One of 'fr', 'ms', 'fr_ms'
clip_shape
out_epsg
'''
def map_inundation(
        aoi_bb=None,
        huc_list=None,
        rem_file=None,
        catchment_file=None,
        catchment_poly=None,
        hydro_table=None,
        config="fr",
        fimrun_fr_data_folder = "/home/user/Desktop/cahaba/data_local/",
        fimrun_ms_data_folder = "/home/user/Desktop/cahaba/data_local/",
        output_name = "/home/user/Desktop/cahaba/data_out/test/",
        forecast_source="NWM2.1_retro",
        forecast_start_time="202010011200",
        forecast_end_time="202010011200",
        mask_type=None,
        hucs=None,
        hucs_layerName=None,
        subset_hucs=None,
        num_workers=1,
        aggregate=False,
        inundation_raster=None,
        inundation_polygon=None,
        depths=None,
        out_raster_profile=None,
        out_vector_profile=None,
        quiet=False
            ):
    
    def normalize_huclist(inlist,length=8):
        # Load normalization data (a list of all huc12's)
        huc_data_path = '/home/user/Desktop/cahaba/tools/aux_data/huc_features.gpkg'
        huc12 = gpd.read_file(huc_data_path, layer = 'huc12')
        huc12_values = huc12.values

        # Find short hucs in give list
        short_inlist = [huc for huc in inlist if len(huc) < length]

        # remove those small ones
        inlist=[x for x in inlist if len(x) != length]

        # Find and append all valid hucs within that list
        for i in short_inlist:
            inlist.extend(list(filter(lambda x: x.startswith(i), huc12.HUC12.values)))

        # Truncate long hucs
        inlist = [item[:length] for item in inlist]

        return set(inlist)

    @ray.remote
    def project_huc(in_raster_path,epsg="3857",resample_type="nearest",config_mode='fr'):
        # grab huc id
        os_path_str = str(in_raster_path)
        os_path_str = os_path_str.split(os.path.sep)
        os_path_huc = os_path_str[-2]

        # grab raster type
        raster_name = os_path_str[-1][:-4]

        # Warp raster to final crs
        input_raster = gdal.Open(in_raster_path)
        output_raster = os.path.join(output_name,config_mode,"tmp",os_path_huc+raster_name+"_proj.tif")
        warp = gdal.Warp(output_raster,input_raster,dstSRS='EPSG:'+epsg,resampleAlg=resample_type)
        warp = None # Closes the files
        return

    def merge_src(src_files,config_mode):
        with open(os.path.join(output_name,config_mode,"base_data","master_src.csv"), 'wb') as outfile:
            for i, fname in enumerate(src_files):
                with open(fname, 'rb') as infile:
                    if i != 0:
                        infile.readline()  # Throw away header on all but first file
                    # Block copy rest of file from input to output without parsing
                    shutil.copyfileobj(infile, outfile)
        return

    def gather_retro_20_data(desired_comid,time_start,time_end):
        #Based on: https://registry.opendata.aws/nwm-archive/
        import re
        import xarray as xr
        import fsspec
        import numpy as np
        import pandas as pd
        import fiona
        import s3fs
        import datetime
        from dask.distributed import Client
        client = Client()

        url = 's3://noaa-nwm-retro-v2-zarr-pds'
        ds = xr.open_zarr(fsspec.get_mapper(url, anon=True),consolidated=True)
        streamflow_array = ds['streamflow'].sel(time='2017-08-29 00:00')
        streamflow_pd = streamflow_array.to_pandas()
        streamflow_pd = pd.DataFrame({'feature_id':streamflow_pd.index, 'discharge':streamflow_pd.values})
        #streamflow_pd = streamflow_pd[streamflow_pd.feature_id.isin(catch_indices)]
        streamflow_pd.to_csv("/home/user/Desktop/2017_08_29_00.csv", index=False)
        client.shutdown()
        return

    def gather_retro_21_data(desired_comid,time_start):
        #https://noaa-nwm-retrospective-2-1-pds.s3.amazonaws.com/model_output/1982/198212312000.CHRTOUT_DOMAIN1.comp
        urllib.request.urlretrieve(
            "https://noaa-nwm-retrospective-2-1-pds.s3.amazonaws.com/model_output/"+time_start[0: 4]+"/"+time_start+".CHRTOUT_DOMAIN1.comp",
            os.path.join(output_name,config,"CHRTOUT_DOMAIN1.comp"))
        flows_nc = xr.open_dataset(os.path.join(output_name,config,"CHRTOUT_DOMAIN1.comp"))
        flows_df = flows_nc.to_dataframe()
        flows_df.reset_index(inplace=True)
        flows_df = flows_df.drop(columns=["time","reference_time","crs","latitude","longitude","order","elevation","q_lateral","velocity","qSfcLatRunoff","qBucket","qBtmVertRunoff"])
        flows_df = flows_df[flows_df.feature_id.isin(desired_comid)]
        flows_df = flows_df.rename(columns={"streamflow": "discharge"})
        #flows_df = flows_df.astype({'streamflow':"str","discharge":"float"}).dtypes
        convert_dict = {'feature_id': str,'discharge': float}
        flows_df = flows_df.astype(convert_dict)
        return flows_df

    def gather_nwm_forecast(desired_comid,forecast_source,time_steps):
        return
    
    def process_base_data_for_config():
        if(config in ["fr","fr_ms"]):
            os.makedirs(os.path.join(output_name,"fr","base_data"), exist_ok=True)
            os.makedirs(os.path.join(output_name,"fr","tmp"), exist_ok=True)
            
            print('Reprojecting fr rem data')
            rem_file_list = [os.path.join(fimrun_fr_data_folder, huc, "rem_zeroed_masked.tif") for huc in hucs_to_process]
            ray.get([project_huc.remote(huc,resample_type="cubic",config_mode='fr') for huc in rem_file_list])
            print('Merging fr rem data')
            rem_proj_file_list = glob.glob(os.path.join(output_name,'fr',"tmp","*"+"rem"+"*"+"_proj.tif"))
            rem_vrt = gdal.BuildVRT(os.path.join(output_name,'fr',"tmp","rem.vrt"),rem_proj_file_list)
            gdal.Translate(os.path.join(output_name,'fr',"base_data","rem.tif"),rem_vrt)
            rem_vrt = None
            
            print('Reprojecting fr catch data')
            catch_file_list = [os.path.join(fimrun_fr_data_folder, huc, "gw_catchments_reaches_filtered_addedAttributes.tif") for huc in hucs_to_process]
            ray.get([project_huc.remote(huc,config_mode='fr') for huc in catch_file_list])
            print('Merging fr catch data')
            catch_proj_file_list = glob.glob(os.path.join(output_name,'fr',"tmp","*"+"catchments"+"*"+"_proj.tif"))
            catch_vrt = gdal.BuildVRT(os.path.join(output_name,'fr',"tmp","catchments.vrt"),catch_proj_file_list)
            gdal.Translate(os.path.join(output_name,'fr',"base_data","catchments.tif"),catch_vrt)
            catch_vrt = None
            
            # merge SRC
            print('Merging fr rating curve data')
            src_file_list = [os.path.join(fimrun_fr_data_folder,huc,"hydroTable.csv") for huc in hucs_to_process]
            merge_src(src_file_list,config_mode='fr')
            print("AOI base data for fr processed")
            print("")

        if(config in ["ms","fr_ms"]):
            os.makedirs(os.path.join(output_name,"ms","base_data"), exist_ok=True)
            os.makedirs(os.path.join(output_name,"ms","tmp"), exist_ok=True)
            
            print('Reprojecting ms rem data')
            rem_file_list = [os.path.join(fimrun_ms_data_folder, huc, "rem_zeroed_masked.tif") for huc in hucs_to_process]
            ray.get([project_huc.remote(huc,resample_type="cubic",config_mode='ms') for huc in rem_file_list])
            print('Merging ms rem data')
            rem_proj_file_list = glob.glob(os.path.join(output_name,"ms","tmp","*"+"rem"+"*"+"_proj.tif"))
            rem_vrt = gdal.BuildVRT(os.path.join(output_name,"ms","tmp","rem.vrt"),rem_proj_file_list)
            gdal.Translate(os.path.join(output_name,"ms","base_data","rem.tif"),rem_vrt)
            rem_vrt = None

            print('Reprojecting ms catch data')
            catch_file_list = [os.path.join(fimrun_ms_data_folder, huc, "gw_catchments_reaches_filtered_addedAttributes.tif") for huc in hucs_to_process]
            ray.get([project_huc.remote(huc,config_mode='ms') for huc in catch_file_list])
            print('Merging ms catch data')
            catch_proj_file_list = glob.glob(os.path.join(output_name,"ms","tmp","*"+"catchments"+"*"+"_proj.tif"))
            catch_vrt = gdal.BuildVRT(os.path.join(output_name,"ms","tmp","catchments.vrt"),catch_proj_file_list)
            gdal.Translate(os.path.join(output_name,"ms","base_data","catchments.tif"),catch_vrt)
            catch_vrt = None

            # merge SRC
            print('Merging ms rating curve data')
            src_file_list = [os.path.join(fimrun_ms_data_folder, huc, "hydroTable.csv") for huc in hucs_to_process]
            merge_src(src_file_list,config_mode='ms')
            print("AOI base data for ms processed")
            print("")

        if(config == "fr_ms"):
            print("Merging fr_ms rem")
            rem_datasets_in_order = [os.path.join(output_name,"ms","base_data","rem.tif"),
                                     os.path.join(output_name,"fr","base_data","rem.tif")]
            src_files_to_mosaic = []
            for rem_dataset in rem_datasets_in_order:
                src = rasterio.open(rem_dataset)
                src_files_to_mosaic.append(src)
            mosaic, out_trans = merge(src_files_to_mosaic)
            out_meta = src.meta.copy()
            out_meta.update({"driver": "GTiff",
                             "height": mosaic.shape[1],
                             "width": mosaic.shape[2],
                             "transform": out_trans
                            }
                           )
            with rasterio.open(os.path.join(output_name,"fr_ms","base_data","rem.tif"), "w", **out_meta) as dest:
                dest.write(mosaic)
                
            print("Merging fr_ms catchments")
            catch_datasets_in_order = [os.path.join(output_name,"ms","base_data","catchments.tif"),
                                       os.path.join(output_name,"fr","base_data","catchments.tif")]
            src_files_to_mosaic = []
            for catch_dataset in catch_datasets_in_order:
                src = rasterio.open(catch_dataset)
                src_files_to_mosaic.append(src)
            mosaic, out_trans = merge(src_files_to_mosaic)
            out_meta = src.meta.copy()
            out_meta.update({"driver": "GTiff",
                             "height": mosaic.shape[1],
                             "width": mosaic.shape[2],
                             "transform": out_trans
                            }
                           )
            with rasterio.open(os.path.join(output_name,"fr_ms","base_data","catchments.tif"), "w", **out_meta) as dest:
                dest.write(mosaic)
                
            print("Merging fr_ms rating curve data")
            src_tables_in_order = [os.path.join(output_name,"ms","base_data","master_src.csv"),
                                   os.path.join(output_name,"fr","base_data","master_src.csv")]
            merge_src(src_tables_in_order,config_mode='fr_ms')

        return
    
    def discharge_to_stage_table(single_value_flows):
        hydro_table_file = full_file_list[2]
        hydro_table = pd.read_csv(hydro_table_file, dtype={'HydroID':str,'feature_id':str,'stage':float,'discharge_cms':float})
        forecast_table = single_value_flows

        #hydro_table = hydro_table[hydro_table["LakeID"] == -999]  # Subset hydroTable to include only non-lake catchments.
        hydro_table = pd.merge(hydro_table, forecast_table, on='feature_id', how='left')
        hydro_table.set_index(['HydroID'],inplace=True)

        # initialize dictionary
        catchment_stages_dict = typed.Dict.empty(types.int32,types.float64)

        # interpolate stages
        #for hid,sub_table in hydro_table.groupby(by='feature_id'):
        for hid,sub_table in hydro_table.groupby(level='HydroID'):    
            interpolated_stage = np.interp(sub_table.loc[:,'discharge'].unique(),
                                           sub_table.loc[:,'discharge_cms'],sub_table.loc[:,'stage'])
            # add this interpolated stage to catchment stages dict
            h = round(interpolated_stage[0],4)

            hid = types.int32(hid) ; h = types.float64(h)
            catchment_stages_dict[hid] = h

        catchment_stages_dict.update({0:np.nan}) 
        return catchment_stages_dict
    
    def stage_to_rWSE_raster(catch_reclass_dict):
        in_file = full_file_list[1]

        with rasterio.open(in_file) as src:
            # Read as numpy array
            array = src.read()
            profile = src.profile
            profile["dtype"] = "float64"

            # Reclassify in a single operation using broadcasting
            reclass_array = np.vectorize(catch_reclass_dict.__getitem__)(array)

        with rasterio.open(os.path.join(output_name, config,"outputs","rWSE.tif"), 'w', **profile) as dst:
            dst.write(reclass_array)
        return 

    def map_innundation_depths():
        rwse = rxr.open_rasterio(os.path.join(output_name, config,"outputs","rWSE.tif"), masked=True).squeeze()
        rem = rxr.open_rasterio(os.path.join(output_name, config,"base_data","rem.tif"), masked=True).squeeze()

        water_depth = rwse - rem
        water_depth.rio.to_raster(os.path.join(output_name, config,"outputs","Flooding.tif"))
        return 
    
    # Test if folder paths valid
    os.makedirs(os.path.join(output_name,config,"base_data"), exist_ok=True)
    os.makedirs(os.path.join(output_name,config,"tmp"), exist_ok=True)
    needed_mapping_files = os.path.join(output_name,config,"base_data")
    
    # If we need to generate the preprocessed data for this AOI
    if(True):
        # Figure out what hucs need to be run
        if((aoi_bb is not None) & (huc_list is not None)):
            print('Too many AOI parameters provided, defaulting to huc list')
            # Maybe print/generate a map?
        if(aoi_bb is not None):
            print("Subsetting bounding box to relevent huc 8 units")
            huc_data_path = '/home/user/Desktop/cahaba/tools/aux_data/huc_features.gpkg'
            huc8 = gpd.read_file(huc_data_path, layer = 'huc8')
            hucs_to_process = huc8.cx[aoi_bb[0]:aoi_bb[1],aoi_bb[2]:aoi_bb[3]].HUC8.values
        elif(huc_list is not None):
            # edge case: huc value not really a huc?
            if(all(len(x) == 8 for x in huc_list)):
                hucs_to_process = huc_list
            else:
                print("Normalizing provided hucs to relevent huc 8 units")
                hucs_to_process = normalize_huclist(huc_list,8)
        
        # If we need to generate the base data for this AOI
        # run_fim_by_unit()
        
    print("Mapping "+str(len(hucs_to_process))+" huc 8's")
        
    file_names_list = ['rem.tif', 'catchments.tif', 'master_src.csv']
    full_file_list = [os.path.join(output_name, config,"base_data",file) for file in file_names_list]
    if(all(list(map(os.path.isfile,full_file_list)))==False):
        process_base_data_for_config()
        
    print("Base mapping data collected - Grabing forecasts")
    os.makedirs(os.path.join(output_name,config,"outputs"), exist_ok=True)
    
    # Load forecasts
    if(forecast_source=="NWM2.1_retro"):
        hydrotable_feature_id = pd.read_csv(full_file_list[2])
        desired_comids = hydrotable_feature_id['feature_id'].unique()
        flows = gather_retro_21_data(desired_comids,forecast_start_time)
        
    print("Transforming discharge to stage")
    catch_reclass_dict = discharge_to_stage_table(flows)
    print("Reclassifing catchment raster")
    stage_to_rWSE_raster(catch_reclass_dict)
    print("Mapping flood depths")
    map_innundation_depths()
    
    # Data cleanup
    
    return print("Inundation written to "+os.path.join(output_name,config,"outputs"))

#map_inundation(aoi_bb=[-80.6,-79.84,32.42,33.32],dataoutput_folder = "harvey")
map_inundation(#huc_list=["12040201","12040203"],
               aoi_bb=[-80.6,-79.84,32.42,33.32],
               #config="fr",
               config="fr_ms",
               fimrun_fr_data_folder = "/home/user/Desktop/cahaba/data_local/",
               fimrun_ms_data_folder = "/home/user/Desktop/cahaba/data_local/",
               output_name = "/home/user/Desktop/cahaba/data_out/test/",
               forecast_source="NWM2.1_retro",
               forecast_start_time="199210240900")


#src_file_list = [os.path.join(datafolder, huc, "hydroTable.csv") for huc in hucs]
