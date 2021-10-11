#!/usr/bin/env python3

import os
import sys
sys.path.append('/foss_fim/src')
from multiprocessing import Pool
import argparse
from utils.reproject_dem import reproject_dem
from utils.shared_functions import update_raster_profile
import gdal
import rasterio
import numpy as np
from rasterio.warp import calculate_default_transform, reproject, Resampling
from pyproj.crs import CRS

#from utils.shared_variables import PREP_PROJECTION, PREP_PROJECTION_CM
PREP_PROJECTION_CM = 'PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-96.0],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_origin",23.0],UNIT["Meter",1.0],VERTCS["NAVD_1988",VDATUM["North_American_Vertical_Datum_1988"],PARAMETER["Vertical_Shift",0.0],PARAMETER["Direction",1.0],UNIT["Centimeter",0.01]]]'
PREP_PROJECTION = 'PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Reproject Elevation rasters and update profile')
    parser.add_argument('-dem_dir','--dem-dir', help='DEM filename', required=True,type=str)
    parser.add_argument('-j','--number-of-jobs',help='Number of processes to use. Default is 1.',required=False, default="1",type=int)
    parser.add_argument('-nodata','--nodata-val', help='DEM nodata value', required=False,type=float,default=-9999.0)
    parser.add_argument('-block','--blocksize', help='DEM blocksize', required=False,type=int,default=512)
    parser.add_argument('-keep','--keep-intermediate', help='keep intermediate files', required=False,type=bool,default=True)

    args = vars(parser.parse_args())

    dem_dir            = args['dem_dir']
    number_of_jobs     = args['number_of_jobs']
    nodata_val         = args['nodata_val']
    blocksize          = args['blocksize']
    keep_intermediate  = args['keep_intermediate']

    reproject_procs_list = []

    for huc in os.listdir(dem_dir):
        raster_dir = os.path.join(dem_dir,huc)
        elev_cm = os.path.join(raster_dir, 'elev_cm.tif')
        elev_cm_proj = os.path.join(raster_dir, 'elev_cm_proj.tif')
        
        print('reprojecting raster')
        input_raster = gdal.Open(elev_cm)
        #output_raster = os.path.join(output_name,config_mode,"tmp",os_path_huc+raster_name+"_proj.tif")
        warp = gdal.Warp(elev_cm_proj,input_raster,dstSRS=PREP_PROJECTION_CM,resampleAlg="cubic")
        warp = None # Closes the files
        
        print('changing elevation')
        dem_cm = rasterio.open(elev_cm_proj)
        print('grab nodata')
        no_data = dem_cm.nodata
        print('read data')
        data = dem_cm.read(1)
        print('cm to m calc')
        dem_m = np.where(data == int(no_data), -9999.0, (data/100).astype(rasterio.float32))
        print('profile')
        dem_m_profile = dem_cm.profile.copy()
        print('write')
        with rasterio.open(os.path.join(raster_dir, 'elev_m.tif'), "w", **dem_m_profile, BIGTIFF='YES') as dest:
            dest.write(dem_m, indexes = 1)
            
    print('fin')
