## Flood extent stack rasters and vectors
# This notebook implements two inundation products: 
# flood extent polygons and flood depth grids. Given a set of stage heights, 
# a flood extent polygon output denotes the extent of inundated areas as 
# polygons corresponding to a stage height. A flood depth grid, on the other 
# hand, is the raster whose cell value is $stage\_height - hand$. Only cells 
# with non-negative flood depth values (inundated) are output as data cells.  

import os, sys, time
import numpy as np
import pandas as pd
import shapely
import shapely.geometry
import shapely.errors
import fiona
import rasterio
import rasterio.mask
import rasterio.features
import affine

from matplotlib import pyplot
from matplotlib import cm

from collections import OrderedDict # for output schema

import rtree
import math
import csv

# create flood extent (binary) polygons
def floodstack(f_hand=None, f_hfmap=None, f_catch=None, stage_start=1.0, num_stages=80, stage_len=1.0, odir=None, prefix=None):
    # load hand
    t0 = time.perf_counter()
    # raster_file = ddir + '/' + 'rem_zeroed_masked.tif'
    src = rasterio.open(f_hand)
    hand = src.read(1)
    
    # read hydroid-feature_id mapping
    # hfmap_file = ddir + '/' + 'crosswalk_table.csv'
    df_hfmap = pd.read_csv(f_hfmap)
    hfmap = {}
    for index, row in df_hfmap.iterrows():
        hydroid = row['HydroID']
        fid = row['feature_id']
        if hydroid in hfmap:
            if fid != hfmap[hydroid]:
                print("Error, hydroid ", hydroid, " has two feature_ids", fid, hfmap[hydroid])
        else:
            hfmap[hydroid] = fid
    
    # load catchment raster
    # catch_file = ddir + '/' + 'gw_catchments_reaches_filtered_addedAttributes.tif'
    r_catch = rasterio.open(f_catch)
    catch = r_catch.read(1)
    t1 = time.perf_counter()
    print('TIME', 'input_reading', (t1-t0), 'seconds', flush=True)

    #stage = 1.0 # 3 feet
    m2ft = 3.280839895
    oneft = 1.0 / m2ft
    results = []
    t0 = time.perf_counter()
    for i in range(num_stages):
        t00 = time.perf_counter()
        
        stage = (stage_start + i * stage_len) * oneft # stage_start & stage_len in feet; stage is in meter
        stageft = int(stage_start + i * stage_len) # in feet, integer
        results.extend(floodstage(hand, hfmap, catch, stage, stageft, r_catch.transform))
    
        t01 = time.perf_counter()
        print('TIME', 'stage(stageft)', stage, stageft, (t01-t00), 'seconds', flush=True)
    t1 = time.perf_counter()
    print('TIME', 'allstages', num_stages, (t1-t0), 'seconds', flush=True)
    
    # output
    t0 = time.perf_counter()
    # o_flood_extent_file = ddir + '/' + 'floodstack' + '/' + hucid + '_' + str(stageft) + '.shp'
    o_flood_extent_file = odir + '/' + prefix + '_floodstack_' + str(int(stage_start)) + '_' + str(num_stages) + '_' + str(int(stage_len)) + '.shp'
    with fiona.open(
        o_flood_extent_file, 'w', 
        driver='ESRI Shapefile',
        crs=r_catch.crs,
        schema={'properties': OrderedDict([
            ('HydroID', 'int'),
            ('feature_id', 'int'),
            ('stage_ft', 'int')
        ]), 'geometry': 'Polygon'}) as dst:
        dst.writerecords(results)
    t1 = time.perf_counter()
    print('TIME', 'output', (t1-t0), 'seconds', flush=True)
    # cleanup
    r_catch.close()
    src.close()

def floodstage(hand, hfmap, catch, stage, stageft, t):
    # create a flood extent raster with a stage value
    flood_extent = hand < stage
    # filter by flood extent
    flooded_catch = catch * flood_extent
    # polygonize 
    results = [
        {'properties': {'HydroID': int(v), 'feature_id': int(hfmap[v]), 'stage_ft': stageft}, 'geometry': s} for i, (s, v) in 
        enumerate(rasterio.features.shapes(flooded_catch, mask=flood_extent, transform=t))
        if v > 0
    ]
    
    return results

if __name__ == '__main__':
    hucid = '12090205'
    ddir = '/dev/shm/fim3/' + hucid
    raster_file = ddir + '/' + 'rem_zeroed_masked.tif'
    hfmap_file = ddir + '/' + 'crosswalk_table.csv'
    catch_file = ddir + '/' + 'gw_catchments_reaches_filtered_addedAttributes.tif'
    #floodstack(f_hand=raster_file, f_hfmap=hfmap_file, f_catch=catch_file, stage_start=1.0, num_stages=3, stage_len=1.0, odir=ddir + '/floodstack', prefix=hucid)
    floodstack(f_hand=raster_file, f_hfmap=hfmap_file, f_catch=catch_file, stage_start=1.0, num_stages=80, stage_len=1.0, odir=ddir + '/floodstack', prefix=hucid)
