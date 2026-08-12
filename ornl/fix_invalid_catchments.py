from skimage import data
from skimage import measure

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import os, sys, time, glob
import pandas as pd
import shapely
import shapely.geometry
import shapely.errors
import shapely.ops
import fiona
import fiona.crs
import rasterio
import rasterio.mask
import rasterio.features
import rasterio.windows
import affine
#import geopandas as gpd

import rtree
import math

import csv

def fix_catch(f):
    g = shapely.geometry.shape(f["geometry"])
    g_list = []
    use_type = 'None'
    source_valid = True
    o_valid = True
    if not shapely.is_valid(g):
        source_valid = False
        g_fix = shapely.make_valid(g)
        g_fix_list = []
        if type(g_fix) == shapely.MultiPolygon:
            use_type = 'MultiPolygon'
            for gg in g_fix.geoms:
                g_fix_list.append(gg)
        elif type(g_fix) == shapely.Polygon:
            use_type = "Polygon"
            g_fix_list.append(g_fix)
        elif type(g_fix) == shapely.GeometryCollection:
            use_type = "GeometryCollection"
            for gg in g_fix.geoms:
                if type(gg) == shapely.Polygon:
                    g_fix_list.append(gg)
                elif type(gg) == shapely.MultiPolygon:
                    for ggg in gg.geoms:
                        g_fix_list.append(ggg)
        # add valid fixed polygons
        for gg in g_fix_list:
            if shapely.is_valid(gg):
                g_list.append(gg)
            else:
                o_valid = False
        if not shapely.is_valid(shapely.geometry.MultiPolygon(g_list)):
            o_valid = False
        # fixed, hydroid, orig_type, o_num_polys, o_split_polys, o_valid
        #print('fixed,', f['properties']['HydroID'], ',', use_type, ',', len(g_list), ',', len(g_fix_list), ',', o_valid)
    else:
        g_list.append(g)
    
    return g_list, use_type, source_valid, o_valid

txcatch = fiona.open('/mnt/CADES_CEPH/pin2flood/catchments_texas.gpkg', mode='r', layer='catchments_texas')

o_fixed_schema = txcatch.schema.copy()
o_fixed_schema['geometry'] = 'MultiPolygon'
o_fixed = fiona.open('/mnt/CADES_CEPH/pin2flood/catchments_texas_fixed.gpkg', 
                     mode='w', 
                     layer='catchments_texas_fixed', 
                     driver = txcatch.driver,
                     crs = txcatch.crs,
                     schema = o_fixed_schema
                    )
o_invalid = fiona.open('/mnt/CADES_CEPH/pin2flood/catchments_texas_invalid.gpkg', 
                     mode='w', 
                     layer='catchments_texas_invalid', 
                     driver = txcatch.driver,
                     crs = txcatch.crs,
                     schema = txcatch.schema
                    )
fixed_count = 0
count = 0
logs = []
for catch in txcatch:
    count += 1
    if count % 2000 == 0:
        print(count, "...")
    g_list, use_type, source_valid, o_valid = fix_catch(catch)
    if not source_valid:
        o_invalid.write(catch) # record invalid original
        try:
            geom = shapely.geometry.MultiPolygon(g_list) #.buffer(0)
            fixed = {
                'geometry': shapely.geometry.mapping(geom),
                'properties': catch['properties']
            }
            o_fixed.write(fixed)
            fixed_count += 1
        except Exception as err:
            print('ERROR', 'hydroid', catch['properties']['HydroID'], "len(g_list):", len(g_list), "valid(geom):", shapely.is_valid(geom))
            break
        logs.append([catch['properties']['HydroID'], use_type, len(g_list), o_valid])
print('fixed', fixed_count, 'out of', count, 'catchments')

o_fixed.close()
o_invalid.close()
txcatch.close()

df_logs = pd.DataFrame(logs, columns=['HydroID', 'use_geom_type', 'o_num_polys', 'o_valid'])
df_logs.to_csv('/mnt/CADES_CEPH/pin2flood/catchments_texas_fixing.log', index=False)

