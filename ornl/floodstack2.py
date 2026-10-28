'''
This module generates a flood stack for a given HUC unit. 

This version implements the following workflow steps:
1. exclude lake catchments, by using the `LakeID` attribute in the rating curve table. If `LakeID != -999`, then the catchment is a lake.
2. generate the depth grid first, by calling my TauDEM command `handinun` in python. `handinun` is C++ code. It generates depth grids for all the stage heights specified in the input stage_id:stage_height table.
3. convert depth grids to binary extent grids. If depth > 0, then the extent grid value is 1, else 0.
4. polygonize the extent grids. For each extent vector, fix invalid polygons which often have queen-type connectivity. The fix uses shapely's `is_valid` and `make_valid` functions. The fixed polygons are written to two formats: a POLYGON geopackage and a MULTIPOLYGON geopackage. A single geopackage is created to include extent vectors for all the stage heights. So only two geopackages are created.

'''

import argparse
import numpy as np
import numpy.ma as ma
import pickle
import matplotlib.pyplot as plt

import os, sys, time, subprocess
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
import json
import pyproj

import rtree
import math
import random

import csv
from collections import OrderedDict

from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp


def fix_polygons(g, sieve_area):
    '''
    fix invalid polygons and return a list of valid polygons. 
    also remove small polygons with area <= sieve_area. area unit must be the same as data projection
    '''
    g_list = []
    if not shapely.is_valid(g):
        g_fix = shapely.make_valid(g)
        g_fix_list = []
        if type(g_fix) == shapely.MultiPolygon:
            for gg in g_fix.geoms:
                g_fix_list.append(gg)
        elif type(g_fix) == shapely.Polygon:
            g_fix_list.append(g_fix)
        elif type(g_fix) == shapely.GeometryCollection:
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
    else: # always return a list of polygons
        if type(g) == shapely.MultiPolygon or type(g) == shapely.GeometryCollection:
            for gg in g.geoms:
                g_list.append(gg)
        else:
            g_list.append(g)
    # remove small polygons 
    return [ g for g in g_list if g.area > sieve_area ]

# polygonization
def extent_polygonize(catch_b, depth_mask, catch_crs_epsg, catch_transform, hydroid2comid, flowlines, flowline_tree, discharge_dict, odir, oepsg, tdir, huc, stage_index, stage_height, sieve_area=10.0, simplify_tolerance=3.0, max_segment_length=4.0, output_singlepolygon=False):

    # create reprojection transformer from catchment grid to output epsg
    transformer = pyproj.Transformer.from_crs(
        pyproj.CRS('EPSG:'+str(catch_crs_epsg)),
        pyproj.CRS('EPSG:'+str(oepsg)),
        always_xy=True).transform
        
    t0 = time.perf_counter()
    #depth[~depth_mask] = catch_b # depth grid --> extent grid
    extent_shp = rasterio.features.shapes(catch_b,mask= ~depth_mask,connectivity=8,transform=catch_transform)

    # post processing rules:
    # 1. since we use queen connectivity, we need to fix invalid polygons.
    # 2. we merge polygons into a multipolygon. if the multipolygon area is less than 10 square meters, we ignore it.
    # 3. we check if this multipolygon intersects with any flowline. if not, we ignore it. this is to remove isolated flood patches.
    # 4. we use Douglas-Peucker algorithm to simplify the multipolygon. we use a tolerance of 3 meters.
    polygons_by_hydroid = {}
    count_invalid_hydro = 0
    count_invalid_geom = 0
    count_invalid_geom_fixed = 0
    count = 0
    count_in = 0
    count_polygons = 0
    count_multipolys = 0
    for i,(geom, hydroid) in enumerate(extent_shp):
        count_in += 1
        if hydroid <= 0:
            count_invalid_hydro += 1
            continue
        g = shapely.geometry.shape(geom)
        # make polygons valid
        count_polygons += 1
        g_list = fix_polygons(g, sieve_area=0.0) # we will filter out small polygons later
        if len(g_list) == 0 :
            continue
        # merge polygons at hydroid catchment level
        mp = shapely.MultiPolygon(g_list)
        mp = mp.buffer(0) # merge adjacent polygons, returns a multipolygon or a polygon
        if type(mp) == shapely.Polygon:
            mp = shapely.MultiPolygon([mp])
        # use rtree to check if this multipolygon intersects with any flowline. if not, we ignore it. this is to remove isolated flood patches.
        intersects = False
        for j in flowline_tree.intersection(mp.bounds):
            if mp.intersects(flowlines[j]):
                intersects = True
                break
        if not intersects:
            continue
        
        # simplify the multipolygon using Douglas-Peucker algorithm, 
        mp = mp.simplify(simplify_tolerance, preserve_topology=True)
        # fix again, to be sure
        g_list = fix_polygons(mp, sieve_area)
        if len(g_list) == 0:
            continue
        mp = shapely.MultiPolygon(g_list)
        mp = mp.buffer(0) # merge adjacent polygons, returns a multipolygon or a polygon
        if type(mp) == shapely.Polygon:
            mp = shapely.MultiPolygon([mp])

        if hydroid not in polygons_by_hydroid:
            polygons_by_hydroid[hydroid] = []
        polygons_by_hydroid[hydroid].append(mp)


    # if a hydroid has multiple multipolygons, we merge them into a single multipolygon
    mp_by_hydroid = {}
    for hydroid, mp_list in polygons_by_hydroid.items():
        if len(mp_list) == 1:
            mp_by_hydroid[hydroid] = mp_list[0]
        else:
            # extract all polygons from the multipolygon list and merge them into a single multipolygon
            all_polygons = []
            for mp in mp_list:
                for poly in mp.geoms:
                    all_polygons.append(poly)
            mp = shapely.MultiPolygon(all_polygons)
            mp = mp.buffer(0) # merge adjacent polygons, returns a multipolygon or a polygon
            if type(mp) == shapely.Polygon:
                mp = shapely.MultiPolygon([mp])
            mp_by_hydroid[hydroid] = mp
            print(f'POLYGONIZE stageindex_{stage_index} hydroid {hydroid} has {len(mp_list)} multipolygons, merged into a single multipolygon with {len(mp.geoms)} polygons.')

    rows = []
    for hydroid, mp in mp_by_hydroid.items():
        count_multipolys += 1
        # added to output list
        rows.append({
            'geometry': json.loads(shapely.to_geojson(mp)),
            'properties': OrderedDict([
                ('ID', int(count)),
                ('nextgen_id', int(hydroid)),
                ('flow', int(discharge_dict[hydroid])),
                ('COMID', int(hydroid2comid[hydroid])),
                ('StageIndex', int(stage_index)),
                ('StageHeight', float(stage_height))
            ])
        })
        count += 1
    
    extent_meta = {
        'crs' : fiona.crs.from_epsg(catch_crs_epsg), # hhand.crs , 
        # 'crs' : fiona.crs.from_epsg(oepsg), # output epsg
        'driver' : 'GPKG', # to avoid int32 limitation in shapefile
        'layer' : 'floodstack_'+huc+'_stageindex_'+str(stage_index),
        'schema': {
            'geometry' : 'MultiPolygon',
            'properties' : OrderedDict([
                ('ID' , 'int'),
                ('nextgen_id' , 'int'),
                ('flow' , 'int'),
                ('COMID' , 'int'),
                ('StageIndex' , 'int'),
                ('StageHeight' , 'float')
            ])
        }
    }
    # extent_fname = odir +'/extent_vector/'+'extent_'+huc+'_stageindex_'+str(stage_index)+f'_EPSG_{catch_crs_epsg}.gpkg'
    extent_fname = odir +'/extent_vector/'+'floodstack_'+huc+f'_EPSG_{catch_crs_epsg}.gpkg'
    with fiona.open(extent_fname,'w',**extent_meta) as ovector: # this adds a new layer to the same geopackage file
        ovector.writerecords(rows)

    if output_singlepolygon:
        # now output a Polygon geopackage, by splitting the multipolygons into polygons and copying attributes. we will recreate IDs
        extent_meta['schema']['geometry'] = 'Polygon'
        extent_fname = odir +'/extent_vector/'+'extent_'+huc+'_stageindex_'+str(stage_index)+f'_EPSG_{catch_crs_epsg}_polygon.gpkg'
        sp_count = 0
        with fiona.open(extent_fname,'w',**extent_meta) as ovector:
            for row in rows:
                mp = shapely.geometry.shape(row['geometry'])
                for poly in mp.geoms:
                    ovector.write({
                        'geometry': json.loads(shapely.to_geojson(poly)),
                        'properties': OrderedDict([
                            ('ID', int(sp_count)),
                            ('nextgen_id', int(row['properties']['nextgen_id'])),
                            ('flow', int(row['properties']['flow'])),
                            ('COMID', int(row['properties']['COMID'])),
                            ('StageIndex', int(row['properties']['StageIndex'])),
                            ('StageHeight', float(row['properties']['StageHeight']))
                        ])
                    })
                    sp_count += 1

    # reproject to output epsg:
    omp_by_hydroid = {}
    for hydroid, mp in mp_by_hydroid.items():
        # segmentize the multipolygon to preserve shape during reprojection. we use a segment length of 4 meters.
        mp = mp.segmentize(max_segment_length=4.0) # assumption: this step should not produce invalid polygons.
        # reporject the multipolygon to output epsg, which is geographic coordinates.
        mp = shapely.ops.transform(transformer, mp)
        # simplify the multipolygon using Douglas-Peucker algorithm, we use a tolerance that is the same as parameter simplify_tolerance, which is 3 meters. we need to convert the tolerance to degrees. 1 degree = 111319.9 meters at the equator. so tolerance in degrees = 3 / 111319.9
        tolerance_degrees = simplify_tolerance / 111319.9
        mp = mp.simplify(tolerance_degrees, preserve_topology=True)
        # let's fix the multipolygon again, to be sure.
        # we still use the sieve_area parameter, but need to convert from square meters to square degrees. 1 square degree = (111319.9)^2 square meters at the equator. so sieve_area in square degrees = sieve_area / (111319.9^2)
        sieve_area_degrees = sieve_area / (111319.9**2)
        g_list = fix_polygons(mp, sieve_area=sieve_area_degrees)
        if len(g_list) == 0:
            print(f'POLYGONIZE stageindex_{stage_index} hydroid {hydroid} has no valid polygons after reprojection and simplification. Ignoring.')
            continue
        mp = shapely.MultiPolygon(g_list)
        mp = mp.buffer(0) # merge adjacent polygons, returns a multipolygon or a polygon
        if type(mp) == shapely.Polygon:
            mp = shapely.MultiPolygon([mp])

        omp_by_hydroid[hydroid] = mp
    
    orows = []
    for hydroid, mp in omp_by_hydroid.items():
        orows.append({
            'geometry': json.loads(shapely.to_geojson(mp)),
            'properties': OrderedDict([
                ('ID', int(count)),
                ('nextgen_id', int(hydroid)),
                ('flow', int(discharge_dict[hydroid])),
                ('COMID', int(hydroid2comid[hydroid])),
                ('StageIndex', int(stage_index)),
                ('StageHeight', float(stage_height))
            ])
        })

    # change crs of the output to the output epsg
    extent_meta['crs'] = fiona.crs.from_epsg(oepsg)
    extent_fname = odir +'/extent_vector/'+'floodstack_'+huc+f'_EPSG_{oepsg}.gpkg'
    with fiona.open(extent_fname,'w',**extent_meta) as ovector: # this adds a new layer to the same geopackage file
        ovector.writerecords(orows)

    print(f'POLYGONIZE stageindex_{stage_index} total polygons/multipolygons after polygonization:', count_in)
    print(f'POLYGONIZE stageindex_{stage_index} total polygons:', count_polygons)
    if output_singlepolygon:
        print(f'POLYGONIZE stageindex_{stage_index} total output SinglePolygon features:', sp_count)
    print(f'POLYGONIZE stageindex_{stage_index} total output generated MultiPolygons:', count)
    print(f'POLYGONIZE stageindex_{stage_index} total output MultiPolygon features:', count_multipolys)
    print(f'POLYGONIZE stageindex_{stage_index} total output stage index:', stage_index)
    print(f'POLYGONIZE stageindex_{stage_index} total output stage height:', stage_height)
    sys.stdout.flush()

    t1 = time.perf_counter()
    print('TIME', 'polygonization.stageindex_'+str(stage_index), 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
    sys.stdout.flush()

def generate_stage_height_table(stages_min = 0, stages_max = 25, stages_interval = 0.3048):
    """
    Generate a table of stage heights in meters.
    """
    stages_max = stages_max + stages_interval
    stages = np.round(np.arange(stages_min,stages_max,stages_interval),4)

    return stages

# python floodstack2.py 12110111                     # serial
# python floodstack2.py 12110111 --parallel --num_workers 8
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Flood inundation map stack generation for a given HUC unit.')
    parser.add_argument('huc', help='HUC code')
    parser.add_argument('--parallel', action='store_true')
    parser.add_argument('--num_workers', type=int, default=8)
    args = parser.parse_args()

    parallel_processing = args.parallel
    num_workers = args.num_workers
    print("PARALLEL PROCESSING:", parallel_processing, 'NUM_WORKERS:', num_workers if parallel_processing else 1)

    t00 = time.perf_counter()

    huc = args.huc
    handdir = f'/mnt/CADES_CEPH/pin2flood/prod/{huc}'
    if not os.path.exists(handdir):
        print('ERROR: handinun output directory does not exist:', handdir)
        sys.exit(1)
    odir = f'/mnt/CADES_CEPH/pin2flood/floodstacks/{huc}'
    for subdir in ['depth_grid','extent_grid','extent_vector']:
        os.makedirs(os.path.join(odir, subdir), exist_ok=True)
    tdir = f'/mnt/CADES_CEPH/pin2flood/cache/floodstacks/{huc}'
    os.makedirs(tdir, exist_ok=True)

    f_braided_lookup = '/mnt/CADES_CEPH/fim3/cahaba/ornl/braided_streams/braided_streams_lookup_minor2main.pickle' # not used.
    # f_hhand = os.path.join(handdir, 'rem_zeroed_masked_healed.tif') # original HAND raster
    f_hhand = os.path.join(odir, f'{huc}_surface_healed_hand.tif') # use this hand, which adds the DSM-DEM diff to the healed HAND. This is the final HAND raster to use for flood stack generation.
    f_catchgrid = os.path.join(handdir, 'gw_catchments_reaches_filtered_addedAttributes.tif') # catchment raster
    f_hydrotable = os.path.join(handdir, 'hydroTable_rp_bf_lmtdischarge_cda.csv') # rating curve table
    f_flowlines = os.path.join(handdir, 'demDerived_reaches_split_filtered_addedAttributes_crosswalked.gpkg') # flowline shapefile

    ## how does Cahaba create the stage height list? refer to `notebooks/stage_table.ipynb`
    # in `params_calibrated.env`: min=0m; max=25m; interval=0.3048m (1ft). 
    # in `make_stages_and_catchlist.py`: 
    #.   stages_max = stages_max + stages_interval
    #.   stages = np.round(np.arange(stages_min,stages_max,stages_interval),4)
    #    for stage in stages:
    #       f.write("{}\n".format(stage))

    # read flowlines using fiona. we only need the linestring geometry.
    flowlines = []
    with fiona.open(f_flowlines, 'r') as src:
        for feature in src:
            geom = shapely.geometry.shape(feature['geometry'])
            flowlines.append(geom)
    # build rtree index for the flowlines in the tmp dir.
    flowline_tree = rtree.index.Index()
    for i, line in enumerate(flowlines):
        flowline_tree.insert(i, line.bounds)

    t0 = time.perf_counter()
    d_hydrotable = pd.read_csv(f_hydrotable, usecols=['FATSGTID', 'feature_id', 'stage', 'discharge_cms', 'LakeID'])
    d_hydrotable.rename(columns={'FATSGTID': 'HydroID', 'feature_id': 'COMID'}, inplace = True)
    comids_h = np.unique(d_hydrotable['COMID'].to_numpy())
    comids = comids_h
    d_data = d_hydrotable
    print('num of unique COMIDs in combined data:', comids.size)
    t1 = time.perf_counter()
    print('TIME', 'load_hydrodata', 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
    sys.stdout.flush()
    
    hydroids = np.unique(d_data['HydroID'].to_numpy())
    print('num of unique HydroIDs:', hydroids.size)
    sys.stdout.flush()
    
    hydroid2comid = {}
    for hydroid in hydroids:
        comid = d_data[d_data['HydroID'] == hydroid].iloc[[0]]['COMID'].values[0]
        hydroid2comid[hydroid] = comid

    # group hydrotable by HydroID, each group should have 84 rows, corresponding to 84 stage heights. If not, we will abort.
    d_data_by_hydroid = d_data.groupby('HydroID')
    for hydroid, group in d_data_by_hydroid:
        if group.shape[0] != 84:
            print('ERROR: HydroID', hydroid, 'has', group.shape[0], 'rows in hydrotable, expected 84 rows. Aborting.')
            sys.exit(1)
    
    # internal function to get discharge of each hydroid at a given stage height, which is just the row index in d_data_by_hydroid.
    def get_discharge(stage_index):
        discharge_dict = {}
        for hydroid, group in d_data_by_hydroid:
            discharge = group.iloc[[stage_index]]['discharge_cms'].values[0]
            discharge_dict[hydroid] = discharge
        return discharge_dict

    # create a list of hydroids which are lakes, i.e., LakeID != -999. These hydroids will be excluded from the flood stack generation.
    lake_hydroids = d_data[d_data['LakeID'] != -999]['HydroID'].unique()
    print('num of lake hydroids:', lake_hydroids.size)
    
    ## load catchment grid
    t0 = time.perf_counter()
    catchgrid = rasterio.open(f_catchgrid)
    catch_b = catchgrid.read(1)
    print('input: catchment raster:', catchgrid.height, 'rows x', catchgrid.width, 'cols', 'nodata:', catchgrid.nodata)
    t1 = time.perf_counter()
    print('TIME', 'load_catchgrid', 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
    sys.stdout.flush()
    # create a catchment mask. the value of catch_b is the HydroID of the catchment. mask those lake catchments.
    catch_mask = np.isin(catch_b, lake_hydroids)
    print('catchment mask: num of lake catchment cells:', np.count_nonzero(catch_mask))

    ## load Healed HAND
    t0 = time.perf_counter()
    hhand = rasterio.open(f_hhand)
    # b = hhand.read(1) # if using explicit nodata mask
    b_hhand = hhand.read(1, masked=True) # as numpy MaskedArray
    print('input: Healed HAND raster:', hhand.height, 'rows x', hhand.width, 'cols', 'nodata:', hhand.nodata)
    t1 = time.perf_counter()
    print('TIME', 'load_healed_hand', 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
    sys.stdout.flush()
    data_mask = b_hhand.mask

    inundation_threshold = 0.1 # 10cm as used by Shih-Chieh's group. In HAND, we use 1mm

    # assert that catchment grid and HAND grid have the same shape
    if catchgrid.height != hhand.height or catchgrid.width != hhand.width:
        print('ERROR: catchment grid and HAND grid have different shapes:', catchgrid.height, 'x', catchgrid.width, 'vs', hhand.height, 'x', hhand.width)
        sys.exit(1)

    stages = generate_stage_height_table(stages_min = 0, stages_max = 25, stages_interval = 0.3048)
    depth_meta = hhand.meta.copy()
    depth_meta['compress'] = 'lzw'
    hhand_nodata = hhand.nodata
    catch_crs_epsg = catchgrid.crs.to_epsg()
    catch_transform = catchgrid.transform

    def _serial_processing():
        for stage_index, stage_height in enumerate(stages):
            # if stage_index < len(stages) - 1:
            #     continue # debugging: only process the last stage height for now
            print('===STAGE HEIGHT', stage_index, 'stage_height_meters', stage_height, '===')
            t0 = time.perf_counter()
            depth = stage_height - b_hhand # depth grid
            # we only keep those cells with depth values > inundation_threshold, and not in lake catchments, and not nodata in HAND. all other cells are masked out.
            depth_mask = (depth <= inundation_threshold) | catch_mask | data_mask
            depth[depth_mask] = hhand_nodata
            t1 = time.perf_counter()
            print('TIME', 'depth_grid.stageindex_'+str(stage_index), 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
            sys.stdout.flush()

            # flood depth grid output
            t0 = time.perf_counter()
            depth_fname = odir +'/depth_grid/'+'depth_'+huc+'_stageindex_'+str(stage_index)+'.tif'
            with rasterio.open(depth_fname, 'w', **depth_meta, BIGTIFF='YES') as oraster:
                oraster.write(depth, indexes=1)
            t1 = time.perf_counter()
            print('TIME', 'output_depth_grid.stageindex_'+str(stage_index), 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
            sys.stdout.flush()

            # depth_mask is also the binary extent grid: True - not inundated or excluded; False - inundated. We can use depth_mask to polygonize the extent grid.

            t0 = time.perf_counter()
            discharge_dict = get_discharge(stage_index)
            extent_polygonize(catch_b, depth_mask, catch_crs_epsg, catch_transform, hydroid2comid, flowlines, flowline_tree, discharge_dict, odir, 4326, tdir, huc, stage_index, stage_height, sieve_area = 10.0, simplify_tolerance = 3.0)
            t1 = time.perf_counter()
            print('TIME', 'extent_polygonization.stageindex_'+str(stage_index), 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
            sys.stdout.flush()

    def _process_stage(stage_index, stage_height):
        print('===STAGE HEIGHT', stage_index, 'stage_height_meters', stage_height, '===')
        t0 = time.perf_counter()
        depth = stage_height - b_hhand
        depth_mask = (depth <= inundation_threshold) | catch_mask | data_mask
        depth[depth_mask] = hhand_nodata
        t1 = time.perf_counter()
        print('TIME', 'depth_grid.stageindex_'+str(stage_index), 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
        sys.stdout.flush()

        t0 = time.perf_counter()
        depth_fname = odir +'/depth_grid/'+'depth_'+huc+'_stageindex_'+str(stage_index)+'.tif'
        with rasterio.open(depth_fname, 'w', **depth_meta, BIGTIFF='YES') as oraster:
            oraster.write(depth, indexes=1)
        t1 = time.perf_counter()
        print('TIME', 'output_depth_grid.stageindex_'+str(stage_index), 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
        sys.stdout.flush()

        t0 = time.perf_counter()
        discharge_dict = get_discharge(stage_index)
        extent_polygonize(catch_b, depth_mask, catch_crs_epsg, catch_transform, hydroid2comid, flowlines, flowline_tree, discharge_dict, odir, 4326, tdir, huc, stage_index, stage_height, sieve_area=10.0, simplify_tolerance=3.0)
        t1 = time.perf_counter()
        print('TIME', 'extent_polygonization.stageindex_'+str(stage_index), 'time_in_seconds', "{:.2f}".format(t1-t0), flush=True)
        sys.stdout.flush()

    if not parallel_processing:
        _serial_processing()
    else:
        # parallel processing using ProcessPoolExecutor
        # reverse the stages list to process larger stage heights first.
        stages = stages[::-1]
        ctx = mp.get_context('fork')  # Use 'fork' start method for Unix/Linux-like systems
        with ProcessPoolExecutor(max_workers=num_workers, mp_context=ctx) as executor:
            futures = []
            for i, stage_height in enumerate(stages):
                stage_index = len(stages) - 1 - i  # reverse index to match the original order
                futures.append(executor.submit(_process_stage, stage_index, stage_height))
            for future in futures:
                future.result()  # wait for all tasks to complete

    # cleanup work
    flowline_tree.close()
    # close opened input raster files
    catchgrid.close()
    hhand.close()

    t11 = time.perf_counter()
    print('TIME', 'TOTAL_TIME', 'time_in_seconds', "{:.2f}".format(t11-t00), flush=True)
    sys.stdout.flush()