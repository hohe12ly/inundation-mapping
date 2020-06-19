#!/usr/bin/env python3

import numpy as np
from numba import njit, typeof, typed, types
import concurrent.features 
import rasterio
import argparse
import json


def inundate(remFileName,catchmentsFileName,forecast_fileName,src_fileName,cross_walk_table_fileName,
                inundation_raster_fileName=None,inundation_polygon_fileName=None,depths_fileName=None,stages_fileName=None,num_workers=4):

    # make a catchment,stages numba dictionary
    catchmentStagesDict = __make_catchment_stages_dictionary(forecast_fileName,src_fileName,cross_walk_table_fileName)
    
    rem = rasterio.open(remFileName)
    catchments = rasterio.open(catchmentsFileName)

    profile = rem.profile
    #profile.update(blockxsize=256, blockysize=256, tiled=True)

    depths = rasterio.open(outfile, "w", **profile):

    # Materialize a list of destination block windows
    # that we will use in several statements below.
    windows = [window for ij, window in depths.block_windows()]

    # This generator comprehension gives us raster data
    # arrays for each window. Later we will zip a mapping
    # of it with the windows list to get (window, result)
    # pairs.
    rem_gen = (rem.read(window=window) for window in windows)
    catchments_gen = (catchments.read(window=window) for window in windows)
    depths_gen = (depths.read(window=window) for window in windows)

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=num_workers)
    for window, result in zip(windows, executor.map(__compute_depths, rem_gen,catchments_gen,depths_gen)):
        depths.write(result, window=window)

    executor.done()
    rem.close()
    catchments.close()
    depths.close()

def __compute_depths(rem,catchments,catchmentStagesDict,depths):
    
    shape = rem.shape
    # flatten
    rem = rem.ravel()
    catchments = catchments.ravel()
    depths = rem.copy()
    
    depths = __make_depths_grid(rem,catchments,catchmentStagesDict,depths)

    depths = depths.reshape(shape)

    return(depths)

@njit
def __make_depths_grid(rem,catchments,catchmentStagesDict,depths):

    for i,(r,cm) in enumerate(zip(rem,catchments)):
        if cm in catchmentStagesDict:
            depth = catchmentStagesDict[cm] - r
            depths[i] = max(depth,0)

    return(depths)

def __make_catchment_stages_dictionary(forecast_fileName,src_fileName,cross_walk_table_fileName):
    """ test """

    forecast = pd.read_csv(forecast_fileName, dtype={'feature_id' : int , 'discharge' : float})
    # forecast = forecast.astype({'feature_id' : int , 'discharge' : float})

    with open(src_fileName,'r') as f:
        src = json.load(f)

    cross_walk_table = pd.read_csv(cross_walk_table_fileName, dtype={'feature_id' : int , 'HydroID' : int})

    # hydroIDs = np.unique(list(src.keys()))[0]
    catchmentStagesDict = typed.Dict.empty(types.int32,types.float64)

    number_of_forecast_points = len(forecast)

    for _,rows in tqdm(forecast.iterrows(),total=number_of_forecast_points):
        discharge = rows['discharge']
        fid = int(rows['feature_id'])

        # discharge = rows[1]
        # fid = rows[0]
        matching_hydroIDs = cross_walk_table['HydroID'][cross_walk_table['feature_id'] == fid]

        for hid in matching_hydroIDs:

            stage_list = np.array(src[str(hid)]['stage_list'])
            q_list = np.array(src[str(hid)]['q_list'])
            indices_that_are_lower = list(q_list < discharge)

            # print(indices_that_are_lower)
            is_index_last = indices_that_are_lower[-1]

            if is_index_last:
                h = stage_list[-1]

                hid = types.int32(hid) ; h = types.float32(h)
                catchmentStagesDict[hid] = h

                continue

            index_of_lower = np.where(indices_that_are_lower)[0][-1]
            index_of_upper = index_of_lower + 1

            Q_lower = q_list[index_of_lower]
            h_lower = stage_list[index_of_lower]

            Q_upper = q_list[index_of_upper]
            h_upper = stage_list[index_of_upper]

            # linear interpolation
            h = h_lower + (discharge - Q_lower) * ((h_upper - h_lower) / (Q_upper - Q_lower))

            hid = types.int32(hid) ; h = types.float32(h)
            catchmentStagesDict[hid] = h

    return(catchmentStagesDict)


if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser(description='Relative elevation from pixel based watersheds')
    parser.add_argument('-r','--rem', help='DEM to use within project path', required=True)
    parser.add_argument('-c','--catchments',help='Basins polygons to use within project path',required=True)
    parser.add_argument('-f','--forecast',help='Discharges CSV file',required=True)
    parser.add_argument('-s','--src',help='SRC CSV file',required=True)
    parser.add_argument('-i','--inundation-raster',help='Inundation Raster',required=False,default=None)
    parser.add_argument('-p','--inundation-polygon',help='Inundation polygon',required=False,default=None)
    parser.add_argument('-d','--depths',help='Depths raster',required=False,default=None)
    parser.add_argument('-g','--stages',help='Stages raster',required=False,default=None)
    parser.add_argument('-w','--crosswalk-table',help='Cross-walk table csv',required=False,default=None)
    # parser.add_argument('-d','--catchment-stages',help='Pixel based watersheds raster to use within project path',required=True)

    # extract to dictionary
    args = vars(parser.parse_args())

    # rename variable inputs
    remFileName = args['rem']
    catchmentsFileName = args['catchments']
    forecast_fileName = args['forecast']
    src_fileName = args['src']
    inundation_raster_fileName = args['inundation_raster']
    inundation_polygon_fileName = args['inundation_polygon']
    depths_fileName = args['depths']
    stages_fileName = args['stages']
    cross_walk_table_fileName = args['crosswalk_table']
    # catchmentStageDictFileName = args['catchment_stages']

    inundate(remFileName,catchmentsFileName,forecast_fileName,src_fileName,cross_walk_table_fileName,
             inundation_raster_fileName,inundation_polygon_fileName,
             depths_fileName,stages_fileName)
