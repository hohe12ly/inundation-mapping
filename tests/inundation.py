#!/usr/bin/env python3

import numpy as np
import pandas as pd
from numba import njit, typeof, typed, types
import argparse
from raster import Raster
import gc
from tqdm import tqdm
from subprocess import run
from os import remove
from os.path import splitext,isfile,basename
import json


def inundateREM(remFileName,catchmentsFileName,forecast_fileName,src_fileName,cross_walk_table_fileName,
                inundation_raster_fileName=None,inundation_polygon_fileName=None,depths_fileName=None,stages_fileName=None):

    # stageGrid and flat catchments
    catchments = Raster(catchmentsFileName)
    catchmentsBoolean = catchments.array != catchments.ndv
    flat_catchments = catchments.array[catchmentsBoolean].ravel()
    flat_stage_grid = flat_catchments.copy()
    flat_stage_grid = flat_stage_grid.astype(np.float32)
    flat_stage_grid[:] = 0

    del catchments
    gc.collect()
    
    # make a catchment,stages numba dictionary
    catchmentStagesDict = __make_catchment_stages_dictionary(forecast_fileName,src_fileName,cross_walk_table_fileName)

    # make a flat stages grid
    flat_stage_grid = __make_stage_grid(flat_stage_grid,catchmentStagesDict,flat_catchments)

    del flat_catchments, catchmentStagesDict
    gc.collect()

    # build stages raster
    rem = Raster(remFileName)
    stage_grid = rem.copy()
    stage_grid.array = stage_grid.array.astype(np.float32)
    stage_grid.ndv = -9999
    stage_grid.array[:] = stage_grid.ndv
    stage_grid.array[catchmentsBoolean] = flat_stage_grid

    del flat_stage_grid
    gc.collect()

    if stages_fileName is not None:
        stage_grid.writeRaster(stages_fileName)

    print("Make depths grid")
    depths_grid = stage_grid.copy()

    depths_grid.array[catchmentsBoolean] = stage_grid.array[catchmentsBoolean] - rem.array[catchmentsBoolean]

    del rem, stage_grid
    gc.collect()

    print("Make inundation grid")
    inundation_grid = depths_grid.copy()
    inundation_grid.array = inundation_grid.array.astype(np.int32)
    inundation_grid.ndv = 0
    inundation_grid.array[:] = inundation_grid.ndv
    inundatedLocations = np.logical_and(catchmentsBoolean,depths_grid.array > 0)
    nonInundatedLocations = np.logical_and(catchmentsBoolean,depths_grid.array <= 0)
    inundation_grid.array[nonInundatedLocations] = 1
    inundation_grid.array[inundatedLocations] = 2

    if inundation_raster_fileName is None:
        inundation_raster_fileName = "inundation_temp.tif"
    inundation_grid.writeRaster(inundation_raster_fileName)
    

    if inundation_polygon_fileName is not None:
        print("Polygonize inundation raster")

        driverDictionary = {'.gpkg' : 'GPKG','.geojson' : 'GeoJSON','.shp' : 'ESRI Shapefile'}
        driver = driverDictionary[splitext(inundation_polygon_fileName)[1]]

        inundation_polygon_layerName = splitext(basename(inundation_polygon_fileName))[0]

        if isfile(inundation_polygon_fileName):
            remove(inundation_polygon_fileName)

        calcCommand = ['gdal_calc.py','--type=Int32','--overwrite','-A',inundation_raster_fileName,'-B',catchmentsFileName,'--calc=\"B*(A>1)\"','--NoDataValue=0','--outfile=\"inundation_catchments.tif\"']
        calcCommand = 'gdal_calc.py --type=Int32 --overwrite -A {} -B {} --calc="B*(A>1)" --NoDataValue=0 --outfile="inundation_catchments.tif"'.format(inundation_raster_fileName,catchmentsFileName)
        polyCommand = ['gdal_polygonize.py','-8','-f',driver,'inundation_catchments.tif',inundation_polygon_fileName,inundation_polygon_layerName,'HydroID']

        run(calcCommand,shell=True)
        run(polyCommand)

        remove("inundation_catchments.tif")
        if inundation_raster_fileName == 'inundation_temp.tif':
            remove(inundation_raster_fileName)

    del inundatedLocations, nonInundatedLocations, inundation_grid
    gc.collect()

    if depths_fileName is not None:
        print("Zero out depths grid")
        negativeBoolean = np.logical_and(catchmentsBoolean,depths_grid.array < 0)
        depths_grid.array[negativeBoolean] = depths_grid.ndv
        depths_grid.writeRaster(depths_fileName)

    # flat_stage_grid = flat_stage_grid.reshape(catchments.nrows,catchments.ncols)

@njit
def __make_stage_grid(flat_stage_grid,catchmentStagesDict,flat_catchments):

    for i,cm in enumerate(flat_catchments):
        if cm in catchmentStagesDict:
            flat_stage_grid[i] = catchmentStagesDict[cm]


    return(flat_stage_grid)

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

    inundateREM(remFileName,catchmentsFileName,forecast_fileName,src_fileName,cross_walk_table_fileName,
                inundation_raster_fileName,inundation_polygon_fileName,
                depths_fileName,stages_fileName)
