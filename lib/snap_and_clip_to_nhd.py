#!/usr/bin/env python3

import geopandas as gpd
from pandas import concat
import argparse
from os.path import splitext
from shapely.strtree import STRtree
from shapely.geometry import Point
from shapely.ops import nearest_points

def subset_vector_layers(hucCode,projection,nwm_headwaters_fileName,nwm_streams_fileName,nwm_lakes_fileName,nwm_catchments_fileName,wbd_fileName,subset_nwm_streams_fileName,subset_nwm_lakes_fileName,subset_nwm_headwaters_fileName,subset_nwm_catchments_fileName,dissolveLinks=False):

    # loading files
    print("Loading files")
    hucUnitLength = len(str(hucCode))

    wbd = gpd.read_file(wbd_fileName)

    # query nhd+HR streams for HUC code
    print("Subsetting NWM Streams for HUC{} {}".format(hucUnitLength,hucCode))
    nwm_streams = gpd.read_file(nwm_streams_fileName)
    nwm_streams = nwm_streams.loc[nwm_streams.intersects(wbd.geometry[0]),:]
    nwm_streams.reset_index(drop=True,inplace=True)
    crossing_nwm_streams = nwm_streams.loc[nwm_streams.crosses(wbd.geometry[0]),:]
    nwm_streams = gpd.overlay(nwm_streams,wbd,how='intersection')
    nwm_streams.reset_index(drop=True,inplace=True)
    nwm_streams.to_file(subset_nwm_streams_fileName,driver=getDriver(subset_nwm_streams_fileName),index=False)
    del nwm_streams
    
    # find intersecting nwm_headwaters and the inflowing streams too
    print("Subsetting NWM Headwaters for HUC{} {}".format(hucUnitLength,hucCode))
    nwm_headwaters = gpd.read_file(nwm_headwaters_fileName)
    nwm_headwaters = nwm_headwaters.loc[nwm_headwaters.intersects(wbd.geometry[0]),:]
    nwm_headwaters.reset_index(drop=True,inplace=True)
    nwm_headwaters.drop(columns='ORIG_FID',inplace=True)
    
    crossing_nwm_streams = crossing_nwm_streams.explode()
    crossing_nwm_streams.reset_index(drop=True,inplace=True)
    
    crossing_nwm_hw_points_dictionary = {'geometry' : [None] * len(crossing_nwm_streams) , 'ID' :  [None] * len(crossing_nwm_streams)}
    for i,r in crossing_nwm_streams.iterrows():
        g = r['geometry'] ; ID = r['ID']
        crossing_nwm_hw_points_dictionary['geometry'][i] = Point(list(zip(*nearest_points(g,wbd.geometry[0])[0].coords.xy))[0])
        crossing_nwm_hw_points_dictionary['ID'][i] = ID

    crossing_nwm_hw_points = gpd.GeoDataFrame(crossing_nwm_hw_points_dictionary,crs=nwm_headwaters.crs,geometry='geometry')
    nwm_headwaters = concat([nwm_headwaters,crossing_nwm_hw_points],ignore_index=True)
    nwm_headwaters.reset_index(drop=True,inplace=True)
    nwm_headwaters.to_file(subset_nwm_headwaters_fileName,driver=getDriver(subset_nwm_headwaters_fileName),index=False)
    del nwm_headwaters, crossing_nwm_streams, crossing_nwm_hw_points_dictionary, crossing_nwm_hw_points
    
    # find intersecting lakes
    print("Subsetting NWM Lakes for HUC{} {}".format(hucUnitLength,hucCode))
    nwm_lakes = gpd.read_file(nwm_lakes_fileName)
    nwm_lakes = nwm_lakes.loc[nwm_lakes.intersects(wbd.geometry[0]),:]
    nwm_lakes.reset_index(drop=True,inplace=True)
    nwm_lakes.to_file(subset_nwm_lakes_fileName,driver=getDriver(subset_nwm_lakes_fileName),index=False)
    del nwm_lakes

    # find intersecting nwm_headwaters
    print("Subsetting NWM Catchments for HUC{} {}".format(hucUnitLength,hucCode))
    nwm_catchments = gpd.read_file(nwm_catchments_fileName)
    nwm_catchments = nwm_catchments.loc[nwm_catchments.intersects(wbd.geometry[0]),:]
    nwm_catchments.reset_index(drop=True,inplace=True)
    nwm_catchments.to_file(subset_nwm_catchments_fileName,driver=getDriver(subset_nwm_headwaters_fileName),index=False)
    del nwm_catchments

def getDriver(fileName):

    driverDictionary = {'.gpkg' : 'GPKG','.geojson' : 'GeoJSON','.shp' : 'ESRI Shapefile'}
    driver = driverDictionary[splitext(fileName)[1]]

    return(driver)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Relative elevation from pixel based watersheds')
    parser.add_argument('-d','--hucCode', help='DEM to use within project path', required=True,type=int)
    parser.add_argument('-p','--projection', help='DEM to use within project path', required=True)
    parser.add_argument('-w','--nwm-headwaters', help='DEM to use within project path', required=True)
    parser.add_argument('-s','--nwm-streams',help='Basins polygons to use within project path',required=True)
    parser.add_argument('-l','--nwm-lakes', help='DEM to use within project path', required=True)
    parser.add_argument('-m','--nwm-catchments', help='DEM to use within project path', required=True)
    parser.add_argument('-u','--wbd',help='Basins polygons to use within project path',required=True)
    parser.add_argument('-c','--subset-streams',help='Basins polygons to use within project path',required=True)
    parser.add_argument('-a','--subset-lakes',help='Basins polygons to use within project path',required=True)
    parser.add_argument('-t','--subset-nwm-headwaters',help='Basins polygons to use within project path',required=True)
    parser.add_argument('-n','--subset-catchments',help='Basins polygons to use within project path',required=True)
    parser.add_argument('-o','--dissolve-links',help='Basins polygons to use within project path',action="store_true",default=False)

    args = vars(parser.parse_args())

    hucCode = args['hucCode']
    projection = args['projection']
    nwm_headwaters_fileName = args['nwm_headwaters']
    nwm_streams_fileName = args['nwm_streams']
    nwm_lakes_fileName = args['nwm_lakes']
    nwm_catchments_fileName = args['nwm_catchments']
    wbd_fileName = args['wbd']
    subset_nwm_streams_fileName = args['subset_streams']
    subset_nwm_lakes_fileName = args['subset_lakes']
    subset_nwm_headwaters_fileName = args['subset_nwm_headwaters']
    subset_nwm_catchments_fileName = args['subset_catchments']
    dissolveLinks = args['dissolve_links']

    subset_vector_layers(hucCode,projection,nwm_headwaters_fileName,nwm_streams_fileName,nwm_lakes_fileName,nwm_catchments_fileName,wbd_fileName,subset_nwm_streams_fileName,subset_nwm_lakes_fileName,subset_nwm_headwaters_fileName,subset_nwm_catchments_fileName,dissolveLinks)
