#!/usr/bin/env python3

import geopandas as gpd
import pandas as pd
import argparse
from os import makedirs
from os.path import join,splitext
from tqdm import tqdm
from shapely.geometry import mapping,shape,box
from shapely.errors import TopologicalError
from shapely.wkb import dumps,loads
import pygeos as pg

def __vprint(message,verbose=False):
    if verbose:
        print(message)

def __gdf_to_pg(gdf):

    gdf = gdf.explode().reset_index(level=1,drop=True)
    wkb_list = [g.wkb for g in gdf.geometry]
    pg_geom = pg.from_wkb(wkb_list)
    
    return(pg_geom)

def __pg_to_gdf(pg_geom,crs):

    wkb_array = pg.to_wkb(pg_geom)
    
    gdf = gpd.GeoDataFrame({'geometry' : [loads(g) for g in wkb_array]},
                           crs=crs,
                           geometry='geometry')

    return(shapely_list)

def __roundGeometry_pg(pg_geom,precision):
    
    round_pg = lambda g : round(g,precision)
    rounded_pg = pg.apply(pg_geom,round_pg)
    
    return(rounded_pg)


def __getDiffVal(validation,analysisExtents):
    
    try:
        diffVal = gpd.overlay(analysisExtents,validation,how='difference')
    except TopologicalError:
        return(None)
    diffVal = diffVal.reset_index(drop=True)
    diffVal = diffVal.loc[~diffVal.is_empty,:]
    diffVal = diffVal.reset_index(drop=True)

    return(diffVal)

def __getVectorDriver(filename):

    driverDictionary = {'.gpkg' : 'GPKG','.geojson' : 'GeoJSON','.shp' : 'ESRI Shapefile'}
    driver = driverDictionary[splitext(filename)[1]]

    return(driver)

def __simplify_gdf(gdf,simplify_tolerance):

    geomSimplified = gdf.simplify(simplify_tolerance)
    geomSimplified = geomSimplified.rename('geometry')
    gdf = gpd.GeoDataFrame({'geometry':geomSimplified},crs=gdf.crs,geometry='geometry')

    return(gdf)

def __fixGeometry_gdf(gdf):

    invalidGeometriesBoolean = ~gdf.geometry.is_valid
    if invalidGeometriesBoolean.any():
        geomBuffered = gdf.geometry.buffer(0)
        geomBuffered = geomBuffered.rename('geometry')
        return(gpd.GeoDataFrame({'geometry' : geomBuffered}, crs=gdf.crs,geometry='geometry'))
    else:
        return(gdf)

def __roundGeometry_gdf(gdf, precision=4):

    def around(coords, precision):
        result = []
        try:
            return round(coords, precision)
        except TypeError:
            for coord in coords:
                result.append(around(coord, precision))
        return(result)

    newGeoms = [None] * len(gdf)
    for i,geom in enumerate(gdf.geometry):
        geojson = mapping(geom)
        geojson['coordinates'] = around(geojson['coordinates'],precision)
        newGeoms[i] = shape(geojson)

    gdf = gpd.GeoDataFrame({'geometry' : newGeoms},crs=gdf.crs,geometry='geometry')

    return (gdf)

def __fishnet_gdf(gdf, threshold,verbose):

    bounds = gdf.geometry.total_bounds
    
    xmin = int(bounds[0] // threshold)
    xmax = int(bounds[2] // threshold)
    ymin = int(bounds[1] // threshold)
    ymax = int(bounds[3] // threshold)
    ncols = int(xmax - xmin + 1)
    nrows = int(ymax - ymin + 1)
    
    result = []
    for i in tqdm(range(xmin, xmax+1),disable=(not verbose)):
        for j in range(ymin, ymax+1):
            b = box(i*threshold, j*threshold, (i+1)*threshold, (j+1)*threshold)
            g = gdf.geometry.intersection(b)
            nonEmptyRows = ~g.is_empty
            nonEmptyGeometries = g.loc[nonEmptyRows]
            for ii in nonEmptyGeometries:
              result.append(ii)

    gdf = gpd.GeoDataFrame({ 'geometry' : result }, crs = gdf.crs , geometry = 'geometry')

    return(gdf)

def __fishnet_pg(geom, threshold,verbose):

    bounds = pg.get_coordinates(pg.boundary(geom))
    
    xmin = int(bounds[0,:] // threshold)
    xmax = int(bounds[2,:] // threshold)
    ymin = int(bounds[1,:] // threshold)
    ymax = int(bounds[3,:] // threshold)
    ncols = int(xmax - xmin + 1)
    nrows = int(ymax - ymin + 1)
    
    result = np.array([])
    for i in tqdm(range(xmin, xmax+1),disable=(not verbose)):
        for j in range(ymin, ymax+1):
            b = pg.box(i*threshold, j*threshold, (i+1)*threshold, (j+1)*threshold)
            g = pg.intersection(geom,b)
            result = np.append(result,g)

    return(result)


def __preprocess_forecasts(crossSections,flows,analysisExtents,crossSections_layerName=None,verbose=True):
    
    if isinstance(crossSections,str): crossSections = gpd.read_file(crossSections,layer=crossSections_layerName,mask=analysisExtents)
    if isinstance(flows,str): flows = gpd.read_file(flows,mask=analysisExtents)

    intersections = gpd.overlay(crossSections,flows,how='intersection')

    orig_flows = [('10yr','E_Q_10PCT'),('100yr','E_Q_01PCT'),('500yr','E_Q_0_2PCT')]
    dischargeMultiplier = 0.3048 ** 3

    forecasts = [None] * len(new_flows)
    for i,(nf,of) in enumerate(flows,orig_flows):
        forecast = intersections[['feature_id',of]]
        forecast = forecast.rename(columns={of : 'discharge'})
        forecast = forecast.astype({'feature_id' : int , 'discharge' : float})

        forecast = forecast.groupby('feature_id').median()
        forecast = forecast.reset_index(level=0)

        forecast['discharge'] = forecast['discharge'] * dischargeMultiplier
        #forecast.to_csv("forecast_{}_{}.csv".format(hucCode,nf),index=False)
        forecasts[i] = forecast

    return(forecasts,crossSections)

def __preprocess_extents(projection,analysisExtents,analysisExtents_layername=None,exclusionMask=None,exclusionMask_layername=None,split_threshold=None,verbose=True):

    __vprint('Analysis Extents ...',verbose)

    # load files
    __vprint('  Loading',verbose)
    if isinstance(analysisExtents,str): analysisExtents = gpd.read_file(analysisExtents,analysisExtents_layername)
    if (exclusionMask is not None) & isinstance(exclusionMask,str): exclusionMask = gpd.read_file(exclusionMask,mask=analysisExtents,layer=exclusionMask_layername) 

    # project
    __vprint('  Projecting',verbose)
    analysisExtents = analysisExtents.to_crs(projection)
    if exclusionMask is not None: exclusionMask = exclusionMask.to_crs(projection)

    # explode
    __vprint('  Exploding',verbose)
    analysisExtents = analysisExtents.explode().reset_index(level=1,drop=True)
    if exclusionMask is not None: exclusionMask = exclusionMask.explode().reset_index(level=1,drop=True)

    # fix geometries
    __vprint('  Buffering',verbose)
    analysisExtents = __fixGeometry_gdf(analysisExtents)
    if exclusionMask is not None: exclusionMask = __fixGeometry_gdf(exclusionMask)

    # split
    __vprint('  Splitting',verbose)
    if split_threshold is not None: analysisExtents = __fishnet_gdf(analysisExtents,split_threshold,verbose)

    # remove exclusion mask from analysisExtents and clip predicted and validation 
    if exclusionMask is not None:
        __vprint('  Removing exclusion mask',verbose)
        exclusionMask = gpd.overlay(analysisExtents,exclusionMask,how='intersection')
        analysisExtents = gpd.overlay(analysisExtents,exclusionMask,how='difference')
        analysisExtents = analysisExtents.explode().reset_index(level=1,drop=True)

    return(analysisExtents)


def __preprocess_validation(projection,validation,analysisExtents,test_case_level,split_threshold=None,simplify_tolerance=None,geometry_precision=2,verbose=True):

    __vprint('Validation: test case level ... {}'.format(test_case_level),verbose)
    
    # load files
    __vprint('  Loading',verbose)
    if isinstance(validation,str): validation = gpd.read_file(validation,mask=analysisExtents)

    # project
    __vprint('  Projecting',verbose)
    validation = validation.to_crs(projection)
     
    # Convert to Pygeos
    validation = __gdf_to_pg(validation)
    analysisExtents = __gdf_to_pg(analysisExtents)
   
     # round precisons
    __vprint('  Rounding',verbose)
    #validation = __roundGeometry_gdf(validation,geometry_precision)
    validation = __roundGeometry_pg(validation,geometry_precision)

    # simplify
    __vprint('  Simplifying',verbose)
    if simplify_tolerance is not None: validation = pg.simplify(validation,simplify_tolerance)

    # explode
    #__vprint('  Exploding',verbose)
    #validation = validation.explode().reset_index(level=1,drop=True)
    
    # fix geometries
    __vprint('  Making valid',verbose)
    validation = pg.make_valid(validation)

    # split
    __vprint('  Splitting',verbose)
    if split_threshold is not None: validation = __fishnet_pg(validation,split_threshold,verbose)
    
    # remove non-polygons
    #validation = validation.loc[validation.geometry.geom_type == 'Polygon',:]

    # Remove empties 
    #__vprint('  Removing empties',verbose)
    #validation = validation.loc[~validation.is_empty]
    #validation = validation.reset_index(drop=True)

    # round precisons
    __vprint('  Rounding',verbose)
    #validation = __roundGeometry_gdf(validation,geometry_precision)
    validation = __roundGeometry_pg(validation,geometry_precision)

    # fix geometries
    __vprint('  Making valid geometry',verbose)
    #validation = __fixGeometry_gdf(validation)
    validation = pg.make_valid(validation)

    # tree building
    valTree = pg.STRtree(validation)

    # clip validation to analysisExtents
    __vprint('  Clipping validation',verbose)

    #if not validation.geometry.intersects(analysisExtents).all():
    #    try:
    #        validation = gpd.overlay(validation,analysisExtents,how='intersection')
    #    except TopologicalError:
    #        return(None)

    #    validation = validation.explode().reset_index(drop=True)
    query = valTree.query_bulk(analysisExtents,predicate='intersects')
    validation = pg.intersection(analysisExtents[query[0,:]],validation[query[1,:]])

    # convert to gdf
    validation = __pg_to_gdf(validation,projection)

    return(validation)

def preprocess_test_case(validation,analysisExtents,crossSections,flows,crossSections_layerName,projection,test_case_directory,test_case_name,test_case_levels,exclusionMask=None,split_threshold=None,simplify_tolerance=None,geometry_precision=2,verbose=True):

    __vprint('Preprocess test case ... {}'.format(test_case_name),verbose)

    # forecasts
    forecasts,crossSections = __preprocess_forecasts(cs,flows,crossSections_layerName=crossSections_layerName,verbose=verbose)
    
    # analysis extents
    analysisExtents = __preprocess_extents(projection,analysisExtents,exclusionMask=exclusionMask,split_threshold=split_threshold,verbose=verbose)

    if isinstance(validation, str) | isinstance(validation,gpd.GeoDataFrame):
        validation = __preprocess_validation(projection,validation,analysisExtents,test_case_levels,split_threshold=split_threshold,simplify_tolerance=simplify_tolerance,geometry_precision=geometry_precision,verbose=verbose)
    else:
        validation = [__preprocess_validation(projection,val,analysisExtents,tc_level,split_threshold=split_threshold,simplify_tolerance=simplify_tolerance,geometry_precision=geometry_precision,verbose=verbose) for tc_level,val in zip(test_case_levels,validation)]
    
    __vprint('Validation difference ...',verbose)

    ## diff validation
    if isinstance(validation,list):
        diffVal = [__getDiffVal(val,analysisExtents) for val in validation]
    else:
        diffVal = __getDiffVal(validation,analysisExtents)

    # build file paths and names
    # validation
    validation_directory = join(test_case_directory,test_case_name,'validation')
    makedirs(validation_directory, exist_ok=True)

    analysisExtents_filename = join(validation_directory,'analysisExtents.gpkg')
    crossSections_filename = join(validation_directory,'crossSections.gpkg')
    forecast_filename_template = join(validation_directory,'forecast_{}.gpkg')
    validation_filename_template = join(validation_directory,'validation_{}.gpkg')
    diffVal_filename_template = join(validation_directory,'diffVal_{}.gpkg')
    
	# write to file
    __vprint('Writing to files ...',verbose)
    analysisExtents.to_file(analysisExtents_filename,driver=__getVectorDriver(analysisExtents_filename)) 
    crossSections.to_file(crossSections_filename,driver=__getVectorDriver(crossSections_filename))
   
    for f in forecasts:
        forecast_filename = forecast_filename_template.format(l)
        f.to_file(forecast_filename,driver=__getVectorDriver(forecast_filename)) 

    if isinstance(validation,list):
        for l,dv in zip(test_case_levels,diffVal):
            if dv is not None:
                diffVal_filename = diffVal_filename_template.format(l)
                dv.to_file(diffVal_filename,driver=__getVectorDriver(diffVal_filename))
            else:
                print("Not writing {test_case_name}, {l} due to Topology error")
        for l,v in zip(test_case_levels,validation):
            if v is not None:
                validation_filename = validation_filename_template.format(l)
                v.to_file(validation_filename,driver=__getVectorDriver(validationi_filename))
    else:
        diffVal_filename = diffVal_filename_template.format(l)
        diffVal.to_file(diffVal_filename,driver=__getVectorDriver(diffVal_filename))
        
        validation_filename = validation_filename_template.format(l)
        validation.to_file(validation_filename,driver=__getVectorDriver(validation_filename))

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Preprocess and add test case to library')
    parser.add_argument('-v','--validation', help='DEM to use within project path', required=True, nargs='*')
    parser.add_argument('-a','--analysis-extents',help='Pixel based watersheds raster to use within project path',required=True)
    parser.add_argument('-c','--cross_sections',help='Pixel based watersheds raster to use within project path',required=True)
    parser.add_argument('-w','--flows',help='Pixel based watersheds raster to use within project path',required=True)
    parser.add_argument('-y','--cross_section_layer',help='Pixel based watersheds raster to use within project path',required=True)
    parser.add_argument('-p','--projection',help='Output REM raster',required=True)
    parser.add_argument('-d','--test-case-directory',help='Output REM raster',required=True)
    parser.add_argument('-n','--test-case-name',help='Output REM raster',required=True)
    parser.add_argument('-f','--test-case-levels',help='Output REM raster',required=True, nargs='*')
    parser.add_argument('-e','--exclusion-mask',help='Output REM raster',required=False,default=None)
    parser.add_argument('-t','--split-threshold',help='Output REM raster',required=False,default=None,type=int)
    parser.add_argument('-l','--simplify-tolerance',help='Output REM raster',required=False,default=None,type=float)
    parser.add_argument('-r','--geometry-precision',help='Output REM raster',required=False,default=2,type=int)
    parser.add_argument('-q','--quiet',help='Output REM raster',required=False,action='store_false',default=True)
    
    # extract to dictionary
    args = vars(parser.parse_args())

    # rename variable inputs
    validation = args['validation']
    analysisExtents = args['analysis_extents']
    crossSections = args['cross_sections']
    flows = args['flows']
    crossSections_layerName = args['cross_section_layer']
    projection = args['projection']
    test_case_directory = args['test_case_directory']
    test_case_name = args['test_case_name']
    test_case_levels = args['test_case_levels']
    exclusionMask = args['exclusion_mask']
    split_threshold = args['split_threshold']
    simplify_tolerance = args['simplify_tolerance']
    geometry_precision = args['geometry_precision']
    verbose = args['quiet']
    
    preprocess_test_case(validation,analysisExtents,crossSections,flows,crossSections_layerName,projection,test_case_directory,test_case_name,test_case_levels,exclusionMask=exclusionMask,split_threshold=split_threshold,simplify_tolerance=simplify_tolerance,geometry_precision=geometry_precision,verbose=verbose)
