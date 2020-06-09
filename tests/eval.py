#!/usr/bin/env python3

import geopandas as gpd
import pandas as pd
import argparse
from os import makedirs
from os.path import join,splitext
from tqdm import tqdm
from shapely.geometry import mapping,shape,boxi
from shapely.errors import TopologicalError:

def __vprint(message,verbose=False):
    if verbose:
        print(message)

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

def __preprocess_cross_sections(crossSections,flows,test_case_level,crossSections_layerName=None,verbose=True):
    
    if isinstance(crossSections,str): crossSections = gpd.read_file(crossSections,crossSections_layerName)
    

def __preprocess_extents(projection,analysisExtents,exclusionMask=None,split_threshold=None,verbose=True):

    __vprint('Analysis Extents ...',verbose)

    # load files
    __vprint('  Loading',verbose)
    if isinstance(analysisExtents,str): analysisExtents = gpd.read_file(analysisExtents)
    if (exclusionMask is not None) & isinstance(exclusionMask,str): exclusionMask = gpd.read_file(exclusionMask,mask=analysisExtents) 

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
    
    # round precisons
    __vprint('  Rounding',verbose)
    validation = __roundGeometry_gdf(validation,geometry_precision)

    # simplify
    __vprint('  Simplifying',verbose)
    if simplify_tolerance is not None: validation = __simplify_gdf(validation,simplify_tolerance)

    # explode
    __vprint('  Exploding',verbose)
    validation = validation.explode().reset_index(level=1,drop=True)
    
    # fix geometries
    __vprint('  Buffering',verbose)
    validation = __fixGeometry_gdf(validation)

    # split
    __vprint('  Splitting',verbose)
    if split_threshold is not None: validation = __fishnet_gdf(validation,split_threshold,verbose)
    
    # remove non-polygons
    validation = validation.loc[validation.geometry.geom_type == 'Polygon',:]

    # Remove empties 
    __vprint('  Removing empties',verbose)
    validation = validation.loc[~validation.is_empty]
    validation = validation.reset_index(drop=True)
    
    # explode
    __vprint('  Exploding',verbose)
    validation = validation.explode().reset_index(level=1,drop=True)
    
    # round precisons
    __vprint('  Rounding',verbose)
    validation = __roundGeometry_gdf(validation,geometry_precision)
    
    # fix geometries
    __vprint('  Buffering',verbose)
    validation = __fixGeometry_gdf(validation)

    # clip validation to analysisExtents
    __vprint('  Clipping validation',verbose)

    if not validation.geometry.intersects(analysisExtents).all():
        try:
            validation = gpd.overlay(validation,analysisExtents,how='intersection')
        except TopologicalError:
            return(None)

        validation = validation.explode().reset_index(drop=True)

    return(validation)

def preprocess_test_case(validation,analysisExtents,projection,test_case_directory,test_case_name,test_case_levels,exclusionMask=None,split_threshold=None,simplify_tolerance=None,geometry_precision=2,verbose=True):

    __vprint('Preprocess test case ... {}'.format(test_case_name),verbose)

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
    validation_filename_template = join(validation_directory,'validation_{}.gpkg')
    diffVal_filename_template = join(validation_directory,'diffVal_{}.gpkg')

	# write to file
    __vprint('Writing to files ...',verbose)
    analysisExtents.to_file(analysisExtents_filename,driver=__getVectorDriver(analysisExtents_filename)) 
    
    if isinstance(validation,list):
        for l,dv in zip(test_case_levels,diffVal):
            diffVal_filename = diffVal_filename_template.format(l)
            dv.to_file(diffVal_filename,driver=__getVectorDriver(diffVal_filename))
        for l,v in zip(test_case_levels,validation):
            validation_filename = validation_filename_template.format(l)
            v.to_file(validation_filename,driver=__getVectorDriver(validation_filename))
    else:
        diffVal_filename = diffVal_filename_template.format(l)
        diffVal.to_file(diffVal_filename,driver=__getVectorDriver(diffVal_filename))
        
        validation_filename = validation_filename_template.format(l)
        validation.to_file(validation_filename,driver=__getVectorDriver(validation_filename))

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Preprocess and add test case to library')
    parser.add_argument('-v','--validation', help='DEM to use within project path', required=True, nargs='*')
    parser.add_argument('-a','--analysis-extents',help='Pixel based watersheds raster to use within project path',required=True)
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
    projection = args['projection']
    test_case_directory = args['test_case_directory']
    test_case_name = args['test_case_name']
    test_case_levels = args['test_case_levels']
    exclusionMask = args['exclusion_mask']
    split_threshold = args['split_threshold']
    simplify_tolerance = args['simplify_tolerance']
    geometry_precision = args['geometry_precision']
    verbose = args['quiet']
    
    preprocess_test_case(validation,analysisExtents,projection,test_case_directory,test_case_name,test_case_levels,exclusionMask=exclusionMask,split_threshold=split_threshold,simplify_tolerance=simplify_tolerance,geometry_precision=geometry_precision,verbose=verbose)
