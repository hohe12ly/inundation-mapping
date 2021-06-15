# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 11:27:44 2021

@author: Trevor.Grout
"""
import shutil
from pathlib import Path
import re
import pandas as pd
import geopandas as gpd
import os
import fiona
import argparse
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import rasterio.mask
import numpy as np


HUC_DIRECTORY
WORKSPACE
PREP_PROJECTION = 'PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'



def preprocess_benchmark_static(benchmark_raster, reference_raster, out_raster_path = None):
    '''
    This function will preprocess a benchmark dataset for purposes of evaluating FIM output. A benchmark dataset will be transformed using properties (CRS, resolution) from an input reference dataset. The benchmark raster will also be converted to a boolean (True/False) raster with inundated areas (True or 1) and dry areas (False or 0). 

    Parameters
    ----------
    benchmark_raster : STRING
        Path to benchmark raster (e.g. BLE elevation or WSE elevations).
    reference_raster : STRING
        Path to reference raster (e.g. output from inundation.py).
    out_raster_path : STRING, optional
        Path to the output raster. The default is None.

    Returns
    -------
    boolean_benchmark : Numpy Array
        The preprocessed benchmark array.
    profile : STRING
        Raster profile information for the preprocessed benchmark array (required for writing to output dataset).

    '''
    #Open and read raster and benchmark rasters
    reference = rasterio.open(reference_raster)
    benchmark = rasterio.open(benchmark_raster)
    benchmark_arr = benchmark.read(1)    
    #Set no data value (keep it integer so we can write to int8 raster)
    nodata_value = -2147483648
    
    #Determine the new transform and dimensions of reprojected/resampled raster.
    new_transform, new_width, new_height = calculate_default_transform(benchmark.crs, reference.crs, benchmark.width, benchmark.height, *benchmark.bounds, resolution = reference.res)

    #Define an empty array that is same dimensions as output by the "calculate_default_transform" command. 
    benchmark_projected = np.empty((new_height,new_width), dtype=np.int32)

    #Reproject and resample the benchmark dataset.
    reproject(benchmark_arr, 
              destination = benchmark_projected,
              src_transform = benchmark.transform, 
              src_crs = benchmark.crs,
              src_nodata = benchmark.nodata,
              dst_transform = new_transform, 
              dst_crs = reference.crs,
              dst_nodata = nodata_value,
              dst_resolution = reference.res,
              resampling = Resampling.bilinear)

    #Convert entire depth grid to boolean (1 = Flood, 0 = No Flood)
    boolean_benchmark = np.where(benchmark_projected != nodata_value, 1, 0)

    #Update profile (data type, NODATA, transform, width/height).
    profile = reference.profile
    profile.update(transform = new_transform)
    profile.update(dtype = rasterio.int8)
    profile.update(nodata = 2) #Update NODATA to some integer so we can keep int8 datatype. There are no NODATA in the raster dataset.
    profile.update (width = new_width)
    profile.update(height = new_height)

    #Write out preprocessed benchmark array to raster if path is supplied
    if out_raster_path is not None:    
        with rasterio.Env():    
            #Write out reassigned values to raster dataset.
            with rasterio.open(out_raster_path, 'w', **profile) as dst:
                dst.write(boolean_benchmark.astype('int8'),1)   
    return boolean_benchmark.astype('int8'), profile

def write_ifc_flow_file(ifc_xs_layer, nwm_geodatabase):
    #Change these if needed
    nwm_stream_layer_name = 'RouteLink_FL'
    nwm_feature_id_field ='ID'
        
    #Search for the layer that has 'XS' in the ble geodatabase. There should be only one and this is the 1D XS used in modeling. Read the layer into a geopandas dataframe.
    xs_layer = gpd.read_file(ifc_xs_layer)
    huc = ifc_xs_layer.parent.parent.name
    units = xs_layer.flow_units.unique().item()
    if units == 'CMS':
        #No conversion units already CMS
        dischargeMultiplier = 1
    elif units == 'CFS':
        #Convert from CFS to CMS
        dischargeMultiplier = 0.3048 ** 3
    
    #Search for the layer that has nwm_stream_layer_name (default 'RouteLink_FL' from NWMv2.1) in the nwm geodatabase. Read in this layer using the bounding box option based on the extents of the BLE XS layer. According to geopandas documentation, "CRS mis-matches are resolved if given a GeoSeries or GeoDataFrame." The NWM layer is read in as a geopandas dataframe. 
    [nwm_layer_name] = [i for i in fiona.listlayers(nwm_geodatabase) if nwm_stream_layer_name in i]
    nwm_river_layer = gpd.read_file(nwm_geodatabase, bbox = xs_layer, layer = nwm_layer_name)
    
    #make sure xs_layer is in same projection as nwm_river_layer.
    xs_layer_proj = xs_layer.to_crs(nwm_river_layer.crs)
    
    #Perform an intersection of the BLE layers and the NWM layers, using the keep_geom_type set to False produces a point output.
    intersection = gpd.overlay(xs_layer_proj, nwm_river_layer, how = 'intersection', keep_geom_type = False)
    
    #Create the flow forecast files
    #define fields containing flow (typically these won't change for BLE)
    flow_fields = ['DSCH_50PCT','DSCH_20PCT','DSCH_10PCT','DSCH_4PCT','DSCH_2PCT','DSCH_1PCT','DSCH_05PCT','DSCH_02PCT']
    
    #define return period associated with flow_fields (in same order as flow_fields). These will also serve as subdirectory names.
    return_period = ['2yr','5yr','10yr','25yr','50yr','100yr','200yr','500yr']
        
    #Write individual flow csv files
    for i,flow in enumerate(flow_fields):
        #Write dataframe with just ID and single flow event
        forecast = intersection[[nwm_feature_id_field,flow]]
    
        #Rename field names and re-define datatypes
        forecast = forecast.rename(columns={nwm_feature_id_field :'feature_id',flow : 'discharge'})
        forecast = forecast.astype({'feature_id' : int , 'discharge' : float})
    
        #Calculate median flow per feature id
        forecast = forecast.groupby('feature_id').median()
        forecast = forecast.reset_index(level=0)
    
        #Convert CFS to CMS
        forecast['discharge'] = forecast['discharge'] * dischargeMultiplier
    
        #Set paths and write file
        output_dir = WORKSPACE/f'{huc}_ifc'/return_period[i]
        output_dir.mkdir(parents = True, exist_ok = True)
        forecast.to_csv(output_dir /f"ble_huc_{huc}_flows_{return_period[i]}.csv" ,index=False) 
    
    
def get_hec_ras_flows(flow_file):
    '''
    Retrieves flows from HEC-RAS flow file.

    Parameters
    ----------
    flow_file : str
        Path to flow file.

    Returns
    -------
    all_flows : pandas DataFrame
        Dataframe of flows and other information from HEC-RAS flow file.

    '''

    all_flows = pd.DataFrame()
    with open(flow_file) as f:
        #For each line in file
        for line in f:
            
            #Get number of profiles        
            if line.startswith('Number of Profiles='):
                # determine the number of profiles
                [profile_key, profiles] = line.split('=') 
                profiles = profiles.strip()
            
            #Get the profile names
            if line.startswith('Profile Names='):
                [profile_name_key, *profile_names] = re.split('=|,',line)
                profile_names = [name.strip() for name in profile_names]
                
                #flow lines only contain 10 flow values before wrapping to next line. This script
                #only retrieves the first line of flows. Only flows with the first 10 profile 
                #names are retrieved.
                if len(profile_names) > 10:
                    profile_names = profile_names[0:9]
                
            #Get the flow values
            if line.startswith('River Rch & RM='):
                #Get line with river/reach/xs information
                split_line = re.split('=|,', line)
                [trash, river, reach, station] = [part.strip() for part in split_line]
                #flow values are on next line, with individual flow values taking 8 characters.
                flows_line = next(f).strip('\n')            
                flow_slot = 8 #flow profiles every 8 characters
                flows = [flows_line[i:i+flow_slot].strip() for i in range(0,len(flows_line),flow_slot)]
                
                #Write flow values to pandas df along with river/reach/xs.
                flow_df = pd.DataFrame(data = flows).T
                flow_df.columns = profile_names
                flow_df[profile_key] = profiles
                flow_df['river'] = river
                flow_df['reach'] = reach
                flow_df['station'] = station
                
                #Append flow information
                all_flows = all_flows.append(flow_df, ignore_index = True)
    return all_flows

def get_model_info(project_file):
    file = Path(project_file)
    with open(file) as f:
        contents = f.read().splitlines()
        [plan_file_ext] = [i.split('=')[-1] for i in contents if i.startswith('Current Plan=')]
        model_unit_system = contents[3]
    return model_unit_system, plan_file_ext
 
def get_model_files(plan_file):
    file = Path(plan_file)
    with open(file) as f:
        contents = f.read().splitlines()
        [geom_file_ext] = [i.split('=')[-1] for i in contents if i.startswith('Geom File=')]
        [flow_file_ext] = [i.split('=')[-1] for i in contents if i.startswith('Flow File=')]
    return geom_file_ext, flow_file_ext
           
#Get spatial data
def get_xs(source, geodatabase):
    geodatabase = Path(geodatabase)
    if source == 'ifc':
        gdb_gpd = gpd.read_file(geodatabase, layer = 'XSCutlines')
        gdb_gpd = gdb_gpd.filter(items= ['RiverCode','ReachCode','ProfileM', 'geometry'])
    elif source == 'usace':
        gdb_gpd = gpd.read_file(geodatabase, layer = 'S_XS')
        gdb_gpd = gdb_gpd.filter(items= ['WTR_NM','STREAM_STN','geometry'])
    return gdb_gpd
def get_river(source, geodatabase, spatial_dir):
    geodatabase = Path(geodatabase)
    if source == 'ifc':
        gdb_gpd = gpd.read_file(geodatabase, layer = 'River2D')
        gdb_gpd = gdb_gpd.filter(items= ['RiverCode','ReachCode','geometry'])
        stem = gdb_gpd.RiverCode.unique().item()
    elif source == 'usace':
        gdb_gpd = gpd.read_file(geodatabase, layer = 'S_XS')
        gdb_gpd = gdb_gpd.filter(items= ['WTR_NM','STREAM_STN','geometry'])
    
    #Reproject to FIM projection and export data to shapefile
    gdb_gpd.to_crs(PREP_PROJECTION)
    #Write to file
    gdb_gpd.to_file(spatial_dir / f'{stem}_river.shp') 
def assign_xs_flows(source, flow_file, project_file, geodatabase, workspace):
    #Get flows from HEC-RAS model
    flows = get_hec_ras_flows(flow_file)
    #Get flow units
    units = {'SI Units':'CMS', 'English Units':'CFS'}
    model_units, trash = get_model_info(project_file)
    flow_units = units.get(model_units)
    #Get spatial XS data
    xs = get_xs(source, geodatabase)

    #Merge flows to XS
    if source == 'ifc':        
        #Join flows with XS
        xs['ProfileM'] = xs['ProfileM'].astype(str)
        joined = xs.merge(flows, left_on = 'ProfileM', right_on = 'station')
        joined.drop(columns = ['river','reach','station'], inplace = True)
        joined['flow_units'] = flow_units
        stem = joined.RiverCode.unique().item()
    elif source == 'usace':
        xs['STREAM_STN'] = xs['STREAM_STN'].astype(str)
        joined = xs.merge(flows, left_on = 'STREAM_STN', right_on = 'station')
        joined.drop(columns = ['river','reach','station'], inplace = True)
        joined['flow_units'] = flow_units
        stem = joined.WTR_NM.unique().item()
    
    joined.rename(columns = {'Number of Profiles':'Num_Profi'}, inplace = True)          
    #Reproject to FIM projection and export data to shapefile
    joined.to_crs(PREP_PROJECTION)
    #Write to file
    joined.to_file(workspace / f'{stem}_xs.shp')
###############################################################################
###############################################################################
#Step 1: Find all mdb and convert to gdb
###############################################################################
mdbs = list(HUC_DIRECTORY.rglob("*.mdb"))
if mdbs:
    print('Use ARCMAP TO CONVERT MDB TO GDB FIRST!')
    #USE ARCMAP TO CONVERT   
    message = '''
    #Paste this script in ARCMAP(NOT ARCPRO)
    import os
    import arcpy
    PATH = r'{0}'
    mdbs = [os.path.join(dp, f) for dp, dn, filenames in os.walk(PATH) for f in filenames if os.path.splitext(f)[1] == '.mdb']
    for mdb in mdbs:
        print(mdb)
        #Write geodatabase
        out_directory = os.path.dirname(mdb)
        gdb_name = os.path.splitext(os.path.basename(mdb))[0]
        arcpy.CreateFileGDB_management(out_folder_path=out_directory, out_name=gdb_name, out_version="CURRENT")
        #Copy mdb layers from mdb to gdb
        spatial_layers = ['S_Fld_Haz_Ar','S_Profil_Basln','S_XS']
        for layer in spatial_layers:
            in_layer = os.path.join(mdb,'FIRM_Spatial_Layers',layer)
            out_layer = os.path.join(out_directory, gdb_name+'.gdb',layer)
            arcpy.Copy_management(in_layer,out_layer)
    '''
    print(message.format(str(HUC_DIRECTORY)))

#Step 2: Copy all gdb files to preprocess directory
###############################################################################
gdb_files = list(HUC_DIRECTORY.rglob('*.gdb'))
for gdb in gdb_files:
    dest = WORKSPACE/'spatial'/gdb.parent.parent.parent.name
    dest.mkdir(parents=True, exist_ok = True)
    shutil.copytree(gdb, dest/gdb.name)

#Step 3: Copy models to preprocess directory
##############################################################################
#Get a unique file associated with model (f01)
flow_files = list(HUC_DIRECTORY.rglob('*.f01'))
#Copy model files
for file in flow_files:
    reach = file.stem
    #Discover and then copy all model files
    model_dir   = file.parent
    model_files = list(model_dir.glob(f'{reach}.*'))
    #Make destination director
    dest_model_dir = WORKSPACE/'models'/reach
    dest_model_dir.mkdir(parents = True, exist_ok = True)
    for model_file in model_files:
        shutil.copy(str(model_file), str(dest_model_dir/model_file.name))

#Step 4: Copy Grids and process (DO THIS IN ARCPRO)
###############################################################################
#Get profile names associated with depth grid names
import arcpy
from pathlib import Path
WORKSPACE = Path('Path/to/workspace')
GRIDS_DIR = WORKSPACE/'grids'
HUC_DIRECTORY=Path('Path/to/source/data')
gdb_files = list(WORKSPACE.rglob('*.gdb'))
all_profile = pd.DataFrame()
for geodatabase in gdb_files:
    geodatabase = geodatabase/'ProfileDefinition'
    data = arcpy.da.TableToNumPyArray(str(geodatabase), ('Profile','Alias'))
    profile = pd.DataFrame(data)
    profile['source'] = str(geodatabase)
    all_profile = all_profile.append(profile, ignore_index = True)
events = all_profile.groupby('Profile')['Alias'].agg('unique')
events.to_csv(GRIDS_DIR/'profiles.csv')
events = events.to_dict()
for profile, aliases in events.items():
    if len(aliases) == 1:
        events.update({profile:aliases.item()})
        print(f'{profile} has {aliases} value only')
#Assuming the aliases are consistent throughout a HUC, preprocess GRIDS
REF_RASTER = Path(r'D:\ifc\ref_ras.tif') #Or path to reference raster
CELL_SIZE = 10
CRS = arcpy.Describe (str(REF_RASTER)).spatialReference
for profile, alias in events.items():        
    depth_grid_name = f'd{profile.lower()}'
    #find all grids with this name
    grids =list(HUC_DIRECTORY.rglob(f'*{depth_grid_name}*'))
    grids = [str(arcgrid) for arcgrid in grids if arcgrid.is_dir()]
    arcpy.management.MosaicToNewRaster(input_rasters=grids, output_location=str(GRIDS_DIR), raster_dataset_name_with_extension=f'{alias.lower()}.tif', coordinate_system_for_the_raster=CRS, pixel_type='32_BIT_FLOAT', cellsize=CELL_SIZE, number_of_bands=1, mosaic_method='MAXIMUM')

#Step 5: Get XS (populated with flows) and rivers.
###############################################################################
#Preprocess IFC data (USE_THIS)

project_files = list((WORKSPACE/'models').rglob('*.prj'))
for project_file in project_files:
    project_file_dir = project_file.parent
    #from project file get the plan file extension
    trash, plan_file_ext = get_model_info(project_file)
    #Get the plan file specified in project file
    plan_file = list(project_file_dir.glob(f'*.{plan_file_ext}'))
    #If plan file exists then get the flow file
    if plan_file:
        [plan_file] = plan_file
        trash, flow_file_ext= get_model_files(plan_file)
        [flow_file] = list(project_file_dir.glob(f'*.{flow_file_ext}'))
    #If the plan file does not exist, find the plan file that is present in directory (assumes 1) and get the flow file specified in that plan file.
    else:
        print(f'using available plan file {project_file_dir}')
        [plan_file] = list(project_file_dir.glob(f'*.p0*'))
        trash, flow_file_ext= get_model_files(plan_file)
        [flow_file] = list(project_file_dir.glob(f'*.{flow_file_ext}'))
    
    #Define model source (usace or ifc)
    source = 'ifc'
    spatial_dir = WORKSPACE/'spatial'/project_file.parent.name 
    geodatabase = spatial_dir/'Hydraulics.gdb'
    #Create XS/River shapefiles
    assign_xs_flows(source, flow_file, project_file, geodatabase, spatial_dir)
    get_river(source, geodatabase,spatial_dir)

#Concatenate XS/River Shapefiles
xs_files = list(WORKSPACE.rglob('*_xs.shp'))
all_xs = gpd.GeoDataFrame()
for xs in xs_files:
    temp = gpd.read_file(xs)
    all_xs = all_xs.append(temp)

river_files = list(WORKSPACE.rglob('*_river.shp'))
all_river = gpd.GeoDataFrame()
for river in river_files:
    temp = gpd.read_file(river)
    all_river = all_river.append(temp)

all_xs.to_file(WORKSPACE/'spatial'/'xs.shp')
all_river.to_file(WORKSPACE/'spatial'/'river.shp')   

#6. Process grids and write flow files
##############################################################################
##############################################################################
#Write Flow Files
###############################################################################
ifc_xs_layer =  WORKSPACE/'spatial'/'xs.shp'
nwm_geodatabase = 
write_ifc_flow_file(ifc_xs_layer, nwm_geodatabase)

###############################################################################
#Preprocess Benchmark Grids
###############################################################################
benchmark_rasters = list(WORKSPACE.rglob('*.tif'))
#define fields containing flow (typically these won't change for BLE)
flow_fields = ['DSCH_50PCT','DSCH_20PCT','DSCH_10PCT','DSCH_4PCT','DSCH_2PCT','DSCH_1PCT','DSCH_05PCT','DSCH_02PCT']
#define return period associated with flow_fields (in same order as flow_fields). These will also serve as subdirectory names.
return_intervals = ['2yr','5yr','10yr','25yr','50yr','100yr','200yr','500yr']
flow_dict = dict(zip(flow_fields, return_intervals))
for raster in benchmark_rasters:
    event = raster.stem
    return_period = flow_dict[event]
    huc = raster.parent.parent.name
    output_path = WORKSPACE/f'{huc}_ifc'/return_period/f'ble_huc_{huc}_extent_{return_period}.tif'
    output_path.parent.mkdir(parents = True, exist_ok = True)    

    preprocess_benchmark_static(raster, REF_RASTER, out_raster_path = output_path)