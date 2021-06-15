#!/usr/bin/env python

from pathlib import Path
import re
import pandas as pd
import geopandas as gpd

PREP_PROJECTION = 'PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'

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
#Preprocess IFC data
#Top level directory
directory = Path('Path to hydraulic models')
subdirs = [subdir for subdir in directory.iterdir() if subdir.is_dir()]
#Define workspace and create if it doesn't exist
workspace = Path('Path to workspace')
workspace.mkdir(exist_ok = True, parents = True)
#Define model source (usace or ifc)
source = 'usace'
for sub in subdirs:
    #Find the project file in the 'Simulations' directory
    [project_file] = [file for file in sub.rglob('*.prj') if file.parent.name == 'Simulations']
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

    #Find the geodatabase in subdir
    [geodatabase] = [gdb for gdb in sub.rglob('*.gdb')]
    assign_xs_flows(source, flow_file, project_file, geodatabase, workspace)


#Single TIME PREPROCESSING USACE DATA
#CONVERT ALL MDB to SHAPEFILE (CANT SEEM TO DO THIS OPEN SOURCE)
#CAN ONLY DO IN ARCMAP as ARCPRO DOES NOT RECOGNIZE MDB!!!!!!!!!
#Get all mdb
import os
import arcpy
PATH = 'Path to HUC 8 datasets
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
        
        
##############################################################################
#IFC Preprocessing
##############################################################################
from pathlib import Path
import geopandas as gpd

HUC_DIRECTORY = 
WORKSPACE = 
WORKSPACE.mkdir(parents = True, exist_ok=True)


#Check if any mdb files are present, if so convert to gdb using arcmap
mdb_files = list(HUC_DIRECTORY.rglob('*.mdb'))
if len(mdb_files) >= 1:
    print("CONVERT to GDB")
#Otherwise get all geb files
gdb_files = list(HUC_DIRECTORY.rglob('*.gdb'))
all_profile = gpd.GeoDataFrame()
for geodatabase in gdb_files:
    #Get all profile information and append
    profile = gpd.read_file(geodatabase, layer = 'ProfileDefinition',dtypes = {'Profile':str,'Alias':str})
    profile['source'] = str(geodatabase)
    all_profile = all_profile.append(profile, ignore_index = True)


##############################################################################
#Depth Grids Preprocessing
##############################################################################     
#Verify that all profiles are the same, they should be do this in FOSS
##############################################################################
#Depth Grids Preprocessing (ARCPRO REQUIRED)
##############################################################################     
#Verify that all profiles are the same, they should be do this in FOSS
from pathlib import Path
import arcpy
import pandas as pd
HUC_DIRECTORY = Path(r'Path/to/HUC/Directory')
WORKSPACE = Path('Path/to/Workspace')
WORKSPACE.mkdir(parents = True, exist_ok=True)
GRIDS_DIR = WORKSPACE/'grids'
GRIDS_DIR.mkdir(exist_ok = True)

gdb_files = list(HUC_DIRECTORY.rglob('*.gdb'))
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

REF_RASTER = Path(r'Path/to/raster')
CELL_SIZE = 10
for profile, alias in events.items():        
    depth_grid_name = f'd{profile.lower()}'
    #find all grids with this name
    grids =list(HUC_DIRECTORY.rglob(f'*{depth_grid_name}*'))
    grids = [str(arcgrid) for arcgrid in grids if arcgrid.is_dir()]
    arcpy.management.MosaicToNewRaster(input_rasters=grids, output_location=str(GRIDS_DIR), raster_dataset_name_with_extension=f'{alias.lower()}.tif', coordinate_system_for_the_raster=str(REF_RASTER), pixel_type='32_BIT_FLOAT', cellsize=CELL_SIZE, number_of_bands=1, mosaic_method='MAXIMUM')
    
    
    
###############################################################################
#Models Preprocessing (Copy model files)
###############################################################################
from pathlib import Path
import shutil
#Get a unique file associated with model (f01)
flow_files = list(HUC_DIRECTORY.rglob('*.f01'))
#Get name of model from flow file
for file in flow_files:
    reach = file.stem
    #Discover and then copy all model files
    model_dir   = file.parent
    model_files = list(model_dir.glob(f'{reach}.*'))
    #Make destination directory
    dest_model_dir = WORKSPACE/'models'/reach
    dest_model_dir.mkdir(parents = True, exist_ok = True)
    for model_file in model_files:
        shutil.copy(str(model_file), str(dest_model_dir/model_file.name))

###############################################################################

#Preprocess IFC data (USE_THIS)
WORKSPACE = Path('Path/to/preprocessing_directory')
project_files = list((WORKSPACE/'models').rglob('*.prj'))
for project_file in project_files:
    #Find the project file in the 'Simulations' directory
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
        


