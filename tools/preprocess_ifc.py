#!/usr/bin/env python

from pathlib import Path
import re
import pandas as pd
import geopandas as gpd

PREP_PROJECTION = 'PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'



flow_file = Path('/Path/to/flow/file')
geodatabase = 
flow_file = 
project_file =
workspace = 



def get_hec_ras_flows(flow_file, units):
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

def get_model_units(project_file):
    file = Path(project_file)
    with open(file) as f:
        model_unit_system = f.read().splitlines()[3]
    return model_unit_system
            
#Get spatial data
def get_xs(geodatabase):
    geodatabase = Path(geodatabase)
    gdb_gpd = gpd.read_file(geodatabase, layer = 'XSCutlines')
    gdb_gpd = gdb_gpd.filter(items= ['RiverCode','ReachCode','ProfileM', 'geometry'])
    return gdb_gpd

def assign_xs_flows(flow_file, geodatabase, workspace):
    #Get flows from HEC-RAS model
    flows = get_hec_ras_flows(flow_file)
    #Get flow units
    units = {'SI Units':'CMS', 'English Units':'CFS'}
    model_units = get_model_units(project_file)
    flow_units = units.get(model_units)
    #Get spatial XS data
    xs = get_xs(geodatabase)

        
    #Join flows with XS
    xs['ProfileM'] = xs['ProfileM'].astype(str)
    joined = xs.merge(flows, left_on = 'ProfileM', right_on = 'station')
    joined.drop(columns = ['river','reach','station'], inplace = True)
    joined['flow_units'] = flow_units
    
    #Reproject to FIM projection and export data to shapefile
    joined.to_crs(PREP_PROJECTION)
    
    stem = joined.RiverCode.drop_duplicates().item
    joined.to_file(workspace / f'{stem}_xs.shp')
