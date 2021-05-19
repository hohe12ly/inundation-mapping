#!/usr/bin/env python

from pathlib import Path
import re
import pandas as pd
import geopandas as gpd

PREP_PROJECTION = 'PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'

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


#Get data from globus using GLOBUS-CLI
import subprocess
from pathlib import Path
import re
import pandas as pd
# For Reference: https://docs.globus.org/cli/reference/ls/

#Global variables
SOURCE_EP = 
DEST_EP = 
DEST_EP_WORKSPACE = '~/Test'

#Get all directories 2 layers deep (hucs/reaches)
bashCommand = f"globus ls {SOURCE_EP}:/ --recursive --recursive-depth-limit 1 --jmespath DATA[*].[name,type,size] --format unix"
process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
output, error = process.communicate()
hucs_reaches = output.splitlines()
hucs_reaches = [a.decode("utf-8") for a in hucs_reaches]
hucs_reaches = [re.split(r'\t+', a) for a in hucs_reaches]

#Filter out hucs_reaches we don't want
hr_df = pd.DataFrame(hucs_reaches, columns = ['path','type', 'size'])
hr_df = hr_df[hr_df.path.str.startswith(('1','0'))& hr_df.path.str.contains('/')]
hr_df = hr_df[hr_df.type != 'file']
#split datasets to IFC and USACE (need to remove phantom _IFC)
usace_df = hr_df[hr_df.path.str.contains('_USACE') & ~ hr_df.path.str.endswith('_IFC')].copy()
ifc_df = hr_df[~hr_df.path.str.contains('_USACE')].copy()

#######################################################################
#Example to transfer USACE data
#usace has small sizes, will download entire contents of these folders
usace_paths = usace_df.path.to_list()
for path in usace_paths:
    bashCommand = f"globus transfer {SOURCE_EP}:/{path} {DEST_EP}:{DEST_EP_WORKSPACE} --recursive"
    process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
    output, error = process.communicate()

#######################################################################
#Example to Transfer IFC data (just model files and hydraulics.gdb)
ifc_df['model_path'] = ifc_df['path'] + '/Hydraulics'
ifc_df['gdb_path'] = ifc_df['model_path'] + '/Hydraulics/Hydraulics.gdb'

model_dirs = ifc_df.model_path.to_list()
gdb_paths = ifc_df.gdb_path.to_list()

#Download models
for model_dir in model_dirs:
    #Get model files
    reach = model_dir.split('/')[-2]
    bashCommand = f"globus ls {SOURCE_EP}:/{model_dir} --filter ~{reach}.* --jmespath DATA[*].[name,type,size] --format unix"
    process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
    output, error = process.communicate()
    model_files = output.splitlines()
    model_files = [a.decode("utf-8") for a in model_files]
    model_files = [re.split(r'\t+', a) for a in model_files]
    model_files_df = pd.DataFrame(model_files, columns = ['name','type','size'])
    #Make sure we have files only
    model_files_df = model_files_df[model_files_df.type == 'file']
   
    #Download model files
    model_file_names = model_files_df.name.to_list()
    for model_file in model_file_names:
        bashCommand = f"globus transfer {SOURCE_EP}:/{model_dir}/{model_file} {dest_ep}:{DEST_EP_WORKSPACE}/{model_dir}/{model_file}"
        process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
        output, error = process.communicate()

#Download Hydraulics GDB
for path in gdb_paths:
    bashCommand = f"globus transfer {SOURCE_EP}:/{path} {dest_ep}:{DEST_EP_WORKSPACE/{path} --recursive"
    process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
    output, error = process.communicate()




#Abandoned this in favor of Globus CLI   (Above)  
###############################################################################    
#Using Globus-sdk python library (STILL NEED TO MODIFY)
import globus_sdk

CLIENT_ID = 'INSERT_CLIENT_ID'

client = globus_sdk.NativeAppAuthClient(CLIENT_ID)
client.oauth2_start_flow()

authorize_url = client.oauth2_get_authorize_url()
print('Please go to this URL and login:\n\n{0}'.format(authorize_url))

# this is to work on Python2 and Python3 -- you can just use raw_input() or
# input() for your specific version
get_input = getattr(__builtins__, 'raw_input', input)
auth_code = get_input(
    'Please enter the code you get after login here: ').strip()
token_response = client.oauth2_exchange_code_for_tokens(auth_code)

globus_auth_data = token_response.by_resource_server['auth.globus.org']
globus_transfer_data = token_response.by_resource_server['transfer.api.globus.org']

# most specifically, you want these tokens as strings
AUTH_TOKEN = globus_auth_data['access_token']
TRANSFER_TOKEN = globus_transfer_data['access_token']

# a GlobusAuthorizer is an auxiliary object we use to wrap the token. In
# more advanced scenarios, other types of GlobusAuthorizers give us
# expressive power
authorizer = globus_sdk.AccessTokenAuthorizer(TRANSFER_TOKEN)
tc = globus_sdk.TransferClient(authorizer=authorizer)

# high level interface; provides iterators for list responses
print("My Endpoints:")
for ep in tc.endpoint_search('INSERT_GLOBUS_DATASET'):
    print("[{}] {}".format(ep["id"], ep["display_name"]))
   
directory=[]
path = '/10280201'
for entry in tc.operation_ls(ep['id'], path = path):
    if entry['type'] == 'dir':
        subpath = f'{path}/{entry["name"]}'
        print (entry['name'], entry['type'])    