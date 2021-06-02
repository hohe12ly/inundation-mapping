#!/usr/bin/env python3

from pathlib import Path
from tools_shared_functions import mainstem_nwm_segs, get_metadata, get_thresholds, get_datum, ngvd_to_navd_ft, get_nwm_segs, flow_data
from dotenv import load_dotenv
import os
import pandas as pd
import geopandas as gpd
import numpy as np


load_dotenv()
#import variables from .env file
API_BASE_URL = os.getenv("API_BASE_URL")
EVALUATED_SITES_CSV = os.getenv("EVALUATED_SITES_CSV")


metadata_url = f'{API_BASE_URL}/metadata'
threshold_url = f'{API_BASE_URL}/nws_threshold'

#Get all possible mainstem segments
print('Getting list of mainstem segments')
#Import list of evaluated sites
list_of_sites = pd.read_csv(EVALUATED_SITES_CSV)['Total_List'].to_list()
#The entire routine to get mainstems is hardcoded in this function.
ms_segs = mainstem_nwm_segs(metadata_url, list_of_sites)


def get_thresh_elevs(sites):
    if not type(sites) == 'list':
        sites = [sites]
        
    #Get all nws_lid sites with datums
    metadata_list, metadata_dataframe = get_metadata(metadata_url, select_by = 'usgs_site_code', selector = sites, upstream_trace_distance = 10, downstream_trace_distance = 10 )
    thresh_elevs_m_dict = {}
    
    for metadata in metadata_list:
        #Only use USGS datums for now
        nws, usgs = get_datum(metadata)    
    
        #Skip certain sites
        if not (usgs.get('datum') or usgs.get('usgs_site_code')) :
            print('skipping')
            continue
        
        #Get datum
        datum = usgs.get('datum')
        if usgs.get('vcs') == 'NGVD29':                
            #Convert NGVD to NAVD if needed
            adj_ft = ngvd_to_navd_ft(datum_info = usgs, region = 'contiguous')
            datum = datum + adj_ft
        
        #Get stages for site
        stages,flows = get_thresholds(threshold_url, select_by='nws_lid', selector=usgs.get('nws_lid'), threshold = 'all')
        threshold_categories = ['action','minor','moderate','major']
        #Check that at least 1 threshold is valid per site
        if not any([stages[threshold] for threshold in threshold_categories]):
            #Skipping because no threshold stages available
            continue
    
        #Convert stages to elevations
        threshold_elevation_m = {threshold: round((stages[threshold] + datum)/3.28084,2) for threshold in threshold_categories if stages[threshold]}
        #Write site and thresholds/elevations to dictionary
        thresh_elevs_m_dict[usgs['usgs_site_code']] = threshold_elevation_m
    return thresh_elevs_m_dict, metadata_list
###############################################################################
#Step 2: Get HAND stage (action water surface elevation - HAND datum) --> use usgs_elev_table.csv to get HAND datum/HydroID
#Path to FIM output
fim_output_dir = Path('/Path/to/fim/output')
workspace = Path('Path/to/workspace')


fim_subdirs = [i for i in fim_output_dir.iterdir() if i.is_dir()]
flood_categories = ['action','minor','moderate','major','record']
#Loop through each folder
for subdir in fim_subdirs:    
    print(f'Working on {subdir}')
    usgs_elev_table = subdir / 'usgs_elev_table.csv'
    hydro_table = subdir / 'hydroTable.csv'
    #Verify tables exist
    if not (usgs_elev_table.exists() and hydro_table.exists()):
        continue
    
    #Get HUC unit
    huc = subdir.name
    
    #Read Tables
    usgs_elev_df = pd.read_csv(usgs_elev_table, dtype = {'location_id':str, 'HydroID':int}, index_col = 'location_id')        
    hydro_table_df = pd.read_csv(hydro_table, dtype = {'HydroID':int,'feature_id':int,'HUC':str})
    
    #Dictionary of HAND datums and HydroIDs
    site_hand_datums = usgs_elev_df['dem_adj_elevation'].to_dict()
    site_hydroid = usgs_elev_df['HydroID'].to_dict()
    #for each location 
    for location in site_hand_datums:
        #get datum and hydroid
        hand_datum_m = site_hand_datums[location]
        hydroid = site_hydroid[location]
        #query appropriate rating curve
        rating_curve = hydro_table_df.query(f'HydroID == {hydroid}').copy()
        rating_curve['elevation_navd88_m'] = rating_curve['stage'] + hand_datum_m
        
        #Step 3: Get HAND flow (Use rating curve to get flow corresponding to HAND stage) --> Use hydroTable.csv
        dictionary, metadata = get_thresh_elevs(location)
        lid = metadata['identifiers']['nws_lid'].lower()
        #Create DataFrame of thresholds/elevations for site
        interpolated_flow_cms_df = pd.DataFrame(dictionary[location].items(), columns = ['Threshold','Elevation_m'])
        #Interpolate HAND flow based on elevation
        interpolated_flow_cms_df['flow_cms'] = np.interp(interpolated_flow_cms_df['Elevation_m'], rating_curve['elevation_navd88_m'], rating_curve['discharge_cms'], left = np.nan, right = np.nan)
        #Create flows dictionary
        flows = interpolated_flow_cms_df.set_index('Threshold')['flow_cms'].to_dict()
                
        #For Location = 07016500 (UNNM7, huc = 07140103)
        #WRDS FLOWS = Action: 8,880.6 cfs, Minor: 11,192.7 cfs, Moderate: 22,131 cfs, Major: 32,573.2 cfs
        #Interpolated Flows (using enforced WSE level) = Action: 3,787.96 cfs, Minor: 6,349.38 cfs, Moderate: 35,299.98 cfs, Major: 80746.83 cfs
        # CSI: Action: 0.274, Minor: 0.317, Moderate: 0.544, Major:0.606
        # FAR: Action: 0.656, Minor: 0.582, Moderate: 0.268, Major: 0.187
        # TPR: Action: 0.581, Minor: 0.565, Moderate: 0.680, Major: 0.704
        
        #Step 4: Apply flow to all NWM segments (10 mi upstream/downstream)
        #Get mainstem segments of LID by intersecting LID segments with known mainstem segments.
        [metadata_dict] = metadata
        segments = get_nwm_segs(metadata_dict)        
        site_ms_segs = set(segments).intersection(ms_segs)
        segments = list(site_ms_segs)  
                
        #Write flow file
        #if no segments, write message and exit out
        if not segments:
            print(f'{lid} no segments')
            continue
        #For each flood category
        for category in flood_categories:
            #Get the flow
            flow = flows.get(category)
            #If there is a valid flow value, write a flow file.
            if flow:
                #round flow to nearest hundredth
                flow = round(flow,2)
                #Create the guts of the flow file.
                flow_info = flow_data(segments,flow, convert_to_cms = False)
                #Define destination path and create folders
                output_file = workspace / huc / lid / category / (f'ahps_{lid}_huc_{huc}_flows_{category}.csv') 
                output_file.parent.mkdir(parents = True, exist_ok = True)
                #Write flow file to file
                flow_info.to_csv(output_file, index = False)                                                                                

# if __name__ == '__main__':
#     #Parse arguments
#     parser = argparse.ArgumentParser(description = '')
#     parser.add_argument('-s', '--source_dir', help = 'Workspace where all source data is located.', required = True)
#     parser.add_argument('-d', '--destination',  help = 'Directory where outputs are to be stored', required = True)
#     args = vars(parser.parse_args())
    