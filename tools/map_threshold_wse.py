# -*- coding: utf-8 -*-
"""
Created on Tue May 25 12:30:34 2021

@author: Trevor.Grout
"""

#Step 1: Get Threshold Water Surface Elevation (action stage + gage datum)
#Get WSE of each stage: 
API_BASE_URL=
metadata_url = f'{API_BASE_URL}/metadata'
threshold_url = f'{API_BASE_URL}/nws_threshold'

#Get all nws_lid sites with datums
conus_list, conus_dataframe = get_metadata(metadata_url, select_by = 'nws_lid', selector = ['all'], must_include = 'nws_data.rfc_forecast_point', upstream_trace_distance = 10, downstream_trace_distance = 10 )
dictionary = {}
for metadata in conus_list:
    #Only use USGS datums for now
    nws, usgs = get_datum(metadata)    
    #If USGS datum is not supplied, skip
    if not usgs.get('datum'):
        continue
    #If USGS site code not supplied, skip
    if not usgs.get('usgs_site_code'):
        continue    
    datum = usgs.get('datum')

    if usgs.get('vcs') == 'NGVD29':                
        continue
        # #Convert NGVD to NAVD if needed
        # adj_ft = ngvd_to_navd_ft(datum_info = usgs, region = 'contiguous')
        # datum = datum + adj_ft
    
    #Get stages for site
    stages,flows = get_thresholds(threshold_url, select_by='nws_lid', selector=usgs.get('nws_lid'), threshold = 'all')
    threshold_categories = ['action','minor','moderate','major']
    #Check that valid threshold stages are supplied
    if not any([stages[threshold] for threshold in threshold_categories]):
        #Skipping because no threshold stages available
        continue
    #Convert stages to elevations
    threshold_elevation = {threshold: stages[threshold] + datum for threshold in threshold_categories}
    #Write site and thresholds/elevations to dictionary
    dictionary[usgs['usgs_site_code']] = threshold_elevation
###############################################################################
#Step 2: Get HAND stage (action water surface elevation - HAND datum) --> use usgs_elev_table.csv to get HAND datum/HydroID
#Path to FIM output
fim_outputs = Path('/Path/to/fim/output')
subdirs = [i for i in fim_outputs.iterdir() if i.is_dir()]
#Loop through each folder
for run in subdirs:
    usgs_elev_table = run / 'usgs_elev_table.csv'
    hydro_table = run / 'hydroTable.csv'

    if usgs_elev_table.exists and hydro_table.exists:
        #Read Tables
        usgs_elev_df = pd.read_csv(usgs_elev_table, dtype = {'location_id':str, 'HydroID':int})        
        hydro_table_df = pd.read_csv(hydro_table, dtype = {'HydroID':int,'feature_id':int,'HUC':str})
        
        #Dictionary of HAND datums and HydroIDs
        site_hand_datums = usgs_elev_df.set_index('location_id')['dem_adj_elevation'].to_dict()
        site_hydroid = usgs_elev_df.set_index('location_id')['HydroID'].to_dict()
        
        #for each location 
        for location in site_hand_datums:
            #get datum and hydroid
            hand_datum = site_hand_datums[location]
            hydroid = site_hydroid[location]
            #query appropriate rating curve
            rating_curve = hydro_table_df.query(f'HydroID == {hydroid}').copy()
            rating_curve['elevation_navd88_ft'] = rating_curve['stage']*3.28084 + hand_datum
            
            
            
            ####
            #!!!TEST THIS SECTION!!!
            #Step 3: Get HAND flow (Use rating curve to get flow corresponding to HAND stage) --> Use hydroTable.csv
            #Create DataFrame of thresholds/elevations for site
            interpolated_flow_cms_df = pd.DataFrame(dictionary[location].items(), columns = ['Threshold','Elevation'])
            #Interpolate HAND flow based on elevation
            interpolated_flow_cms_df['flow_cms'] = np.interp(interpolated_flow_cms_df['Elevation'], rating_curve['elevation_navd88_ft'], rating_curve['discharge'], left = np.nan, right = np.nan)


            #Step 4: Apply flow to all NWM segments (10 mi upstream/downstream)
            #Get mainstems segments (focus on mainstems)
            #Do intersection with NWM segments returned from original site and mainstems
            #Write flow file
            
