#!/usr/bin/env python3

'''
Description:
    1) buffer the MS demDerived_reaches.shp
    2) locate all FR reaches within MS buffer
    3) export FR-->MS subsetted reaches
'''

import sys
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiPoint

fr_dem_flows_fileName     = sys.argv[1]
ms_flows_fileName         = sys.argv[2]
wbd_huc_filename          = sys.argv[3]
ms_dem_flows_fineName     = sys.argv[4]

print('Loading data ...')
fr_dem_flows = gpd.read_file(fr_dem_flows_fileName)
wbd_huc_boundary = gpd.read_file(wbd_huc_filename)
ms_flows = gpd.read_file(ms_flows_fileName)

ms_flows_buffer = ms_flows.copy()
ms_flows_buffer['geometry'] = ms_flows.buffer(10) # 10 meter buffer around ms flow lines
ms_flows_buffer['union'] = 'union' # add dummy column to dissolve all geometries into one
ms_flows_buffer_search = ms_flows_buffer.dissolve(by='union').geometry[0]  # take the single union geometry
#ms_split_flows_subset = gpd.sjoin(fr_split_flows, ms_flows_buffer, how='left', op='within')

ms_dem_flows_subset = fr_dem_flows[fr_dem_flows.within(ms_flows_buffer_search) | (fr_dem_flows.crosses(wbd_huc_boundary))]# & fr_dem_flows.touches(ms_flows_buffer_search))]

ms_dem_flows_subset.to_file(ms_dem_flows_fineName,driver='GPKG',index=False)
