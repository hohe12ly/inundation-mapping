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

fr_split_flows_fileName   = sys.argv[1]
ms_flows_fileName         = sys.argv[2]
ms_split_flows_fineName   = sys.argv[3]

print('Loading data ...')
fr_split_flows = gpd.read_file(fr_split_flows_fileName)
ms_flows = gpd.read_file(ms_flows_fileName)

ms_flows_buffer = ms_flows.copy()
ms_flows_buffer['geometry'] = ms_flows.buffer(10) # 10 meter buffer around lines
ms_flows_buffer['union'] = 'union' # add dummy column to dissolve all geometries into one
ms_flows_buffer_search = ms_flows_buffer.dissolve(by='union').geometry[0]  # take the single union geometry
#ms_split_flows_subset = gpd.sjoin(fr_split_flows, ms_flows_buffer, how='left', op='within')

ms_split_flows_subset = fr_split_flows[fr_split_flows.within(ms_flows_buffer_search)]

ms_split_flows_subset.to_file(ms_split_flows_fineName,driver='GPKG',index=False)
