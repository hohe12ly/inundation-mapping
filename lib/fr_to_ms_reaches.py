#!/usr/bin/env python3

'''
Description:
    1) buffer the MS demDerived_reaches.shp
    2) locate all FR reaches within MS buffer
    3) export FR-->MS subsetted reaches
'''

import sys
import geopandas as gpd

fr_split_flows_fileName   = sys.argv[1]
ms_flows_fileName         = sys.argv[2]
ms_split_flows_fineName   = sys.argv[3]

print('Loading data ...')
fr_split_flows = gpd.read_file(fr_split_flows_fileName)
ms_flows = gpd.read_file(ms_flows_fileName)

ms_flows_buffer['geometry'] = ms_flows.buffer(10) # 10 meter buffer around lines
ms_split_flows_gdf = gpd.sjoin(fr_split_flows, ms_flows_buffer, how='left', op='within')

ms_split_flows_gdf.to_file(ms_split_flows_fineName,driver='GPKG',index=False)