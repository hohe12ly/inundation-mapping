import geopandas as gpd


def write(self,fileName,layer=None):

    """ Gets driver Name from file extension for Geopandas writing """

    driverDictionary = {'.gpkg' : 'GPKG','.geojson' : 'GeoJSON','.shp' : 'ESRI Shapefile'}
    driver = driverDictionary[splitext(fileName)[1]]
    
    self.to_file(fileName, driver=driver, layer=layer, index=False)

    return(driver)


def nhd_plus_hr_merge_vaas(nhd_plus_hr_flow_lines,nhd_plus_hr_vaas,on='NHDPlusID',what='all',nhd_plus_hr_flow_lines_layer_name=None,nhd_plus_hr_vaas_layer_name=None,output_file=None,output_layer_name=None):

    """ Merges two vector files on an attribute """

    # load nhd_plus_hr_flow_lines
    if isinstance(nhd_plus_hr_flow_lines,str):
        nhd_plus_hr_flow_lines = gpd.read_file(nhd_plus_hr_flow_lines,layer=nhd_plus_hr_flow_lines_layer_name)
    elif isinstance(nhd_plus_hr_flow_lines,gpd.GeoDataFrame):
        pass
    else:
        raise TypeError('Pass nhd_plus_hr_flow_lines argument as filepath or GeoDataframe')

    # load vaas
    if isinstance(nhd_plus_hr_vaas,str):
        nhd_plus_hr_vaas = gpd.read_file(nhd_plus_hr_vaas,layer=nhd_plus_hr_vaas_layer_name)
    elif isinstance(nhd_plus_hr_vaas,gpd.GeoDataFrame):
        pass
    else:
        raise TypeError('Pass nhd_plus_hr_vaas argument as filepath or GeoDataframe')
    
    # merging
    if isinstance(what,str):
        if what == 'all':
            nhd_plus_hr_flow_lines = nhd_plus_hr_flow_lines.merge(nhd_plus_hr_vaas,on=on, how='inner')
    if isinstance(what,list):
        nhd_plus_hr_flow_lines = nhd_plus_hr_flow_lines.merge(nhd_plus_hr_vaas[what],on=on, how='inner')
    else:
        raise TypeError('Pass what argument as list of VAAs to merge or "all" str to merge all VAAs')

    if output_file is not None:
        nhd_plus_hr_flow_lines.to_file(output_file,driver=getDriver(output_file),layer=output_layer_name,index=False)

    return(nhd_plus_hr_flow_lines)


def derive_stream_branches(self,toNode_attribute='toNode',fromNode_attribute='fromNode',order_attribute='order_',branch_id_attribute='branchID'):
    return(None)


def buffer_stream_branches(self,buffer_distance,branch_id_attribute='LevelPathI'):
    """ Buffers stream branches by distance """

    new_bids = [None] *len(self) ; new_geoms = new_bids.copy()
    i=0

    for bid in self[branch_id_attribute].unique():
        subset_table = self[self[branch_id_attribute]==bid]
        for _,row in subset_table.iterrows():
            new_geoms[i] = row[self.geometry.name].buffer(buffer_distance)
            new_bids[i] = row[branch_id_attribute]
            i += 1

    buffered_stream_branches = gpd.GeoDataFrame( { branch_id_attribute : new_bids,'geometry': new_geoms},
                                                   crs=self.crs,
                                                   geometry=self.geometry.name)

    return(buffer_stream_branches)

