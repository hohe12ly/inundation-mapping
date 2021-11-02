#!/usr/bin/env python3
from grass_session import Session
import os, shutil 
import grass.script as gscript
import rasterio
import numpy as np
import argparse
from rasterstats import zonal_stats
import geopandas as gpd

def calc_raster_difference(dem_source_filename, dem_conditioned_filename, notch_threshold):

    # Get directory
    directory = os.path.dirname(dem_source_filename)

    # Subtract conditioned dem from source
    dem_source = rasterio.open(dem_source_filename)
    profile = dem_source.profile # get profile for new raster creation later on
    dem_source = dem_source.read(1)
    dem_conditioned = rasterio.open(dem_conditioned_filename).read(1)

    thalweg_difference = dem_source - dem_conditioned

    # Calculate standard deviation if no threshold is supplied
    if not notch_threshold:

        thalweg_difference[thalweg_difference == 0] = np.nan
        notch_threshold = np.nanmean(thalweg_difference) + np.nanstd(thalweg_difference)
        print(f"No notch threshold supplied. Using {notch_threshold} calculated as one standard deviation above the mean.")

    # Filter raster to the notch threshold
    thalweg_difference[thalweg_difference <= notch_threshold] = np.nan

    # Save raster
    thalweg_difference_filename = os.path.join(directory, "dem_thalweg_difference.tif")
    with rasterio.open(thalweg_difference_filename, "w", **profile) as out_raster:
        out_raster.write(thalweg_difference, 1)

    return thalweg_difference_filename

def r_to_vect(input_raster, grass_workspace, min_segment_len, clip2catchments):
    '''
    Runs the r.to.vect GRASS gis tool which given an input raster will produce an output proximity (or distance) and euclidian allocation tool.

    Parameters
    ----------
    input_raster : STR
        Path to input raster. For example, see dem_thalweg_difference.tif.
    grass_workspace : STR
        Path to TEMPORARY directory to store intermediate GRASS data. This directory is deleted upon completion of this function.
    min_segment_len : INT
        Minimum segment length in meters to include in the output. By default, all segments are returned.
    clip2catchments : STR
        Path to HUC8 gpk. For example, see wbd.gpkg.
    
    Returns
    -------
    thalweg_notch.gpkg in the same directory as the input raster

    '''
    
    # Define parent directory of input raster and get input raster name
    input_raster_directory = os.path.dirname(input_raster)
    input_raster_name = os.path.splitext(os.path.basename(input_raster))[0]
    
    # Set up variables for use in GRASS
    grass_gisdb = grass_workspace
    grass_location = 'temporary_location'
    grass_mapset = 'temporary_mapset'
    projected_file = input_raster

    # Start and close PERMANENT session.
    PERMANENT = Session()
    PERMANENT.open(gisdb=grass_gisdb, location=grass_location, create_opts=projected_file)
    PERMANENT.close()

    # Open a temporary session.
    temporary_session = Session()
    temporary_session.open(gisdb=grass_gisdb, location=grass_location, mapset=grass_mapset, create_opts=projected_file)
    
    #Import input raster into temporary session.
    imported_grass_raster = input_raster_name + '@' + grass_mapset
    gscript.run_command('r.in.gdal', input=input_raster, output=imported_grass_raster, quiet=True)

    # r.to.vect tool.
    notch_grass_vector = 'thalweg_notch_init@'+ grass_mapset
    gscript.run_command('r.to.vect', flags=None, input=imported_grass_raster, output=notch_grass_vector, type="line", column="notch_value", quiet=True)

    # Add segment_length column and caculate
    segment_len_col = 'notch_length'
    gscript.run_command('v.db.addcolumn', map=notch_grass_vector, columns=f'{segment_len_col} double precision')
    gscript.run_command('v.to.db', map=notch_grass_vector, option='length', type='line', columns=segment_len_col, units='meters', quiet=True)

    # Variable handler dict - because grass won't let you use the same names for multiple map elements
    grass_vars = grass_var_handler(notch_grass_vector, grass_mapset)

    # Filter out short line segments
    if min_segment_len:
        gscript.run_command('v.extract', input=grass_vars.input(), output=grass_vars.output(), where=f'{segment_len_col} >= {min_segment_len}', type='line', quiet=True, overwrite=True)
    
    # Clip to HUC8 boundary
    if clip2catchments:
        huc8_grass_vector = 'huc'
        gscript.run_command('v.in.ogr', input=clip2catchments, output=huc8_grass_vector, type='boundary')
#        gscript.run_command('v.overlay', ainput=grass_vars.input(), operator='and', binput=huc8_grass_vector, output=grass_vars.output(), olayer='0,1,0', overwrite=True, quiet=True) 
        gscript.run_command('v.overlay', ainput=grass_vars.input(), operator='and', binput=huc8_grass_vector, output=grass_vars.output(), olayer='1,0,0', overwrite=True, quiet=True) 
    
    # Rename the mapset for output
    gscript.run_command('g.rename', vector=f"{grass_vars.input()},thalweg_notch", quiet=True)

    # Export notch vector
    output_notch_filename = os.path.join(input_raster_directory,"thalweg_notch.gpkg")
    gscript.run_command('v.out.ogr', flags='c', input="thalweg_notch", output=output_notch_filename, format='GPKG', quiet=True, type='line', overwrite=True)

    # Close down temporary session and remove temporary workspace.
    temporary_session.close()
    shutil.rmtree(grass_gisdb)

    return output_notch_filename

def apply_zonal_stats(notch_filename, huc8):

    difference_raster_directory = os.path.dirname(notch_filename)
    difference_raster = os.path.join(difference_raster_directory, "dem_thalweg_difference.tif")
    stats = zonal_stats(notch_filename, difference_raster, stats="mean")
    means = [m["mean"] for m in stats]

    # Read in thalweg_notch.gpkg
    thalweg_notch_df = gpd.read_file(notch_filename)
    # Drop unnecessary columns
    thalweg_notch_df.drop(['cat', 'b_cat', 'b_S0', 'b_LengthKm', 'b_LakeID', 'b_areasqkm', 'b_min_thal_elev', 'b_med_thal_elev', 'b_max_thal_elev'], axis=1, inplace=True)
    # Set mean notch depths, huc8, and clean up column names
    thalweg_notch_df["huc8"] = huc8 #* len(thalweg_notch_df.index)
    thalweg_notch_df.a_notch_value = means
    thalweg_notch_df.set_axis(['original_notch_fid',        # Original id of the notch before it was segmented by catchment
                               'mean_notch_depth',          # Average depth of the notch over the length of the catchment
                               'total_notch_length',        # Length of the original notch before it was segmented by catchment
                               'HydroID',                   # HydroID from FIM (obtained from input catchments)
                               'From_Node',                 # obtained from input catchments
                               'To_Node',                   # obtained from input catchments
                               'NextDownID',                # obtained from input catchments
                               'feature_id',                # feature_id from National Water Model (obtained from input catchments)
                               'stream_order',              # obtained from input catchments
                               'geometry',
                               'huc8'],                     # assigned based on current huc
                               axis=1, inplace=True)
    # Overwrite thalweg_notch.gpkg
    thalweg_notch_df.to_file(notch_filename, driver='GPKG', layer='thalweg_notch')
    
    return

class grass_var_handler():
    '''Creates names for intermediate grass layers when there are multiple optional steps'''
    def __init__(self, input_var, grass_mapset):
        self.grass_mapset = grass_mapset
        self.input_var = input_var
        self.output_num = 0
    
    def input(self):
        return self.input_var

    def output(self):
        out_var = f"output_{self.output_num}@{self.grass_mapset}"
        self.output_num += 1
        self.input_var = out_var
        return out_var

def aggregate_notches(fim_dir):
    agg_thalweg_notch_filename = os.path.join(fim_dir, 'thalweg_notch.gpkg')
    for huc in os.listdir(fim_dir):
        print(huc)
        notch_gpkg = os.path.join(fim_dir, huc, "thalweg_notch.gpkg")
        if not os.path.isfile(notch_gpkg): continue

        notch_df = gpd.read_file(notch_gpkg)
        if not len(notch_df.index): continue

        print(f'..{len(notch_df.index)}')
        if os.path.isfile(agg_thalweg_notch_filename):
            notch_df.to_file(agg_thalweg_notch_filename,index=False, driver='GPKG', mode='a')
        else:
            notch_df.to_file(agg_thalweg_notch_filename,index=False, driver='GPKG')
        del notch_df
    return

if __name__ == '__main__':

    #Parse arguments
    parser = argparse.ArgumentParser(description='Calculate thalweg notch segments')
    parser.add_argument('-d', '--dem_meters', help='Input source dem, e.g. dem_meters.tif', required=True)
    parser.add_argument('-f', '--dem_thalweg_cond', help='Input thalweg conditioned dem, e.g. dem_thalwegCond.tif', required=True)
    parser.add_argument('-g', '--grass_workspace', help='Temporary GRASS workspace', required=True)
    parser.add_argument('-t', '--notch_threshold', help='Float in meters. Difference in thalweg conditioning in meters to detect as a notch. Defaults to one standard deviation above the mean', required=False)
    parser.add_argument('-l', '--segment_len', help='Integer in meters. Minimum notch segment length to filter output. If none specified, all segements will be returned', required=False)
    parser.add_argument('-c', '--clip_to_catchments', help='Clip output segments to catchments, e.g. wbd.gpkg', required=False)

    # Extract to dictionary and assign to variables.
    args = vars(parser.parse_args())

    # Rename variable inputs
    dem_meters = args['dem_meters']
    dem_thalweg_cond = args['dem_thalweg_cond']
    grass_workspace = args['grass_workspace']
    notch_threshold = args['notch_threshold']
    segment_len = args['segment_len']
    clip2catchments = args['clip_to_catchments']

    # Typecasting arguments
    if notch_threshold: notch_threshold = float(notch_threshold)
    if segment_len: segment_len = int(segment_len)
    huc8 = os.path.basename(os.path.dirname(grass_workspace))
    
    try:
        # Subtract conditioned thalweg from source dem
        print("Subtracting rasters")
        thalweg_difference_raster = calc_raster_difference(dem_meters, dem_thalweg_cond, notch_threshold)

        # Run r_to_vect
        print("Raster to lines")
        notch_filename = r_to_vect(thalweg_difference_raster, grass_workspace, segment_len, clip2catchments)

        # Compute mean notch depth for each segment
        if clip2catchments:
            apply_zonal_stats(notch_filename, huc8)
    except Exception as ex:
        print(ex)
        home_dir = os.path.dirname(os.path.dirname(grass_workspace))
        with open(os.path.join(home_dir, "logs", "errors.log"), "a") as error_file:
            error_file.write(huc8 + "\n")
        if os.path.isdir(grass_workspace):
            shutil.rmtree(grass_workspace, ignore_errors=True)