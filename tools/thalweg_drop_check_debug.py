#!/usr/bin/env python3

import os
import sys
import geopandas as gpd
sys.path.append('/foss_fim/src')
from shapely.geometry import Point, LineString
import rasterio
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from collections import deque
from os.path import join
from multiprocessing import Pool
from utils.shared_functions import getDriver
from rasterio import features
from reachID_grid_to_vector_points import convert_grid_cells_to_points
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from thalweg_drop_check import compare_thalweg

if __name__ == '__main__':
    compare_thalweg(['/data/outputs/single_pixel_inundation/12090401', 
                 'derived', 
                 'all_points', 
                 '12090401', 
                 '/data/outputs/single_pixel_inundation/12090401/dem_meters.tif', 
                 '/data/outputs/single_pixel_inundation/12090401/dem_lateral_thalweg_adj.tif', 
                 '/data/outputs/single_pixel_inundation/12090401/dem_thalwegCond.tif', 
                 '/data/tools/thalweg_profile_comparison/fim_3_0_22_6/plots/profile_drop_plots_12090401_all_points_derived.png', 
                 '/data/tools/thalweg_profile_comparison/fim_3_0_22_6/spatial_layers/thalweg_elevation_changes_12090401_all_points_derived.gpkg', 
                 '/data/tools/thalweg_profile_comparison/fim_3_0_22_6/spatial_layers/thalweg_elevation_changes_12090401_all_points_derived.csv', 
                 '/data/outputs/single_pixel_inundation/12090401/flows_grid_boolean.tif'])