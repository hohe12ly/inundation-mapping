
import os
import argparse
from multiprocessing import Pool
import multiprocessing
import pandas as pd
from shapely.geometry import Polygon



from tools_shared_variables import TEST_CASES_DIR, MAGNITUDE_DICT, FR_BENCHMARK_CATEGORIES, AHPS_BENCHMARK_CATEGORIES


def sp_gdal_polygonize(args):
    
    benchmark_tif, output_dir, benchmark_type, flows_csv = args[0], args[1], args[2], args[3]
    
    # Open flows_csv and extract flow value and include in output_shapefile name.
    flows_df = pd.read_csv(flows_csv)
    median_flow = str(flows_df['discharge'].median())    
    
    file_handle = os.path.split(benchmark_tif)[1].replace('.tif', '_' + benchmark_type + '_' + median_flow + '_cms.shp')    
    output_shapefile = os.path.join(output_dir, file_handle)
    
    os.system('gdal_polygonize.py -8  {tif} {shape} -q'.format(tif=benchmark_tif, shape=output_shapefile))
    



def pointify_benchmark(output_dir):
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    # Map paths to benchmark tifs.
    test_cases_dir_list = os.listdir(TEST_CASES_DIR)
    
    procs_list = []
    
    for test_case_dir in test_cases_dir_list:
        if "test_cases" in test_case_dir:
            benchmark_type = test_case_dir.split('_')[0]
            validation_data_dir = os.path.join(TEST_CASES_DIR, test_case_dir, 'validation_data_' + benchmark_type)
            validation_data_huc_list = os.listdir(validation_data_dir)
            for huc in validation_data_huc_list:
                try:
                    int(huc)
                except: 
                    continue
                
                magnitude_list = MAGNITUDE_DICT[benchmark_type]
                
                if benchmark_type in AHPS_BENCHMARK_CATEGORIES:
                    lid_dir_list = os.listdir(os.path.join(validation_data_dir, huc))
                    for lid in lid_dir_list:
                        for magnitude in magnitude_list:
                            benchmark_tif = os.path.join(validation_data_dir, huc, lid, magnitude, 'ahps_' + lid + '_huc_' + huc + '_extent_' + magnitude + '.tif')
                            flows_csv = os.path.join(validation_data_dir, huc, lid, magnitude, 'ahps_' + lid + '_huc_' + huc + '_flows_' + magnitude + '.csv')
                            
                            if os.path.exists(benchmark_tif):
                                procs_list.append([benchmark_tif, output_dir, benchmark_type, flows_csv])
                
                if benchmark_type in FR_BENCHMARK_CATEGORIES:
                    for magnitude in magnitude_list:
                        benchmark_tif = os.path.join(validation_data_dir, huc, magnitude, benchmark_type + '_huc_' + huc + '_extent_' + magnitude + '.tif')
                        flows_csv = os.path.join(validation_data_dir, huc, magnitude, benchmark_type + '_huc_' + huc + '_flows_' + magnitude + '.csv')
                        
                        if os.path.exists(benchmark_tif):
                            procs_list.append([benchmark_tif, output_dir, benchmark_type, flows_csv])
                            
    job_number = multiprocessing.cpu_count() - 1
    print("Multiprocessing with " + str(job_number) + " jobs...")    
    
    with Pool(processes=job_number) as pool:
        pool.map(sp_gdal_polygonize, procs_list)


if __name__ == '__main__':
    
     # Parse arguments.
    parser = argparse.ArgumentParser(description='Pointify benchmark data.')
    parser.add_argument('-o','--output-dir',help='Name of directory to output results to',required=True)
    
    args = vars(parser.parse_args())
    
    pointify_benchmark(args['output_dir'])
    
    # Loop through all benchmark files
        # Run gdal_polygonize.py
        # Open the shapefiles, delete all 0 polygons
        # Close shapefile?
        # Convert polygon to point (make new point shapefile)
        # Open point shapefile, delete all odd-numbered rows
    