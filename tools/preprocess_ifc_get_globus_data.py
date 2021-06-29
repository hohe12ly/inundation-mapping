#!/usr/bin/env python
from pathlib import Path
import re
import pandas as pd
import subprocess
import argparse
import time
import fileinput
from time import localtime, strftime

#Global variables
SOURCE_EP = 
DEST_EP = 
WORKSPACE = Path('/Path/to/workspace')

def write_batch_files(HUC8, groups):
    #Create HUC8 column (so we download by HUC8) and then group by huc8
    group= groups.get_group(HUC8)
    ifc_paths = group.path.to_list()
    for ifc_path in ifc_paths:
        #Set paths to files/directories that need to be downloaded
        model_dir = f'{ifc_path}/Hydraulics'
        grid_dir = f'{ifc_path}/Hydraulics/Hydraulics'
        
        #Test whether gdb or mdb are present

        bashCommand = f"globus ls {SOURCE_EP}:/{grid_dir} --filter ~Hydraulics.*db --jmespath DATA[*].[name] --format unix"
        process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
        output, error = process.communicate()
        db_files = output.decode("utf-8").splitlines()
        if 'Hydraulics.mdb' in db_files:
	        db_path = f"{grid_dir}/Hydraulics.mdb {grid_dir}/Hydraulics.mdb\n"
        elif 'Hydraulics.gdb' in db_files:
	        db_path = f"--recursive {grid_dir}/Hydraulics.gdb {grid_dir}/Hydraulics.gdb\n"
                        
        #Get models
        reach = model_dir.split('/')[-2]
        bashCommand = f"globus ls {SOURCE_EP}:/{model_dir} --filter ~{reach}.* --jmespath DATA[*].[name,type,size] --format unix"
        process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
        output, error = process.communicate()
        model_files = output.splitlines()
        model_files = [a.decode("utf-8") for a in model_files]
        model_files = [re.split(r'\t+', a) for a in model_files]
        model_files_df = pd.DataFrame(model_files, columns = ['name','type','size'])
        model_files_df = model_files_df[model_files_df.type == 'file'] 
        model_file_names = model_files_df.name.to_list()
               
        #Get depth grids
        bashCommand = f"globus ls {SOURCE_EP}:/{grid_dir} --jmespath DATA[*].[name] --format unix"
        process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
        output, error = process.communicate()
        grid_files = output.splitlines()
        grid_files = [a.decode("utf-8") for a in grid_files]
        depth_grids = [a for a in grid_files if a.lower().startswith('dp')]
    
        #Get LiDAR
        bashCommand = f"globus ls {SOURCE_EP}:/{ifc_path} --jmespath DATA[*].[name] --format unix"
        process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
        output, error = process.communicate()    
        lidar_files = output.splitlines()    
        lidar_files = [a.decode("utf-8") for a in lidar_files]
        dem_grids = [a for a in lidar_files if a.lower().startswith('1m')]   
           
        #Write Batch File
        BATCH_FILE = WORKSPACE/f'{ifc_path}/batch.txt'
        BATCH_FILE.parent.mkdir(parents = True, exist_ok = True)
        with open(BATCH_FILE,'w') as f:
            #Geodatabase
            f.write('#Geodatabase\n')
            f.write(db_path)
            #Model Files
            f.write('#Model\n')
            for model_file in model_file_names:
                f.write(f"{model_dir}/{model_file} {model_dir}/{model_file}\n")
            #Depth Grids
            f.write('#DepthGrids\n')
            for depth_grid in depth_grids:
                if '.' in depth_grid:
                    f.write(f"{grid_dir}/{depth_grid} {grid_dir}/{depth_grid}\n")
                else:
                    f.write(f"--recursive {grid_dir}/{depth_grid} {grid_dir}/{depth_grid}\n")
            #LiDAR
            f.write('#DEM\n')
            for dem in dem_grids:
                if '.' in dem:
                    f.write(f"{ifc_path}/{dem} {ifc_path}/{dem}\n")
                else:
                    f.write(f"--recursive {ifc_path}/{dem} {ifc_path}/{dem}\n")            

def get_globus(BATCH_FILE):         
    #Get Data              
    bashCommand = f'globus transfer {SOURCE_EP}:/ {DEST_EP}:{WORKSPACE} --batch' 
    batch_input = open(BATCH_FILE,"r")
    process = subprocess.Popen(bashCommand.split(), stdin = batch_input, stdout = subprocess.PIPE)
    output, error = process.communicate()
    batch_input.close()    
    taskid = output.splitlines()[1].decode("utf-8").split(':')[1].strip()    
    message = f'{taskid} will retrieve {BATCH_FILE.parent.name}'
    print(message)
    return message


def get_usace(usace_df):
    #Example to transfer USACE data
    #usace has small sizes, will download entire contents of these folders
    usace_paths = usace_df.path.to_list()
    for path in usace_paths:
        bashCommand = f"globus transfer {SOURCE_EP}:/{path} {DEST_EP}:{WORKSPACE}/{path} --recursive"
        process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
        output, error = process.communicate()


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description = 'Plot and aggregate statistics for benchmark datasets (BLE/AHPS libraries)')
    parser.add_argument('-huc','--huc', help = 'List of HUCs to process', nargs = '+', default = [], required = True)
    # Extract to dictionary and assign to variables
    args = vars(parser.parse_args())
    HUC_LIST = args['huc']
    
    
    #Get all files on GLOBUS (2 layers deep (hucs/reaches))
    print('Finding files...')
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
    
    #USACE
    usace_df = hr_df[hr_df.path.str.contains('_USACE') & ~ hr_df.path.str.endswith('_IFC')].copy()    
    #get_usace(usace_df)
    
    #IFC
    ifc_df = hr_df[~hr_df.path.str.contains('_USACE')].copy()
    ifc_df['huc8'] = ifc_df.path.str.split('/').str[0]
    huc_groups = ifc_df.groupby('huc8')
    
    #Iterate through the user specified HUCS
    for HUC in HUC_LIST:
        print(f"Working on HUC {HUC}")
        #Write Batch File
        print('Writing batch files...')
        write_batch_files(HUC, huc_groups)
        
        #Append files to single batch file
        print('Appending batch files...')
        batch_files = list((WORKSPACE / HUC).rglob('batch.txt')) 
        master_batch_file = WORKSPACE / HUC / 'all_batch.txt'
        with open(master_batch_file, 'w') as file:
            input_lines = fileinput.input(batch_files)
            file.writelines(input_lines)
        
        #Fetch GLOBUS entire HUC
        print('Fetching data...')
        message = get_globus(master_batch_file)
        with open(WORKSPACE/'globus_transfer_info.txt', 'a+') as file:
            file.write(f'{message} transfer time: {strftime("%Y-%m-%d %H:%M:%S%p", localtime())}\n')