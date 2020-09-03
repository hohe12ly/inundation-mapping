import wget
import rasterio
from tqdm import tqdm
import os

download_list='/data/inputs/ned665_20200903_082638.txt'
proj='PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
out_dir='/data/temp/lidar'

# download
with open(download_list) as f:
    for i,path in tqdm(enumerate(f)):
        
        wget.download(url=path,out=out_dir)
        file_name=os.path.basename(path)
        base_file,extension=os.path.splitext(file_name)


        lidar=rasterio.open(os.path.join(out_dir,file_name))
        lidar_proj = rasterio.shutil.copy_files(path,path + base_file + '_proj.' + extension)
        
        lidar_proj=rasterio.warp.reproject(lidar,lidar_proj,resampling='Resampling.bilinear',dst_crs=proj)
        lidar.close()
        
        lidar_window=rasterio.windows.get_data_window(rasterio.lidar_proj.read(1))
        
        if not rasterio.windows.intersect((lidar_window,)):
            lidar_proj.close()
            os.remove(lidar_proj)
            os.remove(lidar)

