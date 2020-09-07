import wget
import rasterio
import rasterio.shutil
import rasterio.warp
import rasterio.mask
import rasterio.crs
import fiona
from tqdm import tqdm
import os
from shapely.geometry import Polygon

download_list='/data/temp/lidar/ned665_20200903_082638.txt'
proj='PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.2572221010042,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_center",23],PARAMETER["longitude_of_center",-96],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
out_dir='/data/temp/lidar'
aoi = '/data/outputs/slope_1eddd72_dem_all_mannings_6/12090301/wbd.gpkg'
buffer_distance = 5000

# open aoi and convert to window
aoi = fiona.open(aoi)
poly = [Polygon(p['geometry']['coordinates'][0][0]) for p in iter(aoi)]
aoi.close()

# make proj
proj = rasterio.crs.CRS.from_string(proj)

# download
with open(download_list,mode='r') as f:
    total_len = sum(1 for _ in f)

with open(download_list,mode='r') as f:
    for url_path in tqdm(f,total=total_len):

        url_path = url_path.split('\n')[0]
        file_name=os.path.basename(url_path)
        base_file,extension=os.path.splitext(file_name)
        full_file_path = os.path.join(out_dir,file_name)

        if not os.path.isfile(file_name):
            wget.download(url=url_path,out=out_dir,bar=None)

        lidar=rasterio.open(full_file_path,'r')
       
        lidar_proj_full_file_path = os.path.join(out_dir,base_file + '_proj' + extension)
        lidar_proj_meta = lidar.meta.copy()
        lidar_proj_meta.update({'crs':proj})
        lidar_proj = rasterio.open(lidar_proj_full_file_path,'w',**lidar_proj_meta)

        lidar_proj_band,_ = rasterio.warp.reproject(rasterio.band(lidar,1),rasterio.band(lidar_proj,1),resampling=1)
        print(lidar_proj_band)
        lidar_proj.write(lidar_proj_band,1)

        lidar.close()
        rasterio.shutil.delete(full_file_path)
        
        try:
            rasterio.mask.raster_geometry_mask(lidar_proj,poly,crop=True)
        except ValueError:
            lidar_proj.close()
            rasterio.shutil.delete(lidar_proj_full_file_path)
        else:
            lidar_proj.close()


