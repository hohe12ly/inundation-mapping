'''
this is a utility to crop and reproject a raster using the projection, bounding box, and dimensions of a reference raster. It is similar to gdalwarp, but it is not as general purpose.
also, the output raster should be geotiff with compression and bigtiff enabled.
'''
from osgeo import gdal
import os, sys, time, argparse
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Crop and reproject a raster using the projection, bounding box, and dimensions of a reference raster.')
    parser.add_argument('--input_raster', help='Input raster file path')
    parser.add_argument('--reference_raster', help='Reference raster file path')
    parser.add_argument('--output_raster', help='Output raster file path')
    args = parser.parse_args()

    t0 = time.perf_counter()

    input_raster = args.input_raster
    reference_raster = args.reference_raster
    output_raster = args.output_raster

    # usal gdal exception handling
    gdal.UseExceptions()

    # Open the input raster
    src_ds = gdal.Open(input_raster)
    if src_ds is None:
        print(f"Failed to open input raster: {input_raster}")
        sys.exit(1)

    # Open the reference raster
    ref_ds = gdal.Open(reference_raster)
    if ref_ds is None:
        print(f"Failed to open reference raster: {reference_raster}")
        sys.exit(1)

    # Get the projection and geotransform from the reference raster
    ref_proj = ref_ds.GetProjection()
    ref_geotransform = ref_ds.GetGeoTransform()
    ref_width = ref_ds.RasterXSize
    ref_height = ref_ds.RasterYSize

    # Create the output raster with the same projection and dimensions as the reference raster
    driver = gdal.GetDriverByName('GTiff')
    # options: compression=LZW, bigtiff=YES, tile size 256x256
    # nodata: carry over the nodata value from the input raster
    out_ds = driver.Create(output_raster, ref_width, ref_height, src_ds.RasterCount, gdal.GDT_Float32, options=['COMPRESS=LZW', 'BIGTIFF=YES', 'TILED=YES', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256'])
    
    if out_ds is None:
        print(f"Failed to create output raster: {output_raster}")
        sys.exit(1)

    out_ds.SetProjection(ref_proj)
    out_ds.SetGeoTransform(ref_geotransform)

    # Carry over the nodata value from the input raster (per band) and
    # initialize the output with it, so pixels outside the input coverage
    # (or sampled from nodata) come out as nodata instead of 0
    for band_idx in range(1, src_ds.RasterCount + 1):
        nodata = src_ds.GetRasterBand(band_idx).GetNoDataValue()
        if nodata is not None:
            out_band = out_ds.GetRasterBand(band_idx)
            out_band.SetNoDataValue(nodata)
            out_band.Fill(nodata)

    # Perform the reprojection and cropping; show progress
    print("Performing reprojection and cropping...")
    # print output projection (EPSG code) and dimensions of the output raster, as well as the bounding box, and resolution
    print(f"Output raster projection: {ref_proj}")
    print(f"Output raster dimensions: {ref_width} x {ref_height}")
    print(f"Output raster bounding box: {ref_geotransform[0]}, {ref_geotransform[3]}, {ref_geotransform[0] + ref_geotransform[1] * ref_width}, {ref_geotransform[3] + ref_geotransform[5] * ref_height}")
    print(f"Output raster resolution: {ref_geotransform[1]} x {ref_geotransform[5]}")   
    gdal.ReprojectImage(src_ds, out_ds, src_ds.GetProjection(), ref_proj, gdal.GRA_Bilinear)

    # Clean up
    src_ds = None
    ref_ds = None
    out_ds = None
    
    t1 = time.perf_counter()
    elapsed_time = t1 - t0

    print(f"Successfully created output raster: {output_raster} in {elapsed_time:.2f} seconds.")