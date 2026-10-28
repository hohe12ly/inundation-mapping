#!/bin/bash
# the workflow steps to crop&reproj DSM, create surface healded HAND, and generate flood polygon stacks
# input: huc code

source $HOME/env.inunmap
h=$1

echo "[`date`] flood stack generation for HUC $h"

t00=`date +%s`
t0=`date +%s`
echo "[`date`] Step 1: crop and reproject road DSM"
state_road_dsm=/mnt/orc17tb/sunsphere/test/texas_road_dsm/o/texas_road_dsm3m.tif
handroot=/mnt/CADES_CEPH/pin2flood/prod
[ ! -d "$handroot/$h" ] && echo "HAND for HUC $h does not exist." && exit 1

idir=$handroot/$h
oroot=/mnt/CADES_CEPH/pin2flood/floodstacks
odir=$oroot/$h
mkdir -p $odir

python gdalwarp_like.py --input_raster $state_road_dsm --reference_raster $idir/dem_meters.tif --output_raster $odir/${h}_road_dsm.tif;
t1=`date +%s`
echo "[`date`] TIMEStep 1 in seconds: `expr $t1 \- $t0`"

echo "[`date`] Step 2: generate surface healed HAND"
t0=`date +%s`
gdal_calc.py --type=Float32 --co "COMPRESS=LZW" --co "BIGTIFF=YES" --co "TILED=YES" --hideNoData --NoDataValue="-9999" -A $odir/${h}_road_dsm.tif -B $idir/dem_meters.tif --calc="(A-B)*(A>=B)*((A-B)<9900)" --outfile=$odir/${h}_dsm_dem_diff.tif
t1=`date +%s`
gdal_calc.py --type=Float32 --co "COMPRESS=LZW" --co "BIGTIFF=YES" --co "TILED=YES" --NoDataValue="-9999"  -A $idir/rem_zeroed_masked_healed.tif -B $odir/${h}_dsm_dem_diff.tif --calc="A+B" --outfile=$odir/${h}_surface_healed_hand.tif
echo "[`date`] TIMEStep 2 in seconds: `expr $t1 \- $t0`"

echo "[`date`] Step 3: generate flood polygon stack"
t0=`date +%s`
python floodstack2.py $h
t1=`date +%s`
echo "[`date`] TIMEStep 3 in seconds: `expr $t1 \- $t0`"

t11=`date +%s`
echo "[`date`] TIMETOTAL in seconds: `expr $t11 \- $t00`"
