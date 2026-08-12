#!/bin/bash
# create a color tif overview for a HUC8 HAND
# cmd: 
# hand-overview.sh 12070101 /home/yxl/test/fim3/region6/12070101.zip 12070101/rem_zeroed_masked.tif /home/yxl/test/fim3/region6-overview
source /srv/sw/env.fim3

curdir=`pwd`
ratio=20 # how many times smaller
colorfile=/srv/nfie-floodmap/test/HAND-blues.clr

huc="$1"
zip=`readlink -f $2`
fpath="$3"
odir=`readlink -f $4`

tdir=/dev/shm/fim3/overview
mkdir -p $tdir
cd $tdir

echo "### extracting $fpath from $zip"
[ ! -f "$fpath" ] && unzip $zip $fpath
[ ! -f "$fpath" ] && echo "$fpath does not exist in $zip" && exit 0

echo "### read size info"
sizeline=`gdalinfo $fpath | grep -E "^Size is " | head -n 1`
w=`echo "$sizeline" | awk -F, '{print $1}' | awk '{print $NF}'`
h=`echo "$sizeline" | awk '{print $NF}'`
echo "raster size is $w x $h"
ow=`expr $w \/ $ratio`
oh=`expr $h \/ $ratio`
echo "overview size will be $ow x $oh"

echo "### reduce input tif to 1/$ratio "
handtif=$odir/ov_hand_${huc}.tif
gdalwarp -r bilinear -ts $ow $oh $fpath $handtif

# gdal, when overlapping color tif together, does not respect alpha,
# so have to generate hand overview vrt, then color it.
#colordd=$odir/ov_${huc}.tif
#echo "### create color tif $colordd"
#gdaldem color-relief $temptif $colorfile $colordd -of GTiff -alpha

[ -f $fpath ] && rm $fpath && rmdir `dirname $fpath`
#[ -f $temptif ] && rm $temptif
cd $curdir
