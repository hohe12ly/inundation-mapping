#!/bin/bash

sdir=/srv/test/cahaba/ornl
rdir=/home/yxl/test/fim3
ddir=$rdir/tx-region6
odir=$rdir/region6-overview

cd $ddir
for z in `ls *.zip`; 
do
  echo "========== $z start ============"
  echo "================================"
  huc=`basename $z .zip`
  $sdir/hand-overview.sh $huc $ddir/$z $huc/rem_zeroed_masked.tif $odir
  echo "================================"
  echo "========== $z end ============"
done

echo "========= region 6 overview start =========="
for f in `ls $odir`:
do
  echo "$odir/$f"
done>$rdir/region6_ov_hand.filelist

source /srv/sw/env.fim3
gdalbuildvrt -input_file_list $rdir/region6_ov_hand.filelist $rdir/region6_ov_hand.vrt

sleep 5

colorfile=/srv/nfie-floodmap/test/HAND-blues.clr
colordd=$rdir/region6_ov_hand_color.tif
gdaldem color-relief $rdir/region6_ov_hand.vrt $colorfile $colordd -of GTiff -alpha
echo "========= region 6 overview end =========="
