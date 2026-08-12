#!/bin/bash
# this is to be executed on orcondo or oramd
# cmd: ./gpkg2fgdb.sh 12030201 /home/yxl/test/fim3/region6 /home/yxl/test/fim3/region6
# source soft env first

# orcondo: 
#source $HOME/sw/softenv
# oramd: 
source /srv/sw/softenv

h=$1
[ -z "$h" ] && echo "usage: $0 huccode input_dir output_dir" && exit 0

ddir="$2"
zipfile=`readlink -f $ddir/$h.zip`
[ ! -f "$zipfile" ] && echo "input $zipfile not accessible" && exit 0

odir="$3"
[ ! -d "$odir" ] && echo "output dir $odir not accessible" && exit 0

tdir=/tmp/gpkg2fdgb.$RANDOM
mkdir -p $tdir
cd $tdir
unzip  $zipfile ${h}/demDerived_reaches_split_filtered_addedAttributes_crosswalked.gpkg ${h}/gw_catchments_reaches_filtered_addedAttributes_crosswalked.gpkg ${h}/nhd_headwater_points_subset.gpkg ${h}/nwm_lakes_proj_subset.gpkg ${h}/nwm_subset_streams.gpkg $h/NHDPlusBurnLineEvent_subset.gpkg $h/wbd_buffered.gpkg $h/nwm_catchments_proj_subset.gpkg # ${h}/wbd8_clp.gpkg unknown geometry type

cd $h
for gpkg in `ls *.gpkg`
do 
  fn=`basename $gpkg .gpkg`
  ogr2ogr -f FileGDB $fn.gdb $gpkg -mapFieldType Integer64=Integer
  zip -q -r $fn.gdb.zip $fn.gdb
done

zip -q $odir/${h}-filegdb.zip *.zip

cd $HOME
rm -fr $tdir

echo "done. output is $odir/${h}-filegdb.zip"

