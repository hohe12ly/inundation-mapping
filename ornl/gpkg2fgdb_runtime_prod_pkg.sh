#!/bin/bash
# this is to be executed on orcondo or oramd
# cmd: ./gpkg2fgdb_runtime.sh 12030201 /home/yxl/test/fim3/region6 /home/yxl/test/fim3/region6
# source soft env first

# orcondo: 
#source $HOME/sw/softenv
# oramd: 
source /srv/sw/softenv

h=$1
[ -z "$h" ] && echo "usage: $0 huccode input_dir output_dir" && exit 0

ddir="$2"
[ ! -d "$ddir" ] && echo "runtime $ddir not accessible" && exit 0
[ ! -d "$ddir/$h" ] && echo "runtime huc $ddir/$h not accessible" && exit 0

odir="$3"
[ ! -d "$odir" ] && echo "output dir $odir not accessible" && exit 0

tdir=/tmp/gpkg2fdgb.$RANDOM
mkdir -p $tdir
cd $tdir

for gpkg in demDerived_reaches_split_filtered_addedAttributes_crosswalked.gpkg gw_catchments_reaches_filtered_addedAttributes_crosswalked.gpkg
do 
  fn=`basename $gpkg .gpkg`
  ogr2ogr -f FileGDB $fn.gdb $ddir/$h/$gpkg -mapFieldType Integer64=Integer
done

zip -q -r $odir/${h}/filegdb.zip *.gdb

cd $HOME
rm -fr $tdir

echo "done. output is $odir/${h}/filegdb.zip"

