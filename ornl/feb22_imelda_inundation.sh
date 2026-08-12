#!/bin/bash
source /srv/test/cahaba/ornl/softenv
source /srv/test/cahaba/ornl/params_calibrated.env
source $srcDir/bash_functions.env
m2f="3.280839895" # meter to foot

datadir=/home/yxl/test/nwm-reanalyais/imelda/feb22
oepsg="2278" # TODO: in foot
tmpdir=/dev/shm/fim3/inun
mkdir -p $tmpdir

tt0=`date +%s`
for h in 12010005 12020003 12020006 12040103 12020007 12030203 12040201 12040203 12040202
do
  t0=`date +%s`
  for f_timestamp in `ls /home/yxl/test/nwm-reanalyais/imelda/feb22/o/12010005/2019-09-*.csv`
  do
    timestamp=`basename $f_timestamp .csv`

    ## flood depth grid
    #python /srv/test/cahaba/tools/inundation.py -r $datadir/$h/rem_zeroed_masked.tif -c $datadir/$h/gw_catchments_reaches_filtered_addedAttributes.tif -b $datadir/$h/gw_catchments_reaches_filtered_addedAttributes_crosswalked.gpkg -t $datadir/$h/hydroTable.csv -f $datadir/o/$h/${timestamp}.csv -i $datadir/o/$h/${timestamp}_inun_bin.tif -p $datadir/o/$h/${timestamp}_inun_bin.gpkg -d $datadir/o/$h/${timestamp}_inun_depth.tif
    # convert to target projection
    gdalwarp -t_srs EPSG:$oepsg $datadir/o/$h/${timestamp}_inun_depth.tif $tmpdir/${h}_${timestamp}_inun_depth_$oepsg.tif
    # gdalwarp does not care about the vertical unit
    gdal_calc.py --quiet --overwrite --co "BIGTIFF=YES" --co "COMPRESS=LZW" -A $tmpdir/${h}_${timestamp}_inun_depth_$oepsg.tif --calc="A*$m2f" --NoDataValue=-9999 --outfile=$datadir/o/$h/${timestamp}_inun_depth_$oepsg.tif

    ## create water surface elev raster
    # to allow nodata participate gdal_calc
    gdal_translate -a_nodata -20000 $datadir/o/$h/${timestamp}_inun_depth.tif $tmpdir/${h}_${timestamp}_inun_depth_tmp.tif
    # calc water surface elev grid
    gdal_calc.py --quiet --overwrite --co "BIGTIFF=YES" -A $datadir/$h/dem_meters.tif -B $tmpdir/${h}_${timestamp}_inun_depth_tmp.tif --calc="(B>=0)*(A+B) + (B<0)*A" --NoDataValue=-9999 --outfile=$tmpdir/${h}_${timestamp}_welev.tif
    # convert to target projection
    gdalwarp -t_srs EPSG:$oepsg $tmpdir/${h}_${timestamp}_welev.tif $tmpdir/${h}_${timestamp}_welev_$oepsg.tif
    # gdalwarp does not care about the vertical unit
    gdal_calc.py --quiet --overwrite --co "BIGTIFF=YES" --co "COMPRESS=LZW" -A $tmpdir/${h}_${timestamp}_welev_$oepsg.tif --calc="A*$m2f" --NoDataValue=-9999 --outfile=$datadir/o/$h/${timestamp}_water_surfave_elev_$oepsg.tif

    ## cleanup
    rm $tmpdir/*.tif

    #exit 0
  done

  ## copy a dem there
  gdal_calc.py --quiet --overwrite --co "BIGTIFF=YES" --co "COMPRESS=LZW" -A $datadir/$h/dem_meters.tif --calc="A*$m2f" --NoDataValue=-9999 --outfile=$datadir/o/$h/dem_meters_$oepsg.tif

  ## copy a hydroTable and flowline geopackage
  cp $datadir/$h/hydroTable.csv $datadir/o/$h/
  ogr2ogr -t_srs EPSG:$oepsg $datadir/o/$h/demDerived_reaches_split_filtered_addedAttributes_crosswalked_$oepsg.gpkg $datadir/$h/demDerived_reaches_split_filtered_addedAttributes_crosswalked.gpkg 

  t1=`date +%s`
  echo "TIME $h `expr $t1 \- $t0` seconds"
  #exit 0
done
tt1=`date +%s`
echo "TIME TOTAL `expr $tt1 \- $tt0` seconds"
