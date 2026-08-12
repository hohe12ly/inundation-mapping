#!/bin/bash
source /srv/sw/env.fim3

# 2021:
# sudo dnf install freetype-devel
# sudo dnf --enablerepo epel install wxGTK3-devel  mesa-libGL-devel mesa-libGLU-devel libXmu-devel libzstd-devel

# 202201:
# conda install -c conda-forge python=3.7 blas lapack fftw freetype

cd /srv/software/grass-7.8.6
# oramd:
swdir=/srv/sw/conda.evns/fim3
# libtiff must use conda's, otherwise, m.nviz.image will not compile
# --with-tiff --with-tiff-includes=$swdir/include --with-tiff-libs=$swdir/lib \

./configure --prefix=/srv/sw/grass78 \
 --with-cxx  --enable-largefile --with-openmp \
 --with-tiff --with-tiff-includes=$swdir/include --with-tiff-libs=$swdir/lib \
 --with-proj-share=$swdir/share/proj --with-proj-includes=$swdir/include --with-proj-libs=$swdir/lib \
 --with-geos=$swdir/bin/geos-config \
 --with-netcdf=$swdir/bin/nc-config  \
 --with-gdal=$swdir/bin/gdal-config \
 --with-sqlite --with-sqlite-includes=$swdir/include --with-sqlite-libs=$swdir/lib \
 --with-cairo --with-cairo-includes=$swdir/include/cairo --with-cairo-libs=$swdir/lib/cairo --with-cairo-ldflags=-lfontconfig \
 --with-freetype --with-freetype-includes=$swdir/include/freetype2 --with-freetype-libs=$swdir/lib \
 --with-blas --with-blas-includes=$swdir/include --with-blas-libs=$swdir/lib \
 --with-lapack --with-lapack-includes=$swdir/include --with-lapack-libs=$swdir/lib \
 --with-fftw --with-fftw-includes=$swdir/include --with-fftw-libs=$swdir/lib \
 --with-wxwidgets=/usr/bin/wx-config \
 --without-mysql \
 --without-postgres \
 --without-odbc 

# --with-freetype --with-freetype-includes=/usr/include/freetype2 \
# --with-blas --with-blas-includes=/usr/include/atlas-x86_64-base/ \
# --with-lapack --with-lapack-includes=/usr/include/atlas-x86_64-base/ \

# make -j16
# make install
# add in /srv/sw/env.fim3:
#  export GRASSBIN=/srv/sw/grass78/bin/grass78
#  export PATH=/srv/sw/grass78/bin:$PATH
#  export LD_LIBRARY_PATH=/srv/sw/grass78/grass78/lib:$LD_LIBRARY_PATH

# test:
# source /srv/sw/env.fim3
# mkdir /tmp/grassdata
# grass78 -c EPSG:5070 /tmp/grassdata/location --text
# in GRASS terminal:
# GRASS 7.8.6 (location):/srv/software/grass-7.8.6 > cd /srv/software/grass-7.8.6
# GRASS 7.8.6 (location):/srv/software/grass-7.8.6 > python -m grass.gunittest.main --location /tmp/grassdata/location --location-type nc
# GRASS 7.8.6 (location):/srv/software/grass-7.8.6 > r.import input=/dev/shm/fim3/outputs/12040201/12040201/dem_meters.tif output=dem # tif will be imported to /tmp/grassdata/location/PERMANENT/fcell/dem
