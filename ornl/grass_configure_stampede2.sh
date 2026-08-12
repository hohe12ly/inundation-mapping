#!/bin/bash
# incomplete. not working

source $WORK/env.conda

cd $WORK/software/grass-7.8.6
# stampede2:
swdir=$WORK/anaconda3/envs/pangeo
./configure --prefix=$WORK/sw/grass78
  --with-freetype \
  --with-cairo --with-cairo-ldflags=-lfontconfig \
  --with-tiff \
  --with-fftw --with-fftw-includes=$TACC_FFTW3_INC --with-fftw-libs=$TACC_FFTW3_LIB \
  --with-proj --with-proj-share=$swdir/share/proj --with-proj-includes=$swdir/include --with-proj-libs=$swdir/lib \
  --with-geos=$swdir/bin/geos-config  \
  --with-netcdf=$swdir/bin/nc-config  \
  --with-gdal=$swdir/bin/gdal-config \
  --with-sqlite --with-sqlite-includes=$swdir/include --with-sqlite-libs=$swdir/lib \
  --with-cxx --enable-largefile --with-python --with-openmp 
  --without-ffmpeg --without-mysql --without-postgres --without-odbc

# old grass compile in 2019:
# ./configure --with-freetype --with-freetype-includes=$WORK/sw/freetype2/include/freetype2 --with-freetype-libs=$WORK/sw/freetype2/lib --with-cairo --with-cairo-includes=$WORK/sw/cairo/include --with-cairo-libs=$WORK/sw/cairo/lib --with-cairo-ldflags=-lfontconfig --with-tiff --with-tiff-includes=$WORK/sw/tiff/include --with-tiff-libs=$WORK/sw/tiff/lib --with-fftw --with-fftw-includes=$TACC_FFTW3_INC --with-fftw-libs=$TACC_FFTW3_LIB  --with-proj-includes=$WORK/sw/proj/include --with-proj-libs=$WORK/sw/proj/lib --with-proj-share=$WORK/sw/proj/share --with-cxx --enable-largefile --with-gdal=$WORK/sw/gdal/bin/gdal-config --with-sqlite --with-sqlite-includes=$WORK/sw/sqlite/include --with-sqlite-libs=$WORK/sw/sqlite/lib --with-python --with-geos=$WORK/sw/geos/bin/geos-config --with-netcdf=$WORK/sw/netcdf4/bin/nc-config --prefix=$WORK/sw/grass

# make -j16
# make install
# add in /srv/sw/env.fim3:
#  export PATH=/srv/sw/grass78/bin:$PATH
#  export LD_LIBRARY_PATH=/srv/sw/grass78/grass78/lib:$LD_LIBRARY_PATH
