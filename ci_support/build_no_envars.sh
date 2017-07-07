#!/usr/bin/env bash

source ci_support/pre_build.sh

prefix=$HOME/miniconda2/
sp_dir=$prefix/lib/python2.7/site-packages

cmake $src_dir   -DENABLE_CONDA=OFF \
                    -DNUMPY_INCLUDE_PATH="$sp_dir/numpy/core/include" \
                    -DBoost_INCLUDE_DIR="$prefix/include" \
                    -DBoost_LIBRARIES="$prefix/lib/libboost_python.${suffix}" \
                    -DFFTW3_INCLUDE_PATH="$prefix/include" \
                    -DFFTW3_LIBRARY="$prefix/lib/libfftw3f.${suffix}" \
                    -DFFTW3d_INCLUDE_PATH="$prefix/include" \
                    -DFFTW3d_LIBRARY="$prefix/lib/libfftw3.${suffix}" \
                    -DFREETYPE_INCLUDE_DIRS="$prefix/include/freetype2" \
                    -DFREETYPE_LIBRARIES="$prefix/lib/libfreetype.${suffix}" \
                    -DFTGL_INCLUDE_PATH="$prefix/include" \
                    -DFTGL_LIBRARY="$prefix/lib/libftgl.${suffix}" \
                    -DGSL_CBLAS_INCLUDE_PATH="$prefix/include" \
                    -DGSL_CBLAS_LIBRARY="$prefix/lib/libgslcblas.${suffix}" \
                    -DGSL_INCLUDE_PATH="$prefix/include" \
                    -DGSL_LIBRARY="$prefix/lib/libgsl.${suffix}" \
                    -DHDF5_INCLUDE_PATH="$prefix/include" \
                    -DHDF5_LIBRARY="$prefix/lib/libhdf5.${suffix}" \
                    -DJPEG_INCLUDE_PATH="$prefix/include" \
                    -DJPEG_LIBRARY="$prefix/lib/libjpeg.${suffix}" \
                    -DPNG_PNG_INCLUDE_DIR="$prefix/include" \
                    -DPNG_LIBRARY_RELEASE="$prefix/lib/libpng.${suffix}" \
                    -DPYTHON_INCLUDE_PATH="$prefix/include/python2.7" \
                    -DPYTHON_LIBRARY="$prefix/lib/libpython2.7.${suffix}" \
                    -DTIFF_INCLUDE_PATH="$prefix/include" \
                    -DTIFF_LIBRARY="$prefix/lib/libtiff.${suffix}" \
                    -DZLIB_LIBRARY="$prefix/lib/libz.${suffix}"

make
make install

export PREFIX="${HOME}"/EMAN2
export SP_DIR="${PREFIX}"/lib

export PYTHONPATH="$SP_DIR:$PYTHONPATH"

ln -s $PREFIX/bin/e2version.py $SP_DIR/e2version.py
ln -s $PREFIX/bin/sxgui.py     $PREFIX/bin/sphire
ln -s $PREFIX/bin/sx.py        $PREFIX/bin/sparx

make test-verbose

source $src_dir/ci_support/post_build.sh
