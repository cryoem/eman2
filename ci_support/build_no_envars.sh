#!/usr/bin/env bash

source ci_support/pre_build.sh

prefix=$HOME/miniconda2/
sp_dir=$prefix/lib/python2.7/site-packages

cmake $src_dir   -DCMAKE_INSTALL_RPATH="$HOME/EMAN2/lib" \
                    -DNUMPY_INCLUDE_PATH="$sp_dir/numpy/core/include" \
                    -DBOOST_INCLUDE_PATH="$prefix/include" \
                    -DBOOST_LIBRARY="$prefix/lib/libboost_python.${suffix}" \
                    -DFFTW3_INCLUDE_PATH="$prefix/include" \
                    -DFFTW3_LIBRARY="$prefix/lib/libfftw3f.${suffix}" \
                    -DFFTW3d_INCLUDE_PATH="$prefix/include" \
                    -DFFTW3d_LIBRARY="$prefix/lib/libfftw3.${suffix}" \
                    -DFREETYPE_INCLUDE_PATH="$prefix/include/freetype2" \
                    -DFREETYPE_LIBRARY="$prefix/lib/libfreetype.${suffix}" \
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
                    -DPNG_INCLUDE_PATH="$prefix/include" \
                    -DPNG_LIBRARY="$prefix/lib/libpng.${suffix}" \
                    -DPYTHON_INCLUDE_PATH="$prefix/include/python2.7" \
                    -DPYTHON_LIBRARY="$prefix/lib/libpython2.7.${suffix}" \
                    -DTIFF_INCLUDE_PATH="$prefix/include" \
                    -DTIFF_LIBRARY="$prefix/lib/libtiff.${suffix}" \
                    -DZLIB_LIBRARY="$prefix/lib/libz.${suffix}"

make
make install

export PYTHONPATH="$HOME/EMAN2/lib:$PYTHONPATH"

source $src_dir/ci_support/post_build.sh
