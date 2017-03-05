#!/usr/bin/env bash

source ci_support/run_circle_pre.sh

export CONDA_BUILD_STATE=BUILD
export PREFIX=$HOME/miniconda2/
export SP_DIR=$PREFIX/lib/python2.7/site-packages

cmake $src_dir -DCMAKE_INSTALL_RPATH="$HOME/EMAN2/lib"
make
make install

source $src_dir/ci_support/run_circle_post.sh
