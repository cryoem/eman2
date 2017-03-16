#!/usr/bin/env bash

source ci_support/pre_build.sh

export CONDA_BUILD_STATE=BUILD
export PREFIX=$HOME/miniconda2/
export SP_DIR=$PREFIX/lib/python2.7/site-packages

cmake $src_dir
make
make install

source $src_dir/ci_support/post_build.sh
