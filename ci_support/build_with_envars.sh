#!/usr/bin/env bash

source ci_support/pre_build.sh

export CONDA_BUILD_STATE=BUILD
export PREFIX=$CONDA_PREFIX
export SP_DIR=$(python -c "import site; print site.getsitepackages()[0]")

cmake $src_dir
make
make install

source $src_dir/ci_support/post_build.sh
