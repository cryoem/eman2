#!/bin/bash

set -xe

if [ "${CONDA_BUILD_STATE}" == "BUILD" ];then
    source "${RECIPE_DIR}"/unset_env_vars.sh
fi

build_dir="${SRC_DIR}/../build_eman"

rm -rf $build_dir
mkdir -p $build_dir
cd $build_dir

cmake $SRC_DIR

make -j${CPU_COUNT}
make install
make test-verbose
