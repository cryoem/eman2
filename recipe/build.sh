#!/bin/bash

set -xe

build_dir="${SRC_DIR}/../build_eman"

rm -rf $build_dir
mkdir -p $build_dir
cd $build_dir

cmake --version

CMAKE_ARGS=""
if [[ ${HOST} =~ .*linux.* ]]; then
    CMAKE_ARGS="${CMAKE_ARGS} -DCMAKE_TOOLCHAIN_FILE=${RECIPE_DIR}/cross-linux.cmake -DENABLE_OPTIMIZE_COMPATIBILITY=ON"
fi

cmake $SRC_DIR ${CMAKE_ARGS}

make -j${CPU_COUNT}
make install
