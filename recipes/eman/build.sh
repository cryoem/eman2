#!/bin/bash

set -xe

build_dir="${SRC_DIR}/../build_eman"

rm -rf $build_dir
mkdir -p $build_dir
cd $build_dir

LDFLAGS=${LDFLAGS/-Wl,-dead_strip_dylibs/}
LDFLAGS=${LDFLAGS/-Wl,-pie/}
CXXFLAGS=${CXXFLAGS/-std=c++17/-std=c++14}

cmake --version
cmake $SRC_DIR -DENABLE_WARNINGS=OFF

make -j${CPU_COUNT}
make install
