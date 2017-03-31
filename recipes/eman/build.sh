#!/bin/bash

build_dir="${SRC_DIR}/../build_eman"

rm -rf $build_dir
mkdir -p $build_dir
cd $build_dir

cmake $SRC_DIR

make -j${CPU_COUNT}
make install

ln -s $PREFIX/bin/e2version.py $SP_DIR/e2version.py
