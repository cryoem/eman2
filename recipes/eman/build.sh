#!/bin/bash

set -xe

unset MACOSX_DEPLOYMENT_TARGET
unset CC
unset CFLAGS
unset CONDA_BACKUP_HOST
unset CPP
unset CPPFLAGS
unset CXX
unset CXXFLAGS
unset DEBUG_CFLAGS
unset DEBUG_CPPFLAGS
unset DEBUG_CXXFLAGS
unset GCC
unset GCC_AR
unset GCC_NM
unset GCC_RANLIB
unset GXX
unset HOST
unset LDFLAGS
unset _PYTHON_SYSCONFIGDATA_NAME

build_dir="${SRC_DIR}/../build_eman"

rm -rf $build_dir
mkdir -p $build_dir
cd $build_dir

cmake $SRC_DIR

make -j${CPU_COUNT}
make install
make test-verbose
