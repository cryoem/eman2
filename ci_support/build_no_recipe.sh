#!/usr/bin/env bash

set -xe

if [ ! -z ${CI} ];then
    source ci_support/setup_conda.sh

    # Following Wiki instructions at
    # http://blake.bcm.edu/emanwiki/EMAN2/COMPILE_EMAN2_ANACONDA
    conda install eman-deps=9 -c cryoem -c defaults -c conda-forge --yes --quiet
fi

# Build and install eman2
rm -vf ${CONDA_PREFIX}/bin/e2*.py

conda info -a
conda list

if [ -z "$JENKINS_HOME" ];then
    CPU_COUNT=4
else
    CPU_COUNT=2
fi

build_dir="../build_eman"
src_dir=${PWD}

rm -rf $build_dir
mkdir -p $build_dir
cd $build_dir

cmake ${src_dir}
make -j${CPU_COUNT}
make install

# Run tests
bash ${src_dir}/tests/run_tests.sh
