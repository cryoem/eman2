#!/usr/bin/env bash

set -xe

MYDIR="$(cd "$(dirname "$0")"; pwd -P)"

bash "${MYDIR}/../tests/future_import_tests.sh"

if [ ! -z ${TRAVIS} ];then
    source ci_support/setup_conda.sh

    # Following Wiki instructions at
    # http://blake.bcm.edu/emanwiki/EMAN2/COMPILE_EMAN2_ANACONDA
    conda install eman-deps=15.0 -c cryoem -c defaults -c conda-forge --yes --quiet
fi

if [ ! -z ${CIRCLECI} ];then
    source ${HOME}/miniconda2/bin/activate root
fi

python -m compileall -q .

# Build and install eman2
rm -vf ${CONDA_PREFIX}/bin/e2*.py

conda info -a
conda list
conda list --explicit

if [ ! -z "$JENKINS_HOME" ];then
    CPU_COUNT=4
else
    CPU_COUNT=2
fi

build_dir="../build_eman"
src_dir=${PWD}

rm -rf $build_dir
mkdir -p $build_dir
cd $build_dir

cmake --version
cmake ${src_dir} -DENABLE_WARNINGS=OFF
make -j${CPU_COUNT}
make install

# Run tests
cd "${src_dir}"
bash tests/run_tests.sh
