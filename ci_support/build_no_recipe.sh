#!/usr/bin/env bash

set -x

if [ ! -z ${CI} ];then
    source ci_support/setup_conda.sh

    # Following Wiki instructions at
    # http://blake.bcm.edu/emanwiki/EMAN2/COMPILE_EMAN2_ANACONDA
    conda install eman-deps=9 -c cryoem -c defaults -c conda-forge --yes --quiet
fi

# Build and install eman2
export CPU_COUNT=2

export SRC_DIR=${PWD}
export PREFIX=${PWD}

rm -vf ${CONDA_PREFIX}/bin/e2*.py

conda info -a
conda list

bash ${SRC_DIR}/recipes/eman/build.sh

# Run tests
bash ${SRC_DIR}/tests/run_tests.sh
