#!/usr/bin/env bash

source ci_support/setup_conda.sh

# Following Wiki instructions at
# http://blake.bcm.edu/emanwiki/EMAN2/COMPILE_EMAN2_ANACONDA
if [ "$(uname -s)" != "Darwin" ];then
    conda install --yes --quiet eman-deps="*"="np18*" -c cryoem -c defaults -c conda-forge
else
    conda install --yes --quiet eman-deps -c cryoem -c defaults -c conda-forge
fi

# Build and install eman2
export SRC_DIR=${PWD}
bash ${SRC_DIR}/recipes/eman/build.sh

# Run tests
e2version.py
e2speedtest.py

mpirun -n 4 $(which python) ${SRC_DIR}/examples/mpi_test.py
bash ${SRC_DIR}/tests/run_prog_tests.sh
python ${SRC_DIR}/tests/test_EMAN2DIR.py
