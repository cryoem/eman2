#!/usr/bin/env bash

set -xe

MYDIR="$(cd "$(dirname "$0")"; pwd -P)"

# 1. Run e2version.py and e2speedtest.py
e2version.py
e2speedtest.py

# 2. Test git-hash consistency
if [ ${CONDA_BUILD:-0} -ne 1 ];then
    bash "${MYDIR}/test_git_hash.sh"
fi

# 3. Existence tests for data files like images, font files, JSON
python "${MYDIR}/test_EMAN2DIR.py"

# 4. Test openmpi
if [ $(whoami) != "root" ] && [ -z ${TRAVIS} ];then
    mpirun -n 4 $(which python) ${MYDIR}/../examples/mpi_test.py
fi

# 5. Run e2*.py -h
bash "${MYDIR}/run_prog_tests.sh"
