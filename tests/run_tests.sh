#!/usr/bin/env bash

set -xe

MYDIR="$(cd "$(dirname "$0")"; pwd -P)"

e2version.py
e2speedtest.py

if [ ${CONDA_BUILD:-0} -ne 1 ];then
    bash "${MYDIR}/test_git_hash.sh"
fi

python "${MYDIR}/test_EMAN2DIR.py"

if [ $(whoami) != "root" ];then
    mpirun -n 4 $(which python) ${MYDIR}/../examples/mpi_test.py
fi
bash "${MYDIR}/run_prog_tests.sh"
