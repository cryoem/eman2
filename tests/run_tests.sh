#!/usr/bin/env bash

set -xe

e2version.py
e2speedtest.py

python tests/test_EMAN2DIR.py

if [ $(whoami) != "root" ];then
    mpirun -n 4 $(which python) ${PREFIX}/examples/mpi_test.py
fi
bash ${SRC_DIR}/tests/run_prog_tests.sh
