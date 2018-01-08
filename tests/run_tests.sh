#!/usr/bin/env bash

set -xe

e2version.py
e2speedtest.py

python tests/test_EMAN2DIR.py

mpirun -n 4 $(which python) ${SRC_DIR}/examples/mpi_test.py
bash ${SRC_DIR}/tests/run_prog_tests.sh
