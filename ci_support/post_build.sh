#!/usr/bin/env bash

set -e

# Run tests
export PATH="$PREFIX/bin:$PATH"
ln -s $PREFIX/bin/e2version.py $SP_DIR/e2version.py

e2version.py
e2speedtest.py
e2display.py -h
mpirun -n 4 $(which python) ${PREFIX}/examples/mpi_test.py
cd ${src_dir}
bash tests/run_prog_tests.sh
