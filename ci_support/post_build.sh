#!/usr/bin/env bash

set -e

# install pydusa
bash -e $src_dir/recipes/eman/install_pydusa.sh

# Run tests
export PATH="$PREFIX/bin:$PATH"

e2version.py
e2speedtest.py
e2display.py -h
mpirun -n 4 $(which python) ${PREFIX}/examples/mpi_test.py
