#!/usr/bin/env bash

# Run tests
export PATH="$PREFIX/bin:$PATH"

e2version.py
e2speedtest.py
e2display.py -h
test -f ${PREFIX}/examples/mpi_test.py
