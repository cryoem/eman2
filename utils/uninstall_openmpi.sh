#!/usr/bin/env bash

set -x

source activate

conda remove openmpi pydusa fftw-mpi --force
