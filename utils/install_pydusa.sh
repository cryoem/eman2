#!/usr/bin/env bash

set -xe

source activate

RECIPE_DIR="${CONDA_PREFIX}/recipes"
conda build ${RECIPE_DIR}/fftw-mpi
conda build ${RECIPE_DIR}/pydusa
conda install pydusa --use-local --yes

# Cleanup
conda build purge
