#!/usr/bin/env bash

set -xe

source activate

RECIPE_DIR="${CONDA_PREFIX}/recipes"
conda build ${RECIPE_DIR}/pydusa
conda remove fftw-mpi --force --yes
conda install pydusa --use-local --yes
conda install fftw-mpi --use-local --yes

conda inspect linkages pydusa
