#!/usr/bin/env bash

set -xe

source activate root

RECIPES_DIR=$(cd $(dirname $0)/../recipes && pwd -P)
conda build ${RECIPES_DIR}/pydusa
conda remove fftw-mpi --force --yes
conda install pydusa --use-local --yes
conda install fftw-mpi --use-local --yes

conda inspect linkages pydusa
