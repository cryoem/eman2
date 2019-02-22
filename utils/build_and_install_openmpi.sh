#!/usr/bin/env bash

# Builds and installs openmpi

set -xe

RECIPES_DIR=$(cd $(dirname $0)/../recipes && pwd -P)
conda build ${RECIPES_DIR}/openmpi -c defaults -c conda-forge -c conda-forge/label/cf201901

CONDA_BLD=$(conda info --root)/conda-bld
if [[ ! -d ${CONDA_BLD} ]]
then
    CONDA_BLD=$(dirname ${CONDA_EXE})/../conda-bld
fi

conda install openmpi -c file://${CONDA_BLD} --override-channels --yes
