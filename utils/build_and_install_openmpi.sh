#!/usr/bin/env bash

# Builds and installs openmpi

set -xe

RECIPES_DIR=$(cd $(dirname $0)/../recipes && pwd -P)
conda build ${RECIPES_DIR}/openmpi -c defaults -c conda-forge -c conda-forge/label/cf201901

CONDA_BLD=$(dirname $(which conda))/../conda-bld
if [[ ! -d ${CONDA_BLD} ]]
then
    CONDA_BLD=$(dirname ${CONDA_EXE})/../conda-bld
fi
if [[ ! -d ${CONDA_BLD} ]]
then
    echo "Conda-bld not found"
    exit 1
fi

conda install openmpi -c file://${CONDA_BLD} --override-channels --yes
