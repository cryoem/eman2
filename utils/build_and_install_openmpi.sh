#!/usr/bin/env bash

# Builds and installs openmpi

set -xe

source activate root

RECIPES_DIR=$(cd $(dirname $0)/../recipes && pwd -P)
conda build ${RECIPES_DIR}/openmpi -c defaults -c conda-forge
conda install openmpi --use-local --yes
