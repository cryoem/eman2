#!/usr/bin/env bash

# Installs Pydusa against user-specified NumPy version

if [ $# -ne 1 ];then
    echo
    echo -e '\033[35m'"  Usage: $(basename ${0}) [NumPy version]"'\033[0m'
    echo
    exit 1
fi

set -xe

source activate root

RECIPES_DIR=$(cd $(dirname $0)/../recipes && pwd -P)
numpy_verison=${1//.}

conda build ${RECIPES_DIR}/pydusa --numpy ${1}
conda remove fftw-mpi --force --yes
conda install pydusa=1.15=np${numpy_verison}_0 --use-local --yes
conda install fftw-mpi --use-local --yes

conda inspect linkages pydusa
