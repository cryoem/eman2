#!/usr/bin/env bash

# Installs Pydusa against user-specified NumPy version

if [ $# -ne 1 ];then
    echo
    echo -e '\033[35m'"  Usage: $(basename ${0}) [NumPy version]"'\033[0m'
    echo
    exit 1
fi

set -xe

RECIPES_DIR=$(cd $(dirname $0)/../recipes && pwd -P)
numpy_version=${1//.}

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

#echo $(conda remove fftw-mpi --force --yes)
conda install fftw-mpi --yes --force-reinstall --override-channels -c file://${CONDA_BLD}
conda install pydusa=1.15=np${numpy_version}* --force-reinstall --yes -c file://${CONDA_BLD} --override-channels

conda inspect linkages pydusa
