#!/usr/bin/env bash

# Installs EMAN against user-specified NumPy version

if [ $# -ne 1 ];then
    echo
    echo -e '\033[35m'"  Usage: $(basename ${0}) [NumPy version]"'\033[0m'
    echo
    exit 1
fi

set -xe

source activate root

numpy_verison=${1//.}

conda install eman2=2.2=np${numpy_verison}py27_0 --use-local
