#!/bin/bash

curl -v -L -O http://ncmi.bcm.edu/ncmi/software/counter_222/software_121/pydusa-1.15es-fftmpi-6__2016_09_07.tgz
tar xzvf pydusa-1.15es-fftmpi-6__2016_09_07.tgz
cd pydusa-1.15es-fftmpi-6
export EMAN2DIR=$SP_DIR

make clean
make distclean

patch Makefile.am "${RECIPE_DIR}"/patch1.diff
patch Makefile.in "${RECIPE_DIR}"/patch2.diff

python "${RECIPE_DIR}"/install_mpi.py
