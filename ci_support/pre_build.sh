#!/usr/bin/env bash

set -x

# Download and install Miniconda
export MINICONDA_URL="https://repo.continuum.io/miniconda"

curl -L -O "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b

# Conda-install packages
source ${HOME}/miniconda2/bin/activate root
conda config --set show_channel_urls true

conda install --yes --quiet qt=4 pyqt=4
conda install --yes --quiet -c conda-forge boost boost-cpp fftw
conda install --yes --quiet -c cryoem ftgl
conda install --yes --quiet bsddb freetype gsl hdf5 ipython jpeg libpng libtiff matplotlib numpy=1.11 pyopengl scikit-learn scipy theano tk cmake

# Build and install eman2
export build_dir=$HOME/build_eman

if [ -e ${HOME}/eman2 ];then
    export src_dir=$HOME/eman2_src
    mv -v $HOME/eman2 $src_dir
else
    export src_dir=${PWD}
fi

rm -rf ${build_dir}
mkdir -p ${build_dir}
cd ${build_dir}
