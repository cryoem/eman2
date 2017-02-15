#!/usr/bin/env bash

set -x

# Download and install Miniconda
export MINICONDA_URL="https://repo.continuum.io/miniconda"
export MINICONDA_FILE="Miniconda2-latest-Linux-x86_64.sh"

curl -L -O "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b

# Conda-install packages
$HOME/miniconda2/bin/conda config --set show_channel_urls true

$HOME/miniconda2/bin/conda install --yes --quiet qt=4 pyqt=4
$HOME/miniconda2/bin/conda install --yes --quiet -c conda-forge boost boost-cpp fftw
$HOME/miniconda2/bin/conda install --yes --quiet -c jmbell ftgl
$HOME/miniconda2/bin/conda install --yes --quiet bsddb freetype gsl hdf5 ipython jpeg libpng libtiff matplotlib numpy=1.11 pyopengl scikit-learn scipy theano tk cmake

# Build and install eman2
export build_dir=$HOME/build_eman
export src_dir=$HOME/eman2_src

mv -v $HOME/eman2 $src_dir

rm -rf ${build_dir}
mkdir -p ${build_dir}
cd ${build_dir}
