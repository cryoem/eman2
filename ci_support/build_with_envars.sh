#!/usr/bin/env bash

set -xe

# Download and install Miniconda
export MINICONDA_URL="https://repo.continuum.io/miniconda"

curl -L -O "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b

# Conda-install packages
source ${HOME}/miniconda2/bin/activate root
conda config --set show_channel_urls true

conda install --yes --quiet qt=4 pyqt=4
conda install --yes --quiet -c conda-forge boost=1.63.* boost-cpp=1.63.* fftw cmake>=3.8
conda install --yes --quiet -c cryoem ftgl
conda install --yes --quiet bsddb freetype gsl hdf5 ipython jpeg libpng libtiff matplotlib numpy=1.11 pyopengl scikit-learn scipy theano tk
conda install --yes --quiet -c cryoem pydusa

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

cmake $src_dir -DENABLE_CONDA=ON
make
make install
make test-verbose

set -e

# Run tests
e2version.py
e2speedtest.py
mpirun -n 4 $(which python) ${src_dir}/examples/mpi_test.py
cd ${src_dir}
bash tests/run_prog_tests.sh
python tests/test_EMAN2DIR.py
