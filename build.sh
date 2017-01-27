#!/bin/bash

FFTW=fftw-3.3.6-pl1
BOOST=boost_1_63_0
FREETYPE=freetype-2.7
JPEG=jpeg-8d
SIP=sip-4.19
CMAKE=cmake-3.7.2
FTGL=ftgl-2.1.3~rc5
PNG=libpng-1.6.28
TIFF=tiff-3.8.2
DB=db-5.3.28
PYQT=PyQt4_gpl_x11-4.12
PYTHON=Python-2.7.8
HDF=hdf5-1.10.0-patch1
QT=qt-everywhere-opensource-src-4.8.7
GSL=gsl-2.3
SSL=openssl-1.0.2j

CURRENT_PATH=${PATH}
CURRENT_PYTHONPATH=${PYTHONPATH}
CURRENT_LDFLAGS=${LDFLAGS}
CURRENT_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
CURRENT_CPPFLAGS=${CPPFLAGS}
CURRENT_CFLAGS=${CFLAGS}

BDIR=${HOME}/eman2src/EMAN2 # Build "into" directory
EXT=${BDIR}/extlib # EMAN2 binary so called extlib

SDIR=${HOME}/eman2src/src # Source directory

THREADS=8

BENV=${HOME}/eman2src/src # Build environment

# Check that CMAKE is installed
command -v foo >/dev/null 2>&1 || { echo >&2 "This script requires CMake. Install and try again."; exit 1; }

# Move to build directory
cd ${SDIR}

# BUILD/INSTALL FFTW
cd ${FFTW}
./configure --enable-static=no --enable-shared=yes --prefix=${EXT}
make -j${THREADS} install
./configure --enable-static=no --enable-shared=yes --enable-float --prefix=${EXT}
make -j${THREADS} install
cd ${SDIR}

# BUILD/INSTALL GSL
cd ${GSL}
./configure --prefix=${EXT} --disable-static --enable-shared
make -j${THREADS} install
cd ${SDIR}

# BUILD/INSTALL BOOST
cd ${BOOST}
./bootstrap.sh --prefix=${EXT} --with-libraries=python,system,filesystem,thread
./b2 install --prefix=${EXT}
cd ${SDIR}

# BUILD/INSTALL FREETYPE
cd ${FREETYPE}
./configure --prefix=${EXT} --enable-shared
make -j${TREADS} install
cd ${SDIR}

# BUILD/INSTALL FTGL
export LDFLAGS="-L${EXT}/lib -L${EXT}/lib/open -lGLU -lGL -lglut -lm"
./configure --prefix=${EXT} --enable-shared
make -j${THREADS} install

# BUILD/INSTALL PNG
cd ${PNG}
./configure --prefix=${EXT}
make -j${TREADS} install
cd ${SDIR}

# BUILD/INSTALL TIFF
cd ${TIFF}
./configure --prefix=${EXT}
make -j${TREADS} install
cd ${SDIR}

# BUILD/INSTALL JPEG
cd ${JPEG}
./configure --prefix=${EXT}
make -j${TREADS} install
cd ${SDIR}

# BUILD/INSTALL QT4
cd ${QT}
./configure --prefix=${EXT}
make -j${TREADS} install
cd ${SDIR}

# BUILD/INSTALL OPENSSL
cd ${SSL}
./config --prefix=${EXT}/ssl shared
make -j${THREADS} install
cd ${SDIR}

# BUILD/INSTALL PYTHON
cd ${PYTHON}

# we use a modified version of the Modules/Setup.dist code, following these instructions:
# http://stackoverflow.com/questions/5937337/building-python-with-ssl-support-in-non-standard-location

export LDFLAGS="-L${EXT}/lib -L${EXT}/lib/open"
export LD_LIBRARY_PATH="${EXT}/lib"
export CPPFLAGS="-I${EXT}/include -I${EXT}/ssl"
export CFLAGS="-I${EXT}/include -I${EXT}/ssl"
./configure --enable-shared --prefix ${EXT} --enable-unicode=ucs4
make -j${THREADS} install

export PATH=${EXT}/bin:${PATH}
export PYTHONPATH=${EXT}/lib/python2.7/site-packages

cd ${SDIR}

# install pip into EMAN2 python environment
python get-pip.py
# fix https related issues
pip install requests[security] --upgrade
pip install pyopenssl ndg-httpsclient pyasn1

# install core EMAN2 python dependencies
pip install ipython
pip install pyopengl pyopengl-accelerate
pip install readline
pip install numpy 
pip install matplotlib

# Install SIP
cd ${SIP}
python configure.py
make install
cd ${SDIR}

# Install PyQt4
cd ${PYQT}
python ./configure.py --confirm-license -e QtCore -e QtGui -e QtOpenGL
make install
cd ${SDIR}

# install accessory EMAN2 python dependencies
pip install scipy 
pip install theano

git clone git@github.com:cryoem/eman2.git

cd ..
mkdir build
cd build

ccmake ../src/eman2
make -j${THREADS}
make install

# need to copy the .config/matplotlibrc for new matplotlib to work properly
# need to copy the eman2-installer script

# need to replace install/local rpaths 
# can probably use build script 

tar -zcvf eman2_2.2_linux64.tar.gz EMAN2

export PATH=${CURRENT_PATH}
export PYTHONPATH=${CURRENT_PYTHONPATH}
export LDFLAGS=${CURRENT_LDFLAGS}
export LD_LIBRARY_PATH=${CURRENT_LD_LIBRARY_PATH}
export CPPFLAGS=$CURRENT_CPPFLAGS}
export CFLAGS=${CURRENT_CFLAGS}
