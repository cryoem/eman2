# Centos 6.3 - x86_64
export TARGET="x86_64-redhat-linux"
export BROOT="/build"
export PREFIX="${BROOT}/local"

# PATHS
export PATH=${PREFIX}/bin:${PATH}
export PKG_CONFIG_PATH=${PREFIX}/lib/pkgconfig:${PREFIX}/python/lib/pkgconfig:${PKG_CONFIG_PATH}
export CMAKE_PREFIX_PATH=${PREFIX}
export PYTHONPATH=${PREFIX}/site-packages:${PYTHONPATH}

# Configure and compile flags
export CFLAGS="-O2 -g -I${PREFIX}/include -I${PREFIX}/python/include/python2.7"
export CXXFLAGS=$CFLAGS
export LDFLAGS="-L${PREFIX}/lib -L${PREFIX}/python/lib"

# local Python
export PATH=${PREFIX}/python/bin:${PATH}
export LD_LIBRARY_PATH=${PREFIX}/python/lib:${LD_LIBRARY_PATH}
