# Mac Snow Leopard
export TARGET="i686-apple-darwin10"
export BROOT="/build"
export PREFIX="${BROOT}/local"

# PATHS
export PATH=${PREFIX}/bin:${PATH}
export PKG_CONFIG_PATH=${PREFIX}/lib/pkgconfig:${PREFIX}/python/lib/pkgconfig:${PKG_CONFIG_PATH}
export CMAKE_PREFIX_PATH=${PREFIX}
export PYTHONPATH=${PREFIX}/site-packages:${PYTHONPATH}

# Configure and compile flags
export CFLAGS="-arch i386 -arch x86_64  -O2 -g -I${PREFIX}/include"
export CXXFLAGS=$CFLAGS
export LDFLAGS="-L${PREFIX}/lib -arch i386 -arch x86_64 -Wl,-headerpad_max_install_names"
