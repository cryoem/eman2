
PROJECT(EMAN)

MARK_AS_ADVANCED(CLEAR CMAKE_VERBOSE_MAKEFILE)
OPTION(CMAKE_VERBOSE_MAKEFILE "if all commands will be echoed to the console during the make" ON)

# SET(EMAN_INSTALL_PREFIX $ENV{HOME}/EMAN2 CACHE PATH "installation prefix")
SET(EMAN_INSTALL_PREFIX   $ENV{EMAN2_PREFIX} CACHE PATH "installation prefix")

SET(CMAKE_INSTALL_PREFIX ${EMAN_INSTALL_PREFIX} CACHE INTERNAL "installation prefix")
MARK_AS_ADVANCED(CLEAR EMAN_INSTALL_PREFIX)

SET(EXECUTABLE_OUTPUT_PATH ${EMAN_BINARY_DIR}/bin CACHE PATH "output directory for building all executables.")
SET(LIBRARY_OUTPUT_PATH ${EMAN_BINARY_DIR}/lib CACHE PATH "output directory for building all libraries.")
SET(EMAN_EXECUTABLE_PATH "${EXECUTABLE_OUTPUT_PATH}")
SET(EMAN_LIBRARY_PATH "${LIBRARY_OUTPUT_PATH}")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

OPTION(ENABLE_TIFF "if enable TIFF support" ON)
OPTION(ENABLE_HDF "if enable HDF support" ON)
OPTION(ENABLE_PNG "if enable PNG support" ON)
OPTION(ENABLE_PYTHON "if enable boost python support" ON)


IF(CMAKE_SYSTEM MATCHES "IRIX.*")
    SET(EMAN_CXX_FLAGS "-g -n32 -w2 -fullwarn  -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
#    SET(EMAN_CXX_FLAGS "-O2 -n32 -woffall -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
ELSE(CMAKE_SYSTEM MATCHES "IRIX.*")
    SET(EMAN_CXX_FLAGS " -g -Wall -Werror" CACHE INTERNAL "C++ compiler flags")
#    SET(EMAN_CXX_FLAGS "-O2 -w" CACHE INTERNAL "C++ compiler flags")
ENDIF(CMAKE_SYSTEM MATCHES "IRIX.*")


SET(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER})
SET(CMAKE_ORIG_FLAGS ${CMAKE_CXX_FLAGS})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EMAN_CXX_FLAGS} ${OPT_FLAGS}")

MARK_AS_ADVANCED(CMAKE_C_FLAGS)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS)

SUBDIRS(libEM rt programs)


FIND_PATH(GSL_INCLUDE_PATH gsl/gsl_linalg.h
  $ENV{EMAN2_PREFIX}/include  
  $ENV{GSLDIR}/include $ENV{HOME}/include /usr/local/include
  /usr/include)  

FIND_LIBRARY(GSL_LIBRARY
  NAMES gsl
  PATHS $ENV{EMAN2_PREFIX}/lib
  $ENV{GSLDIR}/lib $ENV{HOME}/lib /usr/lib /usr/local/lib)

FIND_LIBRARY(GSL_CBLAS_LIBRARY
  NAMES gslcblas
  PATHS $ENV{EMAN2_PREFIX}/lib
  $ENV{GSLDIR}/lib $ENV{HOME}/lib /usr/lib /usr/local/lib)

IF(GSL_LIBRARY)
ELSE(GSL_LIBRARY)
  MESSAGE(SEND_ERROR "ERROR: GSL not found. GSL is required to use EMAN.")
ENDIF(GSL_LIBRARY)

INCLUDE_DIRECTORIES(. ${EMAN_SOURCE_DIR} ${EMAN_SOURCE_DIR}/libEM
  ${GSL_INCLUDE_PATH} 
  $ENV{EMAN2_PREFIX}/include
  ) 

IF (ENABLE_TIFF)
  FIND_PATH(TIFF_INCLUDE_PATH tiffio.h
    $ENV{TIFFDIR}/include $ENV{HOME}/include /usr/local/include /usr/include)

  IF(TIFF_INCLUDE_PATH)
    FIND_LIBRARY(TIFF_LIBRARY
      NAMES tiff
      PATHS $ENV{HOME}/lib $ENV{HOME}/EMAN/lib /usr/lib /usr/local/lib)
    IF(TIFF_LIBRARY)
      ADD_DEFINITIONS(-DEM_TIFF)
      INCLUDE_DIRECTORIES(${TIFF_INCLUDE_PATH})
    ELSE(TIFF_LIBRARY)    
      MESSAGE("WARNING: TIFF not found, ENABLE_TIFF is set back to OFF")
      OPTION(ENABLE_TIFF "if enable TIFF support" OFF)
    ENDIF(TIFF_LIBRARY)
  ELSE(TIFF_INCLUDE_PATH)
    OPTION(ENABLE_TIFF "if enable TIFF support" OFF)
    MESSAGE("WARNING: TIFF not found, ENABLE_TIFF is set back to OFF")
  ENDIF(TIFF_INCLUDE_PATH)
ENDIF(ENABLE_TIFF)

IF (ENABLE_HDF)
  FIND_PATH(HDF_INCLUDE_PATH hdf5.h
    $ENV{HDFDIR}/include $ENV{HOME}/include /usr/local/include /usr/include)
  IF(HDF_INCLUDE_PATH)
    FIND_LIBRARY(HDF_LIBRARY
      NAMES hdf5
      PATHS $ENV{HOME}/lib $ENV{HOME}/EMAN/lib /usr/lib /usr/local/lib)
    IF(HDF_LIBRARY)
      ADD_DEFINITIONS(-DEM_HDF5)
      INCLUDE_DIRECTORIES(${HDF_INCLUDE_PATH})
    ELSE(HDF_LIBRARY)    
      MESSAGE("WARNING: HDF5 not found, ENABLE_HDF is set back to OFF")
      OPTION(ENABLE_HDF "if enable HDF support" OFF)
    ENDIF(HDF_LIBRARY)
  ELSE(HDF_INCLUDE_PATH)
    MESSAGE("WARNING: HDF5 not found, ENABLE_HDF is set back to OFF")
    OPTION(ENABLE_HDF "if enable HDF support" OFF)	    
  ENDIF(HDF_INCLUDE_PATH)
ENDIF(ENABLE_HDF)



IF (ENABLE_PNG)
  FIND_PATH(PNG_INCLUDE_PATH png.h
    $ENV{PNGDIR}/include $ENV{HOME}/include /usr/local/include /usr/include)

  IF(PNG_INCLUDE_PATH)
    FIND_LIBRARY(PNG_LIBRARY
      NAMES png
      PATHS $ENV{PNGDIR}/lib $ENV{HOME}/lib $ENV{HOME}/EMAN/lib /usr/lib /usr/local/lib)
    IF(PNG_LIBRARY)
      ADD_DEFINITIONS(-DEM_PNG)
      INCLUDE_DIRECTORIES(${PNG_INCLUDE_PATH})
    ELSE(PNG_LIBRARY)
      MESSAGE("WARNING: PNG not found, ENABLE_PNG is set back to OFF")
      OPTION(ENABLE_PNG "if enable PNG support" OFF)
    ENDIF(PNG_LIBRARY)
  ELSE(PNG_INCLUDE_PATH)
    MESSAGE("WARNING: PNG not found, ENABLE_PNG is set back to OFF")
    OPTION(ENABLE_PNG "if enable PNG support" OFF)
  ENDIF(PNG_INCLUDE_PATH)
ENDIF(ENABLE_PNG)


IF(ENABLE_PYTHON)
  FIND_PATH(PYTHON_INCLUDE_PATH Python.h
    $ENV{HOME}/include/python2.2 /usr/local/include/python2.2 /usr/include/python2.2
    $ENV{HOME}/include/python2.3 /usr/local/include/python2.3 /usr/include/python2.3)
  FIND_PATH(BOOST_INCLUDE_PATH boost/python.hpp
    $ENV{EMAN2_BOOST_DIST}
    $ENV{HOME}/include /usr/local/include /usr/include)

  FIND_LIBRARY(BOOST_LIBRARY
    NAMES boost_python
    PATHS     $ENV{EMAN2_BOOST_DIST} 
    $ENV{EMAN2_PREFIX}/lib
    $ENV{HOME}/lib /usr/lib /usr/local/lib)

  IF(PYTHON_INCLUDE_PATH)
    IF(BOOST_INCLUDE_PATH AND BOOST_LIBRARY)
      SUBDIRS(libpyEM)
    ELSE(BOOST_INCLUDE_PATH AND BOOST_LIBRARY)
      MESSAGE("WARNING: boost_python not found, ENABLE_PYTHON is set back to OFF")
      OPTION(ENABLE_PYTHON "if enable boost python support" OFF)
    ENDIF(BOOST_INCLUDE_PATH AND BOOST_LIBRARY)
  ELSE(PYTHON_INCLUDE_PATH)
    MESSAGE("WARNING: Python.h not found, ENABLE_PYTHON is set back to OFF")
    OPTION(ENABLE_PYTHON "if enable boost python support" OFF)
  ENDIF(PYTHON_INCLUDE_PATH)
ENDIF(ENABLE_PYTHON)
