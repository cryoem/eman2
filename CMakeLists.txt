PROJECT(EMAN)

MARK_AS_ADVANCED(FORCE CMAKE_BACKWARDS_COMPATIBILITY)
MARK_AS_ADVANCED(CLEAR CMAKE_VERBOSE_MAKEFILE)
OPTION(CMAKE_VERBOSE_MAKEFILE "if all commands will be echoed to the console during the make" ON)

SET(EMAN_INSTALL_PREFIX $ENV{HOME}/EMAN2 CACHE PATH "installation prefix")
SET(CMAKE_INSTALL_PREFIX ${EMAN_INSTALL_PREFIX} CACHE INTERNAL "installation prefix")
MARK_AS_ADVANCED(CLEAR EMAN_INSTALL_PREFIX)

SET(EXECUTABLE_OUTPUT_PATH ${EMAN_BINARY_DIR}/bin CACHE PATH "output directory for building all executables.")
SET(LIBRARY_OUTPUT_PATH ${EMAN_BINARY_DIR}/lib CACHE PATH "output directory for building all libraries.")
SET(EMAN_EXECUTABLE_PATH "${EXECUTABLE_OUTPUT_PATH}")
SET(EMAN_LIBRARY_PATH "${LIBRARY_OUTPUT_PATH}")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

OPTION(ENABLE_FFTW2 "if enable fftw 2 support" ON)
OPTION(ENABLE_FFTW3 "if enable fftw 3 support" OFF)
OPTION(ENABLE_TIFF "if enable TIFF support" ON)
OPTION(ENABLE_HDF5 "if enable HDF5 support" ON)
OPTION(ENABLE_PNG "if enable PNG support" ON)
OPTION(ENABLE_PYTHON "if enable boost python support" ON)

OPTION(ENABLE_DEBUG "if enable debug support" ON)

OPTION(ENABLE_ATHLON "if enable ATHLON specific support" OFF)
OPTION(ENABLE_OPTERON "if enable OPTERON/Ath64 specific support" OFF)
OPTION(ENABLE_G5 "if enable G5 CPU specific support" OFF)

IF (ENABLE_DEBUG)
    SET (OPT_FLAGS "-DDEBUG")
ELSE (ENABLE_DEBUG)
    IF (ENABLE_ATHLON)
	SET (OPT_FLAGS "-O3 -march=athlon-xp -mfpmath=sse -msse -funroll-loops -pipe")
    ELSE(ENABLE_ATHLON)
	IF (ENABLE_OPTERON)
	    SET (OPT_FLAGS "-O3 -march=k8 -m64 -mfpmath=sse -ffast-math -msse2 -funroll-loops -pipe")
	ELSE (ENABLE_OPTERON)
	    SET (OPT_FLAGS "-O2")
	ENDIF(ENABLE_OPTERON)
    ENDIF(ENABLE_ATHLON)

    IF(ENABLE_G5)
	SET(OPT_FLAGS "-fast -fPIC -mcpu=970 -mtune=970 -mpowerpc64 -O5" CACHE INTERNAL "C++ compiler flags")
    ENDIF(ENABLE_G5)
    
ENDIF (ENABLE_DEBUG)


IF(CMAKE_SYSTEM MATCHES "IRIX.*")
    SET(EMAN_CXX_FLAGS "-n32 -LANG:std -w2 -fullwarn -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
    #SET(EMAN_CXX_FLAGS "-n32 -woffall -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
ENDIF(CMAKE_SYSTEM MATCHES "IRIX.*")

IF(CMAKE_SYSTEM MATCHES "AIX*")
    SET(EMAN_CXX_FLAGS "" CACHE INTERNAL "C++ compiler flags")
ENDIF(CMAKE_SYSTEM MATCHES "AIX*")

IF(CMAKE_SYSTEM MATCHES "OSF1-V*")
    SET(EMAN_CXX_FLAGS "-w0" CACHE INTERNAL "C++ compiler flags")
ENDIF(CMAKE_SYSTEM MATCHES "OSF1-V*")

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(EMAN_CXX_FLAGS "-Wall -Werror" CACHE INTERNAL "EMAN CXX FLAGS")
    #SET(EMAN_CXX_FLAGS "-Wall" CACHE INTERNAL "EMAN CXX FLAGS")
    SET(LIB64 "/usr/lib64" CACHE INTERNAL "lib64")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

SET(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER})
SET(CMAKE_ORIG_FLAGS ${CMAKE_CXX_FLAGS})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EMAN_CXX_FLAGS} ${OPT_FLAGS}")

MARK_AS_ADVANCED(CMAKE_C_FLAGS)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS)

# for 64-bit large file support
ADD_DEFINITIONS(-D_LARGEFILE_SOURCE)
ADD_DEFINITIONS(-D_FILE_OFFSET_BITS=64)
ADD_DEFINITIONS(-D_LARGE_FILES)

INCLUDE(${CMAKE_ROOT}/Modules/CheckIncludeFiles.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)

MACRO(EMAN_CHECK_FUNCTION FUNCTION VARIABLE)
    CHECK_FUNCTION_EXISTS(${FUNCTION} ${VARIABLE})
    IF(${VARIABLE})
	ADD_DEFINITIONS(-D${VARIABLE})
    ENDIF(${VARIABLE})
ENDMACRO(EMAN_CHECK_FUNCTION)

EMAN_CHECK_FUNCTION(fseek64 HAVE_FSEEK64)
EMAN_CHECK_FUNCTION(fseeko HAVE_FSEEKO)
EMAN_CHECK_FUNCTION(ftell64 HAVE_FTELL64)
EMAN_CHECK_FUNCTION(ftello HAVE_FTELLO)
# end for 64-bit large file support

SUBDIRS(libEM plugins rt programs utils)

SET(INCLUDE_DIRS ${EMAN_SOURCE_DIR}/libEM $ENV{HOME}/include $ENV{FFTWDIR}/include $ENV{GSLDIR}/include /usr/local/include /usr/include /sw/include)
SET(LIB_DIRS $ENV{EMANDIR}/lib $ENV{HOME}/lib ${LIB64} $ENV{FFTWDIR}/lib $ENV{GSLDIR}/lib /usr/local/lib /usr/lib /sw/lib)

MACRO(CHECK_REQUIRED_LIB upper lower header)
    FIND_PATH(${upper}_INCLUDE_PATH ${header} $ENV{${upper}DIR}/include ${INCLUDE_DIRS})
    FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} PATHS $ENV{${upper}DIR}/lib ${LIB_DIRS})

    IF(${upper}_LIBRARY)
    ELSE(${upper}_LIBRARY)
        MESSAGE(SEND_ERROR "ERROR: ${upper} not found. please install ${upper} first!")
    ENDIF(${upper}_LIBRARY)    
ENDMACRO(CHECK_REQUIRED_LIB)

MACRO(CHECK_OPTIONAL_LIB upper lower header)
    FIND_PATH(${upper}_INCLUDE_PATH ${header} $ENV{${upper}DIR}/include ${INCLUDE_DIRS}) 
    
    IF(${upper}_INCLUDE_PATH)
        FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} PATHS $ENV{${upper}DIR}/lib ${LIB_DIRS})
	
        IF(${upper}_LIBRARY)
	    ADD_DEFINITIONS(-DEM_${upper})
	    INCLUDE_DIRECTORIES(${${upper}_INCLUDE_PATH})
        ELSE(${upper}_LIBRARY)    
	    MESSAGE("WARNING: ${upper} not found, ENABLE_${upper} is set back to OFF")
	    OPTION(ENABLE_${upper} "if enable ${upper} support" OFF)
        ENDIF(${upper}_LIBRARY)
	
    ELSE(${upper}_INCLUDE_PATH)
	MESSAGE("WARNING: ${upper} not found, ENABLE_${upper} is set back to OFF")
        OPTION(ENABLE_${upper} "if enable ${upper} support" OFF)
    ENDIF(${upper}_INCLUDE_PATH)
ENDMACRO(CHECK_OPTIONAL_LIB)

IF(ENABLE_FFTW2)
    CHECK_REQUIRED_LIB(FFTW2 sfftw sfftw.h)
    CHECK_REQUIRED_LIB(RFFTW2 srfftw srfftw.h)
    ADD_DEFINITIONS(-DFFTW2)
    SET(FFTW_LIBRARIES ${RFFTW2_LIBRARY} ${FFTW2_LIBRARY})
    SET(FFTW_INCLUDE_PATH ${FFTW2_INCLUDE_PATH})
ENDIF(ENABLE_FFTW2)

IF(ENABLE_FFTW3)
    CHECK_REQUIRED_LIB(FFTW3 fftw3f fftw3.h)
    ADD_DEFINITIONS(-DFFTW3)
    SET(FFTW_LIBRARIES ${FFTW3_LIBRARY})
    SET(FFTW_INCLUDE_PATH ${FFTW3_INCLUDE_PATH})
ENDIF(ENABLE_FFTW3)

CHECK_REQUIRED_LIB(GSL gsl gsl/gsl_linalg.h)
CHECK_REQUIRED_LIB(GSL_CBLAS gslcblas gsl/gsl_linalg.h)

INCLUDE_DIRECTORIES(. ${EMAN_SOURCE_DIR} ${EMAN_SOURCE_DIR}/libEM ${GSL_INCLUDE_PATH}  ${FFTW_INCLUDE_PATH}) 

IF(ENABLE_TIFF)
    CHECK_OPTIONAL_LIB(TIFF tiff tiffio.h)
ENDIF(ENABLE_TIFF)

IF(ENABLE_HDF5)
    CHECK_OPTIONAL_LIB(HDF5 hdf5 hdf5.h)
ENDIF(ENABLE_HDF5)

IF(ENABLE_PNG)
    CHECK_OPTIONAL_LIB(PNG png png.h)
ENDIF(ENABLE_PNG)
    
IF(ENABLE_PYTHON)
    FIND_PATH(PYTHON_INCLUDE_PATH Python.h
	$ENV{HOME}/include/python2.2 /usr/local/include/python2.2 /usr/include/python2.2
	$ENV{HOME}/include/python2.3 /usr/local/include/python2.3 /usr/include/python2.3)

    IF(WIN32)
	FIND_LIBRARY(PYTHON_LIBRARY
	    NAMES python23 python22 python
	    PATHS C:/Python23/libs C:/Python22/libs)
    ENDIF(WIN32)
    
    IF(PYTHON_INCLUDE_PATH)
	CHECK_OPTIONAL_LIB(BOOST boost_python boost/python.hpp)
	IF(BOOST_INCLUDE_PATH AND BOOST_LIBRARY)
	    SUBDIRS(libpyEM)
	ENDIF(BOOST_INCLUDE_PATH AND BOOST_LIBRARY)
    ELSE(PYTHON_INCLUDE_PATH)
	MESSAGE("WARNING: Python.h not found, ENABLE_PYTHON is set back to OFF")
	OPTION(ENABLE_PYTHON "if enable boost python support" OFF)
    ENDIF(PYTHON_INCLUDE_PATH)
    
ENDIF(ENABLE_PYTHON)


MARK_AS_ADVANCED(FORCE CMAKE_BACKWARDS_COMPATIBILITY)
