
PROJECT(EMAN)

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

OPTION(ENABLE_TIF "if enable TIF support" ON)
OPTION(ENABLE_HDF "if enable HDF support" ON)
OPTION(ENABLE_PNG "if enable PNG support" ON)


IF(CMAKE_SYSTEM MATCHES "IRIX.*")
	SET(EMAN_CXX_FLAGS "-g -n32 -w2 -fullwarn" CACHE INTERNAL "C++ compiler flags")
	ADD_DEFINITIONS(-DSGI -DWORDS_BIGENDIAN)
ELSE(CMAKE_SYSTEM MATCHES "IRIX.*")
	SET(EMAN_CXX_FLAGS " -g -Wall -Werror" CACHE INTERNAL "C++ compiler flags")
ENDIF(CMAKE_SYSTEM MATCHES "IRIX.*")

# check if config.h exists, if not create an empty one to meet the original autoconfig requirements
IF(EXISTS ${EMAN_SOURCE_DIR}/config.h)
ELSE(EXISTS ${EMAN_SOURCE_DIR}/config.h)
	EXEC_PROGRAM(touch ARGS ${EMAN_SOURCE_DIR}/config.h)
ENDIF(EXISTS ${EMAN_SOURCE_DIR}/config.h)


SET(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER})
SET(CMAKE_ORIG_FLAGS ${CMAKE_CXX_FLAGS})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EMAN_CXX_FLAGS} ${OPT_FLAGS}")
MARK_AS_ADVANCED(CMAKE_C_FLAGS)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS)

SUBDIRS(libEM rt programs)
INCLUDE_DIRECTORIES(. ${EMAN_SOURCE_DIR} ${EMAN_SOURCE_DIR}/libEM)


# this looks for TIF support
# it will define the following values
#
# TIF_INCLUDE_PATH = where tiffio.h can be found
# TIF_LIBRARY      = the library to link against (tiff)

IF (ENABLE_TIF)
	# this looks for TIF support
	# it will define the following values
	#
	# TIF_LIBRARY 

	FIND_PATH(TIF_INCLUDE_PATH tiffio.h
          $ENV{TIFDIR}/include $ENV{HOME}/include /usr/local/include /usr/include)

	IF(TIF_INCLUDE_PATH)
		FIND_LIBRARY(TIF_LIBRARY
		NAMES tiff
		PATHS $ENV{HOME}/lib $ENV{HOME}/EMAN/lib /usr/lib /usr/local/lib)
		IF(TIF_LIBRARY)
		    ADD_DEFINITIONS(-DEM_TIFF)
		    INCLUDE_DIRECTORIES(${TIF_INCLUDE_PATH})
		ELSE(TIF_LIBRARY)    
		    MESSAGE("WARNING: TIF not found, ENABLE_TIF is set back to OFF")
		    OPTION(ENABLE_TIF "if enable TIF support" OFF)
		ENDIF(TIF_LIBRARY)
	ELSE(TIF_INCLUDE_PATH)
	    MESSAGE("WARNING: TIF not found, ENABLE_TIF is set back to OFF")
	    OPTION(ENABLE_TIF "if enable TIF support" OFF)	    
	ENDIF(TIF_INCLUDE_PATH)
ENDIF(ENABLE_TIF)

IF (ENABLE_HDF)
	# this looks for HDF support
	# it will define the following values
	#
	# HDF_LIBRARY 

	FIND_PATH(HDF_INCLUDE_PATH hdf5.h
          $ENV{HDFDIR}/include $ENV{HOME}/include /usr/local/include /usr/include)
	IF(HDF_INCLUDE_PATH)
		FIND_LIBRARY(HDF_LIBRARY
		NAMES hdf5
		PATHS $ENV{HOME}/lib $ENV{HOME}/EMAN/lib /usr/lib /usr/local/lib)
		IF(HDF_LIBRARY)
		    ADD_DEFINITIONS(-DHDF5)
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
	# this looks for PNG support
	# it will define the following values
	#
	# PNG_LIBRARY 

	FIND_PATH(PNG_INCLUDE_PATH png.h
          $ENV{PNGDIR}/include $ENV{HOME}/include /usr/local/include /usr/include)

	IF(PNG_INCLUDE_PATH)
		FIND_LIBRARY(PNG_LIBRARY
	        NAMES png
	        PATHS $ENV{PNGDIR}/lib $ENV{HOME}/lib $ENV{HOME}/EMAN/lib /usr/lib /usr/local/lib)
		IF(PNG_LIBRARY)
			ADD_DEFINITIONS(-DPNG)
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






