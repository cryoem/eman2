CMAKE_MINIMUM_REQUIRED(VERSION 3.14)
PROJECT(EMAN
		VERSION 2.99.70
		LANGUAGES CXX
		)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_COMPILER_IS_GNUCXX AND "${CMAKE_CXX_COMPILER_VERSION}" LESS 5)
	message(FATAL_ERROR "\nExpected GCC compiler minimum version: 5\nFound GCC version: ${CMAKE_CXX_COMPILER_VERSION}")
endif()

include(cmake/functions.cmake)

OPTION(DEBUG_CMAKE_TARGETS "enable debug output for cmake target properties" OFF)
if(DEBUG_CMAKE_TARGETS)
	set(CMAKE_DEBUG_TARGET_PROPERTIES
		INCLUDE_DIRECTORIES
		COMPILE_DEFINITIONS
		COMPILE_OPTIONS
		COMPILE_FEATURES
		LINK_DIRECTORIES
		LINK_OPTIONS
		SOURCES
		POSITION_INDEPENDENT_CODE
		CONTAINER_SIZE_REQUIRED
		LIB_VERSION
		)
endif()

# for 64-bit large file support
ADD_DEFINITIONS(-D_LARGEFILE_SOURCE)
ADD_DEFINITIONS(-D_FILE_OFFSET_BITS=64)
ADD_DEFINITIONS(-D_LARGE_FILES)

INCLUDE(${CMAKE_ROOT}/Modules/CheckIncludeFiles.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckTypeSize.cmake)

EMAN_CHECK_FUNCTION(fseek64 HAVE_FSEEK64)
EMAN_CHECK_FUNCTION(fseeko HAVE_FSEEKO)
EMAN_CHECK_FUNCTION(ftell64 HAVE_FTELL64)
EMAN_CHECK_FUNCTION(ftello HAVE_FTELLO)
# end for 64-bit large file support

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

MARK_AS_ADVANCED(CLEAR CMAKE_VERBOSE_MAKEFILE)

# Set EMAN_PREFIX
if("$ENV{CONDA_BUILD_STATE}" STREQUAL "BUILD" )
	message("ENV{CONDA_BUILD_STATE}: $ENV{CONDA_BUILD_STATE}")
	if(NOT WIN32)
		set(EMAN_PREFIX $ENV{PREFIX})
	else()
		set(EMAN_PREFIX $ENV{LIBRARY_PREFIX})
	endif()
else()
	find_package(Conda REQUIRED)
		
	if(CONDA_PREFIX)
		if(NOT WIN32)
			set(EMAN_PREFIX ${CONDA_PREFIX})
		else()
			set(EMAN_PREFIX ${CONDA_PREFIX}/Library)
		endif()
	endif()
endif()

set(EMAN_PREFIX_INC      ${EMAN_PREFIX}/include)
set(EMAN_PREFIX_LIB      ${EMAN_PREFIX}/lib)
set(CMAKE_INSTALL_PREFIX ${EMAN_PREFIX} CACHE PATH "installation prefix" FORCE)
set(CMAKE_PREFIX_PATH    ${EMAN_PREFIX} ${EMAN_PREFIX}/..)
cmake_print_variables(CMAKE_PREFIX_PATH)

include(cmake/find_all.cmake)

if("$ENV{CONDA_BUILD_STATE}" STREQUAL "BUILD" )
	set(CMAKE_INSTALL_RPATH "${Python3_SITELIB}")
else()
	set(CMAKE_INSTALL_RPATH "${Python3_SITELIB};${EMAN_PREFIX_LIB}")
endif()

option(ENABLE_CATCH2 "Enable Catch 2 framework for unit testing (internal development only)" OFF)

OPTION(ENABLE_AUTODOC "enable latex/doxygen document generation and etc. " OFF)

OPTION(ENABLE_FFTW3 "enable fftw 3 support (USE THIS)" ON)
OPTION(ENABLE_FFTW_PLAN_CACHING "enable fftw caching" ON)
OPTION(ENABLE_NATIVE_FFT "enable native fft support (for non-GPL use)" OFF)
OPTION(ENABLE_ACML_FFT "enable AMD Core Math Library fft support" OFF)
OPTION(ENABLE_FFT_CACHING "enable FFT Caching" OFF)
MARK_AS_ADVANCED(ENABLE_ACML_FFT)

OPTION(ENABLE_V4L2 "enable V4L2 framegrabbing support" OFF)

OPTION(ENABLE_SPARX_CUDA "enable SPARX CUDA usage" OFF)
OPTION(ENABLE_EMAN_CUDA "enable EMAN2 CUDA usage" OFF)

#flags used for debugging and testing
OPTION(ENABLE_MEMDEBUG "enable memory debug support of objects" OFF)
OPTION(ENABLE_MEMDEBUG2 "enable memory debug support of EMData" OFF)
OPTION(ENABLE_BOUNDS_CHECKING "enable bounds checking support" OFF)
OPTION(ENABLE_RT "enable RT support" ON)
OPTION(ENABLE_WARNINGS "display warnings during compilation" OFF)
OPTION(ENABLE_IOCACHE "enable ImageIO caching" OFF)

#flags for optimization level. You can only turn one of following option to ON, or leave all of them to OFF.
OPTION(ENABLE_DEBUG "enable debug support" OFF)
OPTION(ENABLE_PROFILE "enable symbol generation without debug support" OFF)
OPTION(ENABLE_OPTIMIZE_X86_64 "GCC optimizations for generic X86_64" OFF)
OPTION(ENABLE_OPTIMIZE_MACHINE "Optimize code for this specific computer (typ 5-10% speedup). Binaries may incompatible with some other computers." OFF)
OPTION(ENABLE_OPTIMIZE_PHI "Optimize code for the Intel compiler on the PHI" OFF)
OPTION(ENABLE_OPTIMIZE_WINDOWS_VC "Windows visual c++ specific optimization" OFF)

set(optimization_options
		ENABLE_DEBUG
		ENABLE_PROFILE
		ENABLE_OPTIMIZE_X86_64
		ENABLE_OPTIMIZE_MACHINE
		ENABLE_OPTIMIZE_PHI
		ENABLE_OPTIMIZE_WINDOWS_VC
	)

check_if_options_on_greater_than("${optimization_options}" 1)

IF(CMAKE_SYSTEM_NAME MATCHES "Darwin")
	OPTION(ENABLE_OPTIMIZE_OSX "enable OSX support" ON)
ELSE()
	OPTION(ENABLE_OPTIMIZE_OSX "enable OSX support" OFF)
ENDIF()

if(NOT WIN32)
	set(WARNING_FLAGS_ON  "-W -Wall")
	set(WARNING_FLAGS_OFF "-w")
else()
	set(WARNING_FLAGS_ON  "")
	set(WARNING_FLAGS_OFF "/W0")
endif()

if(ENABLE_WARNINGS)
	add_definitions(${WARNING_FLAGS_ON})
else()
	add_definitions(${WARNING_FLAGS_OFF})
endif()

IF(ENABLE_MEMDEBUG)
	ADD_DEFINITIONS(-DMEMDEBUG)
ENDIF()

IF(ENABLE_MEMDEBUG2)
	ADD_DEFINITIONS(-DMEMDEBUG2)
ENDIF()

IF(ENABLE_V4L2)
	ADD_DEFINITIONS(-DBOUNDS_CHECKING -DENABLE_V4L2)
ENDIF()

IF(ENABLE_BOUNDS_CHECKING)
	ADD_DEFINITIONS(-DBOUNDS_CHECKING)
ENDIF()

IF(ENABLE_DEBUG)
	SET (OPT_FLAGS "-DDEBUG -g -fno-inline")
ELSEIF(ENABLE_PROFILE)
	SET (OPT_FLAGS "-O2 -march=native  -ffast-math -pipe")
	OPTION(ENABLE_DEBUG "enable debug support" OFF)
ELSEIF(ENABLE_OPTIMIZE_OSX)
#	SET (OPT_FLAGS "-O3 -ffast-math")
	SET (OPT_FLAGS "-O3 -D _LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION")
#	SET (OPT_FLAGS "-Os -ffast-math")
ELSEIF(ENABLE_OPTIMIZE_MACHINE)
	SET (OPT_FLAGS "-O5 -march=native  -ffast-math -pipe")
ELSEIF(ENABLE_OPTIMIZE_COMPATIBILITY)
	SET (OPT_FLAGS "-O5 -march=haswell  -ffast-math -pipe")
ELSEIF(ENABLE_OPTIMIZE_X86_64)
	SET (OPT_FLAGS "-O5 -march=k8 -mfpmath=sse  -ffast-math -pipe")
#	SET (OPT_FLAGS "-O5 -mfpmath=both -march=core-avx-i  -ffast-math -pipe")
ELSEIF(ENABLE_OPTIMIZE_WINDOWS_VC)
	SET (OPT_FLAGS "/Ox /DHAVE_INLINE /MP /Gm-")
ELSEIF(ENABLE_OPTIMIZE_PHI)
	SET (OPT_FLAGS "-fast")
#	SET (OPT_FLAGS "-O3 -xHost")
ELSE()
	IF(WIN32)
		SET(OPT_FLAGS "")
	ELSE()
		SET (OPT_FLAGS "-O3 -ffast-math")     # This is the default
	ENDIF()
ENDIF()

IF(CMAKE_SYSTEM MATCHES "IRIX.*")
#	SET(EMAN_CXX_FLAGS "-n32 -LANG:std -w2 -fullwarn -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
	SET(EMAN_CXX_FLAGS "-n32 -woffall -LANG:std -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
ENDIF()

IF(CMAKE_SYSTEM MATCHES "AIX*")
	SET(EMAN_CXX_FLAGS "" CACHE INTERNAL "C++ compiler flags")
ENDIF()

IF(CMAKE_SYSTEM MATCHES "OSF1-V*")
	SET(EMAN_CXX_FLAGS "-w0" CACHE INTERNAL "C++ compiler flags")
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
#	SET(EMAN_CXX_FLAGS " " CACHE INTERNAL "EMAN CXX FLAGS")
	IF(ENABLE_DEBUG)
		SET(EMAN_CXX_FLAGS "-g -Wcast-qual -Wwrite-strings" CACHE INTERNAL "EMAN CXX FLAGS")
		SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g")
		SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g")
	ELSE()
		SET(EMAN_CXX_FLAGS "" CACHE INTERNAL "EMAN CXX FLAGS")
	ENDIF()

	if(${Boost_VERSION} LESS_EQUAL "106900")
		message("Setting _GLIBCXX_USE_CXX11_ABI=0...")
		SET(EMAN_CXX_FLAGS "${EMAN_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0" CACHE INTERNAL "EMAN CXX FLAGS")
	endif()
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EMAN_CXX_FLAGS} ${OPT_FLAGS}")

if(WIN32)
	add_definitions(-DNOMINMAX)
endif()

IF(ENABLE_IOCACHE)
	ADD_DEFINITIONS(-DIMAGEIO_CACHE)
ENDIF()


IF(ENABLE_FFTW_PLAN_CACHING)
	ADD_DEFINITIONS(-DFFTW_PLAN_CACHING)
ENDIF()

IF(ENABLE_FFT_CACHING)
	ADD_DEFINITIONS(-DFFT_CACHING)
ENDIF()

set(fft_options
	ENABLE_FFTW3
	ENABLE_NATIVE_FFT
	ENABLE_ACML_FFT
	)
check_if_options_not_equal_to("${fft_options}" 1)

IF(ENABLE_NATIVE_FFT)
	ADD_DEFINITIONS(-DNATIVE_FFT)
ENDIF()

IF(ENABLE_SPARX_CUDA)
	INCLUDE(${CMAKE_SOURCE_DIR}/libEM/cuda/FindCUDA.cmake)
	CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/libEM)
	CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/libEM/sparx/cuda)
	ADD_DEFINITIONS(-DSPARX_USING_CUDA)
	INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})
	ADD_SUBDIRECTORY(libEM/sparx/cuda)
ENDIF()

IF(ENABLE_EMAN_CUDA)
	INCLUDE(${CMAKE_SOURCE_DIR}/libEM/cuda/FindCUDA.cmake)
	CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/libEM)
	ADD_DEFINITIONS(-DEMAN2_USING_CUDA)
	INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})
	ADD_SUBDIRECTORY(libEM/cuda)
ENDIF()

ADD_SUBDIRECTORY(libEM)
ADD_SUBDIRECTORY(libpyEM)
ADD_SUBDIRECTORY(sparx)
ADD_SUBDIRECTORY(sphire)
ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(programs)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(images)

add_custom_target(PythonFiles
				  COMMAND ${CMAKE_COMMAND} -DCOMPONENT=PythonFiles -P cmake_install.cmake
				  )

if(NOT WIN32)
	install(CODE "execute_process(
					COMMAND ${CMAKE_COMMAND} -E create_symlink ./sp_gui.py ${EMAN_PREFIX}/bin/sphire
					COMMAND ${CMAKE_COMMAND} -E create_symlink ./sp_sx.py ${EMAN_PREFIX}/bin/sp_sparx
					COMMAND ${CMAKE_COMMAND} -E create_symlink ./sx.py    ${EMAN_PREFIX}/bin/sparx
				 )"
			)
endif()


IF(ENABLE_RT)
	include(cmake/testing.cmake)
	ADD_SUBDIRECTORY(rt)
ENDIF()

install(CODE "MESSAGE(\"\nTo install SPHIRE, run 'make install-sphire'.\n\")")
