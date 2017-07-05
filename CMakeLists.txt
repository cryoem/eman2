CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(EMAN)

IF(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
       cmake_policy(SET CMP0002 OLD)
ENDIF(COMMAND cmake_policy)

FUNCTION(message_var var)
	message("${var}: ${${var}}")
endfunction()

MARK_AS_ADVANCED(CLEAR CMAKE_VERBOSE_MAKEFILE)
OPTION(CMAKE_VERBOSE_MAKEFILE "if all commands will be echoed to the console during the make" ON)

SET(CMAKE_MACOSX_RPATH ON)

OPTION(ENABLE_CONDA "enable conda-supported build" ON)

MACRO(CHECK_PYTHON)
	find_package(PythonInterp REQUIRED)
	find_package(PythonLibs   REQUIRED)
	
	set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_DIRS} CACHE PATH "")
	set(PYTHON_LIBRARY      ${PYTHON_LIBRARIES}    CACHE PATH "")
	
	message("PYTHON_EXECUTABLE:   ${PYTHON_EXECUTABLE}")
	message("PYTHON_LIBRARIES:    ${PYTHON_LIBRARIES}")
	message("PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
	message("PYTHON_INCLUDE_PATH: ${PYTHON_INCLUDE_PATH}")
	message("PYTHON_INCLUDE_DIR:  ${PYTHON_INCLUDE_DIR}")
ENDMACRO(CHECK_PYTHON)

# Handle build strategies, 
#        set CONDA_PREFIX and SP_DIR
#        set EMAN_PY_INSTALL_PREFIX and EMAN_INSTALL_PREFIX
#        write eman2dir relative to SP_DIR path
# If it is conda-build build, ENABLE_CONDA is ignored
# If it is not conda-build build, ENABLE_CONDA is checked
# If it is not conda-build build and ENABLE_CONDA is not set, it is legacy build
if(ENABLE_CONDA)
	if("$ENV{CONDA_BUILD_STATE}" STREQUAL "BUILD" )
		if(NOT WIN32)
			set(EMAN_PREFIX $ENV{PREFIX})
		else()
			set(EMAN_PREFIX $ENV{LIBRARY_PREFIX})
		endif()
		
		set(SP_DIR $ENV{SP_DIR})
	elseif(DEFINED ENV{CONDA_PREFIX})
		set(EMAN_PREFIX $ENV{CONDA_PREFIX} CACHE PATH "")
	else()
		# Assuming the active conda environment is on PATH, this finds the path of bin/ in the environment
		find_program(CONDA_EXECUTABLE conda
					 PATHS ${CONDA_ROOT}/bin ${CONDA_ROOT}/Scripts ENV PATH
				 )
		if(CONDA_EXECUTABLE)
			execute_process(COMMAND ${CONDA_EXECUTABLE} info --root
					OUTPUT_VARIABLE out_var
					OUTPUT_STRIP_TRAILING_WHITESPACE
					)
			if(out_var)
				set(CONDA_ROOT ${out_var}            CACHE PATH "")
			else()
				set(CONDA_ROOT CONDA_ROOT-NOTFOUND   CACHE PATH "")
			endif()
		endif()
			
		if(CONDA_ROOT)
			if(WIN32)
				set(EMAN_PREFIX ${CONDA_ROOT}/Library CACHE PATH "")
			else()
				set(EMAN_PREFIX ${CONDA_ROOT} CACHE PATH "")
			endif()
		endif()
	endif()
endif()

CHECK_PYTHON()
if(NOT WIN32)
	set(py_sp_dir_command "import site; print site.getsitepackages()[0]")
else()
	set(py_sp_dir_command "import site; print site.getsitepackages()[1]")
endif()
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "${py_sp_dir_command}"
		OUTPUT_VARIABLE SP_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE
		)
message("Python site-packages: ${SP_DIR}")
		
	SET(EMAN_INSTALL_PREFIX    ${EMAN_PREFIX} CACHE PATH "installation prefix")
	set(EMAN_PY_INSTALL_PREFIX ${SP_DIR})
		
	if(NOT WIN32)
		set(eman2dir_relative_path "../../../")
	else()
		set(eman2dir_relative_path "../../Library/")
	endif()

if(NOT ENABLE_CONDA)
	SET(EMAN_INSTALL_PREFIX    $ENV{HOME}/EMAN2 CACHE PATH "installation prefix")
	set(EMAN_PY_INSTALL_PREFIX ${EMAN_INSTALL_PREFIX}/lib)
	set(eman2dir_relative_path "../")
endif()

if(EMAN_PY_INSTALL_PREFIX)
	file(WRITE ${EMAN_PY_INSTALL_PREFIX}/eman2dir_relative_path_to_sp_dir ${eman2dir_relative_path})
endif()

SET(CMAKE_INSTALL_PREFIX ${EMAN_INSTALL_PREFIX} CACHE INTERNAL "installation prefix")
MARK_AS_ADVANCED(CLEAR EMAN_INSTALL_PREFIX)

set(CMAKE_INSTALL_RPATH ${EMAN_PY_INSTALL_PREFIX})

OPTION(ENABLE_AUTODOC "enable latex/doxygen document generation and etc. " OFF)

OPTION(ENABLE_FFTW2 "enable fftw 2 support (LEGACY)" OFF)
OPTION(ENABLE_FFTW3 "enable fftw 3 support (USE THIS)" ON)
OPTION(ENABLE_FFTW_PLAN_CACHING "enable fftw caching" ON)
OPTION(ENABLE_NATIVE_FFT "enable native fft support (for non-GPL use)" OFF)
OPTION(ENABLE_ACML_FFT "enable AMD Core Math Library fft support" OFF)
OPTION(ENABLE_FFT_CACHING "enable FFT Caching" OFF)
MARK_AS_ADVANCED(ENABLE_ACML_FFT)
OPTION(ENABLE_DJBFFT "enable djbfft support" OFF)
MARK_AS_ADVANCED(ENABLE_DJBFFT)

OPTION(ENABLE_OPTPP "enable Opt++ support" OFF)
OPTION(ENABLE_TAO "enable TAO support" OFF)
MARK_AS_ADVANCED(ENABLE_TAO)

OPTION(ENABLE_TIFF "enable TIFF support" ON)
OPTION(ENABLE_V4L2 "enable V4L2 framegrabbing support" OFF)
OPTION(ENABLE_PNG "enable PNG support" ON)
OPTION(ENABLE_JPEG "enable JPEG support" ON)

OPTION(ENABLE_OPENGL "enable EMAN2 OpenGL usage" ON)
OPTION(ENABLE_SPARX_CUDA "enable SPARX CUDA usage" OFF)
OPTION(ENABLE_EMAN_CUDA "enable EMAN2 CUDA usage" OFF)
OPTION(ENABLE_FTGL "enable EMAN2 FTGL (3D fonts) usage" ON)
OPTION(ENABLE_STATIC_FTGL "If build EMAN with libftgl.a" OFF)
MARK_AS_ADVANCED(ENABLE_STATIC_FTGL)

#flags used for debugging and testing
OPTION(ENABLE_MEMDEBUG "enable memory debug support of objects" OFF)
OPTION(ENABLE_MEMDEBUG2 "enable memory debug support of EMData" OFF)
OPTION(ENABLE_BOUNDS_CHECKING "enable bounds checking support" OFF)
OPTION(ENABLE_RT "enable RT support" ON)

OPTION(ENABLE_IOCACHE "enable ImageIO caching" OFF)
# OPTION(ENABLE_HDFIOCACHE "enable HDF5 file caching" OFF)

#flags for optimization level. You can only turn one of following option to ON, or leave all of them to OFF.
OPTION(ENABLE_DEBUG "enable debug support" OFF)
OPTION(ENABLE_PROFILE "enable symbol generation without debug support" OFF)
OPTION(ENABLE_OPTIMIZE_X86_64 "GCC optimizations for generic X86_64" OFF)
OPTION(ENABLE_OPTIMIZE_MACHINE "Optimize code for this specific computer (typ 5-10% speedup). Binaries may incompatible with some other computers." OFF)
OPTION(ENABLE_OPTIMIZE_PHI "Optimize code for the Intel compiler on the PHI" OFF)
OPTION(ENABLE_OPTIMIZE_WINDOWS_VC "Windows visual c++ specific optimization" OFF)

IF(CMAKE_SYSTEM_NAME MATCHES "Darwin")
	OPTION(ENABLE_OPTIMIZE_OSX "enable OSX support" ON)
ELSE()
	OPTION(ENABLE_OPTIMIZE_OSX "enable OSX support" OFF)
ENDIF()

IF (ENABLE_MEMDEBUG)
	ADD_DEFINITIONS(-DMEMDEBUG -W -Wall)
ENDIF (ENABLE_MEMDEBUG)

IF (ENABLE_MEMDEBUG2)
	ADD_DEFINITIONS(-DMEMDEBUG2 -W -Wall)
ENDIF (ENABLE_MEMDEBUG2)

IF (ENABLE_V4L2)
	ADD_DEFINITIONS(-DBOUNDS_CHECKING -DENABLE_V4L2)
ENDIF (ENABLE_V4L2)

IF (ENABLE_BOUNDS_CHECKING)
	ADD_DEFINITIONS(-DBOUNDS_CHECKING)
ENDIF (ENABLE_BOUNDS_CHECKING)

IF (ENABLE_DEBUG)
    SET (OPT_FLAGS "-DDEBUG -O2 -g -ansi -Wall -W -fno-inline")
ELSEIF(ENABLE_PROFILE)
    SET (OPT_FLAGS "-O2 -ansi -march=native -m64 -ffast-math -pipe")
    OPTION(ENABLE_DEBUG "enable debug support" OFF)
ELSEIF (ENABLE_OPTIMIZE_OSX)
#    SET (OPT_FLAGS "-O3 -ffast-math")
    SET (OPT_FLAGS "-O3")
#    SET (OPT_FLAGS "-Os -ffast-math")
ELSEIF (ENABLE_OPTIMIZE_MACHINE)
    SET (OPT_FLAGS "-O5 -ansi -march=native -m64 -ffast-math -pipe")
ELSEIF (ENABLE_OPTIMIZE_X86_64)
    SET (OPT_FLAGS "-O5 -ansi -march=k8 -mfpmath=sse -m64 -ffast-math -pipe")
#    SET (OPT_FLAGS "-O5 -ansi -mfpmath=both -march=core-avx-i -m64 -ffast-math -pipe")
ELSEIF (ENABLE_OPTIMIZE_WINDOWS_VC)
	SET (OPT_FLAGS "/Ox /DHAVE_INLINE")
ELSEIF (ENABLE_OPTIMIZE_PHI)
	SET (OPT_FLAGS "-fast")
#	SET (OPT_FLAGS "-O3 -xHost")
ELSE()
    IF(WIN32)
       SET(OPT_FLAGS "")
    ELSE(WIN32)
       SET (OPT_FLAGS "-O3 -ansi -ffast-math")     # This is the default
    ENDIF(WIN32)
ENDIF()
#MESSAGE("OPT_FLAGS = " ${OPT_FLAGS})

IF(CMAKE_SYSTEM MATCHES "IRIX.*")
	SET(PLATFORMLIB "/usr/lib32" CACHE INTERNAL "IRIX lib32")
    #SET(EMAN_CXX_FLAGS "-n32 -LANG:std -w2 -fullwarn -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
    SET(EMAN_CXX_FLAGS "-n32 -woffall -LANG:std -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
ENDIF(CMAKE_SYSTEM MATCHES "IRIX.*")

IF(CMAKE_SYSTEM MATCHES "AIX*")
    SET(EMAN_CXX_FLAGS "" CACHE INTERNAL "C++ compiler flags")
ENDIF(CMAKE_SYSTEM MATCHES "AIX*")

IF(CMAKE_SYSTEM MATCHES "OSF1-V*")
    SET(EMAN_CXX_FLAGS "-w0" CACHE INTERNAL "C++ compiler flags")
ENDIF(CMAKE_SYSTEM MATCHES "OSF1-V*")

IF(CMAKE_COMPILER_IS_GNUCXX)
    #SET(EMAN_CXX_FLAGS "-Wall -Werror" CACHE INTERNAL "EMAN CXX FLAGS")
    IF (ENABLE_DEBUG)
    	SET(EMAN_CXX_FLAGS "-g -Wall -W -Wcast-qual -Wwrite-strings" CACHE INTERNAL "EMAN CXX FLAGS")
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g")
    ELSE (ENABLE_DEBUG)
    	SET(EMAN_CXX_FLAGS "-Wall" CACHE INTERNAL "EMAN CXX FLAGS")
    ENDIF (ENABLE_DEBUG)
    SET(PLATFORMLIB "/usr/lib64" CACHE INTERNAL "lib64")
	if(ENABLE_CONDA)
		SET(EMAN_CXX_FLAGS "${EMAN_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0" CACHE INTERNAL "EMAN CXX FLAGS")
	endif()
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

IF(ENABLE_IOCACHE)
	ADD_DEFINITIONS(-DIMAGEIO_CACHE)
ENDIF(ENABLE_IOCACHE)

# IF(ENABLE_HDFIOCACHE)
#     ADD_DEFINITIONS(-DHDFIO_CACHE)
# ENDIF(ENABLE_HDFIOCACHE)

SET(INCLUDE_DIRS
	${EMAN_SOURCE_DIR}/libEM
	$ENV{HOME}/EMAN2/anaconda2/include
	$ENV{HOME}/anaconda2/include
	$ENV{HOME}/include
	$ENV{FFTWDIR}/include
	$ENV{GSLDIR}/include
	/usr/local/include
	/usr/local/cuda/include
	/usr/include
	/sw/include
	/usr/include/freetype2
	)

SET(LIB_DIRS
	$ENV{EMANDIR}/lib
	$ENV{HOME}/EMAN2/anaconda2/lib
	$ENV{HOME}/anaconda2/lib
	$ENV{HOME}/lib
	${PLATFORMLIB}
	$ENV{FFTWDIR}/lib
	$ENV{GSLDIR}/lib
	/usr/local/lib
	/usr/local/cuda/lib
	/usr/lib64
	/usr/lib
	/sw/lib
	)


MACRO(CHECK_REQUIRED_LIB2 upper lower header lower2 header2)
	if(ENABLE_CONDA)
		message("\n### BEGIN ### CHECK_REQUIRED_LIB2 upper lower header lower2 header2")
		message_var(${upper}_INCLUDE_PATH)
		message_var(${upper}_LIBRARY)
	endif()
		
	IF(NOT ${upper}_INCLUDE_PATH)
		FIND_PATH(${upper}_INCLUDE_PATH 
				NAMES ${header} ${header2}
				PATHS $ENV{${upper}DIR}/include ${INCLUDE_DIRS}
				NO_DEFAULT_PATH
				)
    ENDIF()

    IF(${upper}_INCLUDE_PATH)
		FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} ${lower2} PATHS $ENV{${upper}DIR}/lib ${LIB_DIRS})
		message("FIND_LIBRARY: ${upper}_LIBRARY NAMES ${lower} ${lower2} PATHS $ENV{${upper}DIR}/lib ${LIB_DIRS}")
    ENDIF()
		
	IF(NOT ${upper}_INCLUDE_PATH OR NOT ${upper}_LIBRARY)
		MESSAGE(SEND_ERROR "ERROR: ${upper} not found. please install ${upper} first!")
    ENDIF()

	if(ENABLE_CONDA)
		message("### END ### CHECK_REQUIRED_LIB2 upper lower header lower2 header2")
		message_var(${upper}_INCLUDE_PATH)
		message_var(${upper}_LIBRARY)
	endif()
ENDMACRO(CHECK_REQUIRED_LIB2)

MACRO(CHECK_REQUIRED_LIB upper lower header)
    CHECK_REQUIRED_LIB2(${upper} ${lower} ${header} "" "")
ENDMACRO(CHECK_REQUIRED_LIB)

MACRO(CHECK_OPTIONAL_LIB upper lower header)
    FIND_PATH(${upper}_INCLUDE_PATH ${header} $ENV{${upper}DIR}/include ${INCLUDE_DIRS} NO_DEFAULT_PATH)

    IF(${upper}_INCLUDE_PATH)
        FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} PATHS $ENV{${upper}DIR}/lib ${LIB_DIRS})

        IF(${upper}_LIBRARY)
			ADD_DEFINITIONS(-DEM_${upper})
			INCLUDE_DIRECTORIES(${${upper}_INCLUDE_PATH})
        ELSE()
			#MESSAGE("WARNING: ${upper} not found, ENABLE_${upper} is set back to OFF")
			OPTION(ENABLE_${upper} "if enable ${upper} support" OFF)
        ENDIF()

    ELSE()
		#MESSAGE("WARNING: ${upper} not found, ENABLE_${upper} is set back to OFF")
        OPTION(ENABLE_${upper} "if enable ${upper} support" OFF)
    ENDIF()
ENDMACRO(CHECK_OPTIONAL_LIB)

MACRO(CHECK_LIB_ONLY upper lower)
    FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} PATHS $ENV{${upper}DIR}/lib ${LIB_DIRS} NO_DEFAULT_PATH)
	message(STATUS "CHECK_LIB_ONLY upper lower")
	message_var(${upper}_LIBRARY)
ENDMACRO(CHECK_LIB_ONLY)

IF(ENABLE_FFTW2)
    OPTION(ENABLE_NFFT "enable nfft support" OFF)
    OPTION(ENABLE_FFTW3 "enable fftw 3 support" OFF)
    OPTION(ENABLE_NATIVE_FFT "enable native fft support" OFF)
    OPTION(ENABLE_ACML_FFT "enable AMD Core Math Library fft support" OFF)

	IF (WIN32)
		CHECK_REQUIRED_LIB2(FFTW2 FFTW2dll sfftw.h fftw fftw.h)
		CHECK_REQUIRED_LIB2(RFFTW2 RFFTW2dll srfftw.h rfftw rfftw.hx)
	ELSE (WIN32)
		CHECK_REQUIRED_LIB2(FFTW2 sfftw sfftw.h fftw fftw.h)
		CHECK_REQUIRED_LIB2(RFFTW2 srfftw srfftw.h rfftw rfftw.hx)
	ENDIF (WIN32)

	ADD_DEFINITIONS(-DFFTW2)
    SET(FFTW_LIBRARIES ${RFFTW2_LIBRARY} ${FFTW2_LIBRARY})
    SET(FFTW_INCLUDE_PATH ${FFTW2_INCLUDE_PATH})

    IF(ENABLE_NFFT)
    	CHECK_REQUIRED_LIB(FFTW2D fftw fftw.h)
    	CHECK_REQUIRED_LIB(NFFT nfft nfft.h)
    	ADD_DEFINITIONS(-DNFFT)
		INCLUDE_DIRECTORIES(${NFFT_INCLUDE_PATH})
    ENDIF(ENABLE_NFFT)
ENDIF(ENABLE_FFTW2)

IF(ENABLE_FFTW3)
    OPTION(ENABLE_NFFT2 "enable nfft2 support" OFF)
    OPTION(ENABLE_FFTW2 "enable fftw 2 support" OFF)
    OPTION(ENABLE_NATIVE_FFT "enable native fft support" OFF)
    OPTION(ENABLE_ACML_FFT "enable AMD Core Math Library fft support" OFF)

    CHECK_REQUIRED_LIB2(FFTW3 fftw3f fftw3.h libfftw3f-3 "")
    CHECK_REQUIRED_LIB2(FFTW3d fftw3 fftw3.h libfftw3-3  "")
    ADD_DEFINITIONS(-DFFTW3)
    SET(FFTW_LIBRARIES ${FFTW3_LIBRARY} ${FFTW3d_LIBRARY})
    SET(FFTW_INCLUDE_PATH ${FFTW3_INCLUDE_PATH})

	if(ENABLE_CONDA)
		message("FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
		message("FFTW_INCLUDE_PATH: ${FFTW_INCLUDE_PATH}")
	endif()
		
    IF(ENABLE_NFFT2)
    	CHECK_REQUIRED_LIB(FFTW3D fftw3 fftw3.h)
    	CHECK_REQUIRED_LIB(NFFT2 nfft nfft.h)
    	ADD_DEFINITIONS(-DNFFT2)
		INCLUDE_DIRECTORIES(${NFFT2_INCLUDE_PATH})
    ENDIF(ENABLE_NFFT2)
ENDIF(ENABLE_FFTW3)

IF(ENABLE_FFTW_PLAN_CACHING)
	ADD_DEFINITIONS(-DFFTW_PLAN_CACHING)
ENDIF(ENABLE_FFTW_PLAN_CACHING)

IF(ENABLE_ACML_FFT)
	OPTION(ENABLE_FFTW2 "enable fftw 2 support" OFF)
	OPTION(ENABLE_FFTW3 "enable fftw 3 support" OFF)
	OPTION(ENABLE_NATIVE_FFT "enable native fft support" OFF)

	CHECK_REQUIRED_LIB(ACML acml acml.h)
	FIND_LIBRARY(G2C_LIBRARY NAMES g2c PATHS
		/usr/lib64
		/usr/lib
		/usr/local/lib
		$ENV{HOME}/lib )
	ADD_DEFINITIONS(-DACML)
	SET(ACML_LIBRARIES ${ACML_LIBRARY})
    SET(ACML_INCLUDE_PATH ${ACML_INCLUDE_PATH})
ENDIF(ENABLE_ACML_FFT)

IF (ENABLE_FFT_CACHING)
	ADD_DEFINITIONS(-DFFT_CACHING)
ENDIF(ENABLE_FFT_CACHING)

IF(ENABLE_NATIVE_FFT)
	ADD_DEFINITIONS(-DNATIVE_FFT)
	OPTION(ENABLE_FFTW2 "enable fftw 2 support" OFF)
	OPTION(ENABLE_FFTW3 "enable fftw 3 support" OFF)
	OPTION(ENABLE_ACML_FFT "enable AMD Core Math Library fft support" OFF)
ENDIF(ENABLE_NATIVE_FFT)

IF(ENABLE_DJBFFT)
	OPTION(ENABLE_NFFT2 "enable nfft2 support" OFF)
	CHECK_REQUIRED_LIB(DJBFFT fftc4.h fftr4.h)
	ADD_DEFINITIONS(-DDJBFFT)
	SET(DJBFFT_LIBRARIES ${DJBFFT_LIBRARY})
	SET(DJBFFT_INCLUDE_PATH ${DJBFFT_INCLUDE_PATH})
ENDIF(ENABLE_DJBFFT)

IF(ENABLE_OPTPP)
    FIND_PATH(OPTPP_INCLUDE_PATH Opt.h /usr/include /usr/include/opt++
              /usr/local/include /usr/local/include/opt++ $ENV{HOME}/include
             $ENV{HOME}/include/opt++ NO_DEFAULT_PATH)
    FIND_PATH(NEWMAT_INCLUDE_PATH newmat.h /usr/include /usr/include/opt++
              /usr/include/newmat /usr/local/include/newmat
              /usr/local/include /usr/local/include/opt++ $ENV{HOME}/include
             $ENV{HOME}/include/opt++ NO_DEFAULT_PATH)
    FIND_LIBRARY(OPTPP_LIBRARY NAMES opt-linux PATHS /usr/lib /usr/local/lib)
    FIND_LIBRARY(NEWMAT_LIBRARY NAMES newmat-linux PATHS /usr/lib /usr/local/lib)

    ADD_DEFINITIONS(-DOPTPP)
    INCLUDE_DIRECTORIES(${OPTPP_INCLUDE_PATH} ${NEWMAT_INCLUDE_PATH})
ENDIF(ENABLE_OPTPP)

CHECK_REQUIRED_LIB(GSL gsl gsl/gsl_linalg.h)
CHECK_REQUIRED_LIB(GSL_CBLAS gslcblas gsl/gsl_linalg.h)

IF (ENABLE_OPENGL)
	find_package(OpenGL REQUIRED)
	set(GL_LIBRARY       ${OPENGL_gl_LIBRARY} CACHE PATH "")
	set(GLU_LIBRARY      ${OPENGL_LIBRARIES}  CACHE PATH "")
	set(GL_INCLUDE_PATH  ${OPENGL_INCLUDE_DIR})
	set(GLU_INCLUDE_PATH ${OPENGL_INCLUDE_DIR})

	message("OPENGL_LIBRARY:     ${OPENGL_LIBRARY}")
	message("GL_LIBRARY:         ${GL_LIBRARY}")
	message("GLU_LIBRARY:        ${GLU_LIBRARY}")
	message("OPENGL_INCLUDE_DIR: ${OPENGL_INCLUDE_DIR}")
	message("GL_INCLUDE_PATH:    ${GL_INCLUDE_PATH}")
	message("GLU_INCLUDE_PATH:   ${GLU_INCLUDE_PATH}")

	ADD_DEFINITIONS(-DEMAN2_USING_OPENGL)
ENDIF (ENABLE_OPENGL)

IF (ENABLE_FTGL)
	CHECK_REQUIRED_LIB2(FTGL ftgl FTGL/FTGL.h ftgl FTGL/ftgl.h)
	ADD_DEFINITIONS(-DEMAN2_USING_FTGL)
	IF(EXISTS ${FTGL_INCLUDE_PATH}/FTGL/FTGL.h AND COMMAND IF)
		ADD_DEFINITIONS(-DOLD_FTGL)
	ENDIF(EXISTS ${FTGL_INCLUDE_PATH}/FTGL/FTGL.h AND COMMAND IF)
	ADD_SUBDIRECTORY(fonts)
	find_package(Freetype REQUIRED)

	message("FREETYPE_LIBRARIES:    ${FREETYPE_LIBRARIES}")
	message("FREETYPE_INCLUDE_DIRS: ${FREETYPE_INCLUDE_DIRS}")
		
	IF (ENABLE_STATIC_FTGL)
		IF(WIN32)
			CHECK_REQUIRED_LIB(FREETYPE freetype ft2build.h)
		ENDIF(WIN32)
	ENDIF (ENABLE_STATIC_FTGL)
ENDIF (ENABLE_FTGL)

IF(ENABLE_TIFF)
    CHECK_OPTIONAL_LIB(TIFF tiff tiffio.h)
    CHECK_LIB_ONLY(JPEG jpeg)
ENDIF(ENABLE_TIFF)

IF(ENABLE_JPEG)
    CHECK_OPTIONAL_LIB(JPEG jpeg jpeglib.h)
ENDIF(ENABLE_JPEG)

CHECK_OPTIONAL_LIB(HDF5 hdf5 hdf5.h)
CHECK_LIB_ONLY(ZLIB z)
IF (WIN32)
	CHECK_LIB_ONLY(SZLIB libszip)
ENDIF (WIN32)

IF(ENABLE_PNG)
	find_package(PNG REQUIRED)
	message("PNG_LIBRARIES: ${PNG_LIBRARIES}")
	ADD_DEFINITIONS(-DEM_PNG)
ENDIF(ENABLE_PNG)

execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
		OUTPUT_VARIABLE numpy_include_path
		OUTPUT_STRIP_TRAILING_WHITESPACE
		)

FIND_PATH(NUMPY_INCLUDE_PATH numpy/arrayobject.h
	${numpy_include_path}
	$ENV{HOME}/anaconda2/lib/python2.7/site-packages/numpy/core/include
	$ENV{HOME}/EMAN2/anaconda2/lib/python2.7/site-packages/numpy/core/include
	$ENV{PYTHON_ROOT}/lib/python2.7/site-packages/numpy/core/include
	/usr/local/lib/python2.7/site-packages/numpy/core/include
	/usr/lib64/python2.7/site-packages/numpy/core/include
	/usr/lib/python2.7/site-packages/numpy/core/include
    /usr/include 
    NO_DEFAULT_PATH)

message("NUMPY_INCLUDE_PATH: ${NUMPY_INCLUDE_PATH}")

INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_PATH})


# Find Boost
find_package(Boost COMPONENTS python REQUIRED)

message("Boost_LIBRARIES:   ${Boost_LIBRARIES}")
message("Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")

# IF(ENABLE_HDFIOCACHE)
#    IF(APPLE)
#        CHECK_OPTIONAL_LIB(BOOST_THREAD boost_thread boost/thread.hpp)
#        CHECK_OPTIONAL_LIB(BOOST_FILESYSTEM boost_filesystem boost/filesystem.hpp)
#        CHECK_OPTIONAL_LIB(BOOST_SYSTEM boost_system boost/filesystem.hpp)
#    ELSE(APPLE)
#        CHECK_OPTIONAL_LIB(BOOST_THREAD boost_thread-mt boost/thread.hpp)
#        CHECK_OPTIONAL_LIB(BOOST_FILESYSTEM boost_filesystem-mt boost/filesystem.hpp)
#        CHECK_OPTIONAL_LIB(BOOST_SYSTEM boost_system-mt boost/filesystem.hpp)
#    ENDIF(APPLE)
# ENDIF(ENABLE_HDFIOCACHE)

INCLUDE_DIRECTORIES(.
	    ${EMAN_SOURCE_DIR}
	    ${EMAN_SOURCE_DIR}/libEM
	    ${GL_INCLUDE_PATH}
	    ${GSL_INCLUDE_PATH}
	    ${GSL_INCLUDE_PATH}/gsl
	    ${FFTW_INCLUDE_PATH}
	    ${DJBFFT_INCLUDE_PATH}
	    ${ACML_INCLUDE_PATH}
	    ${Boost_INCLUDE_DIR}
	    ${PYTHON_INCLUDE_PATH}
	    ${FREETYPE_INCLUDE_DIRS}
            )

#---------------------------------------------------------------------------------------
if(ENABLE_CONDA)
	message("\nINCLUDE_DIRECTORIES---START")
	get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
	
	foreach(dir ${dirs})
		message(STATUS "INCLUDE_DIR='${dir}'")
	endforeach()
	message("\nINCLUDE_DIRECTORIES---END")
endif()
#---------------------------------------------------------------------------------------

IF(ENABLE_RT)
    ADD_SUBDIRECTORY(rt)
ENDIF(ENABLE_RT)

IF (ENABLE_SPARX_CUDA)
	INCLUDE(${CMAKE_SOURCE_DIR}/libEM/cuda/FindCUDA.cmake)
	CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/libEM)
	CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/libEM/sparx/cuda)
	ADD_DEFINITIONS(-DSPARX_USING_CUDA)
	INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})
	ADD_SUBDIRECTORY(libEM/sparx/cuda)
ENDIF (ENABLE_SPARX_CUDA)

IF (ENABLE_EMAN_CUDA)
	INCLUDE(${CMAKE_SOURCE_DIR}/libEM/cuda/FindCUDA.cmake)
	CUDA_INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/libEM)
	ADD_DEFINITIONS(-DEMAN2_USING_CUDA)
	INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})
	ADD_SUBDIRECTORY(libEM/cuda)
ENDIF (ENABLE_EMAN_CUDA)

IF(ENABLE_TAO)
    ADD_SUBDIRECTORY(tao)
ENDIF(ENABLE_TAO)

ADD_SUBDIRECTORY(libEM)
ADD_SUBDIRECTORY(libpyEM)
ADD_SUBDIRECTORY(pyemtbx)
ADD_SUBDIRECTORY(sparx)
ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(examples)
ADD_SUBDIRECTORY(programs)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(images)
ADD_SUBDIRECTORY(recipes)

if(NOT WIN32)
	if(ENABLE_CONDA)
		install(CODE "execute_process(
						COMMAND ${CMAKE_COMMAND} -E create_symlink ${EMAN_PREFIX}/bin/e2version.py ${SP_DIR}/e2version.py
						COMMAND ${CMAKE_COMMAND} -E create_symlink ${EMAN_PREFIX}/bin/sxgui.py     ${EMAN_PREFIX}/bin/sphire
						COMMAND ${CMAKE_COMMAND} -E create_symlink ${EMAN_PREFIX}/bin/sx.py        ${EMAN_PREFIX}/bin/sparx
					 )"
				)
	endif()
endif()

enable_testing()

add_custom_target(test-verbose
		COMMAND ${CMAKE_COMMAND} -P cmake_install.cmake
		COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release
		)

#add_custom_target(test-verbose
#		COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release
#		DEPENDS Tests
#		)

function(add_test_wrapper name)
	add_test(NAME ${name}
			COMMAND nosetests ${EMAN_INSTALL_PREFIX}/test/rt/pyem/${name}
			)
endfunction()

set(test_methods
		test_imageio.py:TestHdfIO.test_read_image
		test_imageio.py:TestHdfIO.test_write_image
		test_imageio.py:TestHdfIO.test_read_write_hdf
		test_imageio.py:TestPNGIO.test_write_png
		)

foreach(t ${test_methods})
	add_test_wrapper(${t})
endforeach()
# --------------------------------------------------------------------
#set the list of unit test to run
#
SET (TESTS
	test_vec
	test_utils
	test_typeconverter
	test_transform
	test_misc
	test_imageio
	test_emdata
	test_processor
	test_cmp
	test_aligner
	test_averager
	)

#add the custom commands for each unit test program
#
FOREACH (TEST ${TESTS})

	ADD_CUSTOM_COMMAND(
        OUTPUT ${TEST}.PRINT_SPACE
        COMMAND ${CMAKE_COMMAND} ARGS -E echo ""
	)

	#this is the trick to make the target run every time, DoUnitTest suppose is a file never been built
	#
	ADD_CUSTOM_COMMAND(
        OUTPUT ${TEST}.DoUnitTest
		COMMAND ${CMAKE_COMMAND} ARGS -E time ${EMAN_INSTALL_PREFIX}/test/rt/${TEST}.py
		DEPENDS	${TEST}.PRINT_SPACE
	)

	#build a list of all the results
	SET( TEST_RESULTS ${TEST_RESULTS}
		${TEST}.DoUnitTest
	)

ENDFOREACH (TEST)
