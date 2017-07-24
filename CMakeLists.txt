CMAKE_MINIMUM_REQUIRED(VERSION 3.7)
PROJECT(EMAN)

IF(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
       cmake_policy(SET CMP0002 OLD)
ENDIF(COMMAND cmake_policy)

include(cmake/functions.cmake)

OPTION(DEBUG_CMAKE_TARGETS "enable debug output for cmake target properties" OFF)
if(DEBUG_CMAKE_TARGETS)
	set(CMAKE_DEBUG_TARGET_PROPERTIES
		INCLUDE_DIRECTORIES
		COMPILE_DEFINITIONS
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

SET(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

MARK_AS_ADVANCED(CLEAR CMAKE_VERBOSE_MAKEFILE)
OPTION(CMAKE_VERBOSE_MAKEFILE "if all commands will be echoed to the console during the make" ON)

OPTION(ENABLE_CONDA "enable conda-supported build" ON)

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
		find_package(Conda REQUIRED)
			
		if(CONDA_ROOT)
			if(WIN32)
				set(EMAN_PREFIX ${CONDA_ROOT}/Library CACHE PATH "")
			else()
				set(EMAN_PREFIX ${CONDA_ROOT} CACHE PATH "")
			endif()
		else()
			set(CONDA_ROOT CONDA_ROOT-NOTFOUND   CACHE PATH "")
		endif()
	endif()

	set(CMAKE_PREFIX_PATH ${EMAN_PREFIX} ${EMAN_PREFIX}/..)
	message_var(CMAKE_PREFIX_PATH)
endif()

find_package(Python REQUIRED)
find_package(NumPy  REQUIRED)

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

if(ENABLE_CONDA)
	SET(EMAN_INSTALL_PREFIX    ${EMAN_PREFIX} CACHE PATH "installation prefix")
	set(EMAN_PY_INSTALL_PREFIX ${SP_DIR})
	
	set(EMAN_PREFIX_INC ${EMAN_PREFIX}/include)
	set(EMAN_PREFIX_LIB ${EMAN_PREFIX}/lib)
	
	if(NOT WIN32)
		set(eman2dir_relative_path "../../../")
	else()
		set(eman2dir_relative_path "../../Library/")
	endif()
else()
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

OPTION(ENABLE_FFTW3 "enable fftw 3 support (USE THIS)" ON)
OPTION(ENABLE_FFTW_PLAN_CACHING "enable fftw caching" ON)
OPTION(ENABLE_NATIVE_FFT "enable native fft support (for non-GPL use)" OFF)
OPTION(ENABLE_ACML_FFT "enable AMD Core Math Library fft support" OFF)
OPTION(ENABLE_FFT_CACHING "enable FFT Caching" OFF)
MARK_AS_ADVANCED(ENABLE_ACML_FFT)
OPTION(ENABLE_DJBFFT "enable djbfft support" OFF)
MARK_AS_ADVANCED(ENABLE_DJBFFT)

OPTION(ENABLE_OPTPP "enable Opt++ support" OFF)

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

IF(ENABLE_MEMDEBUG)
	ADD_DEFINITIONS(-DMEMDEBUG -W -Wall)
ENDIF()

IF(ENABLE_MEMDEBUG2)
	ADD_DEFINITIONS(-DMEMDEBUG2 -W -Wall)
ENDIF()

IF(ENABLE_V4L2)
	ADD_DEFINITIONS(-DBOUNDS_CHECKING -DENABLE_V4L2)
ENDIF()

IF(ENABLE_BOUNDS_CHECKING)
	ADD_DEFINITIONS(-DBOUNDS_CHECKING)
ENDIF()

IF(ENABLE_DEBUG)
    SET (OPT_FLAGS "-DDEBUG -O2 -g -ansi -Wall -W -fno-inline")
ELSEIF(ENABLE_PROFILE)
    SET (OPT_FLAGS "-O2 -ansi -march=native -m64 -ffast-math -pipe")
    OPTION(ENABLE_DEBUG "enable debug support" OFF)
ELSEIF(ENABLE_OPTIMIZE_OSX)
#    SET (OPT_FLAGS "-O3 -ffast-math")
    SET (OPT_FLAGS "-O3")
#    SET (OPT_FLAGS "-Os -ffast-math")
ELSEIF(ENABLE_OPTIMIZE_MACHINE)
    SET (OPT_FLAGS "-O5 -ansi -march=native -m64 -ffast-math -pipe")
ELSEIF(ENABLE_OPTIMIZE_X86_64)
    SET (OPT_FLAGS "-O5 -ansi -march=k8 -mfpmath=sse -m64 -ffast-math -pipe")
#    SET (OPT_FLAGS "-O5 -ansi -mfpmath=both -march=core-avx-i -m64 -ffast-math -pipe")
ELSEIF(ENABLE_OPTIMIZE_WINDOWS_VC)
	SET (OPT_FLAGS "/Ox /DHAVE_INLINE")
ELSEIF(ENABLE_OPTIMIZE_PHI)
	SET (OPT_FLAGS "-fast")
#	SET (OPT_FLAGS "-O3 -xHost")
ELSE()
    IF(WIN32)
       SET(OPT_FLAGS "")
    ELSE()
       SET (OPT_FLAGS "-O3 -ansi -ffast-math")     # This is the default
    ENDIF()
ENDIF()
#MESSAGE("OPT_FLAGS = " ${OPT_FLAGS})

IF(CMAKE_SYSTEM MATCHES "IRIX.*")
	SET(PLATFORMLIB "/usr/lib32" CACHE INTERNAL "IRIX lib32")
    #SET(EMAN_CXX_FLAGS "-n32 -LANG:std -w2 -fullwarn -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
    SET(EMAN_CXX_FLAGS "-n32 -woffall -LANG:std -LANG:ansi-for-init-scope=ON" CACHE INTERNAL "C++ compiler flags")
ENDIF()

IF(CMAKE_SYSTEM MATCHES "AIX*")
    SET(EMAN_CXX_FLAGS "" CACHE INTERNAL "C++ compiler flags")
ENDIF()

IF(CMAKE_SYSTEM MATCHES "OSF1-V*")
    SET(EMAN_CXX_FLAGS "-w0" CACHE INTERNAL "C++ compiler flags")
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
    #SET(EMAN_CXX_FLAGS "-Wall -Werror" CACHE INTERNAL "EMAN CXX FLAGS")
    IF(ENABLE_DEBUG)
    	SET(EMAN_CXX_FLAGS "-g -Wall -W -Wcast-qual -Wwrite-strings" CACHE INTERNAL "EMAN CXX FLAGS")
	SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g")
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g")
    ELSE()
    	SET(EMAN_CXX_FLAGS "-Wall" CACHE INTERNAL "EMAN CXX FLAGS")
    ENDIF()
    SET(PLATFORMLIB "/usr/lib64" CACHE INTERNAL "lib64")
	if(ENABLE_CONDA)
		SET(EMAN_CXX_FLAGS "${EMAN_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0" CACHE INTERNAL "EMAN CXX FLAGS")
	endif()
ENDIF()

SET(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER})
SET(CMAKE_ORIG_FLAGS ${CMAKE_CXX_FLAGS})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EMAN_CXX_FLAGS} ${OPT_FLAGS}")

MARK_AS_ADVANCED(CMAKE_C_FLAGS)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS)

IF(ENABLE_IOCACHE)
	ADD_DEFINITIONS(-DIMAGEIO_CACHE)
ENDIF()


IF(ENABLE_FFTW3)
    OPTION(ENABLE_NFFT2 "enable nfft2 support" OFF)
    OPTION(ENABLE_NATIVE_FFT "enable native fft support" OFF)
    OPTION(ENABLE_ACML_FFT "enable AMD Core Math Library fft support" OFF)

	find_package(FFTW3 REQUIRED)
		
    IF(ENABLE_NFFT2)
		find_package(NFFT2 REQUIRED)
		target_link_libraries(NFFT2 FFTW3D::FFTW3D)
	ENDIF()
ENDIF()

IF(ENABLE_FFTW_PLAN_CACHING)
	ADD_DEFINITIONS(-DFFTW_PLAN_CACHING)
ENDIF()

IF(ENABLE_ACML_FFT)
	OPTION(ENABLE_FFTW3 "enable fftw 3 support" OFF)
	OPTION(ENABLE_NATIVE_FFT "enable native fft support" OFF)

	find_package(ACML REQUIRED)
ENDIF()

IF(ENABLE_FFT_CACHING)
	ADD_DEFINITIONS(-DFFT_CACHING)
ENDIF()

IF(ENABLE_NATIVE_FFT)
	ADD_DEFINITIONS(-DNATIVE_FFT)
	OPTION(ENABLE_FFTW3 "enable fftw 3 support" OFF)
	OPTION(ENABLE_ACML_FFT "enable AMD Core Math Library fft support" OFF)
ENDIF()

IF(ENABLE_DJBFFT)
	OPTION(ENABLE_NFFT2 "enable nfft2 support" OFF)
		
	find_package(DJBFFT REQUIRED)
ENDIF()

IF(ENABLE_OPTPP)
	find_package(OPTPP REQUIRED)
ENDIF()

IF(ENABLE_OPENGL)
	include(cmake/OpenGL.cmake)
ENDIF()

IF(ENABLE_FTGL)
	find_package(FTGL REQUIRED)

	ADD_SUBDIRECTORY(fonts)
	include(cmake/Freetype.cmake)
ENDIF()

IF(ENABLE_TIFF)
	include(cmake/TIFF.cmake)
ENDIF()

IF(ENABLE_JPEG)
    include(cmake/JPEG.cmake)
ENDIF()

include(cmake/HDF5.cmake)
IF(WIN32)
	find_package(SZLIB REQUIRED)
ENDIF()

IF(ENABLE_PNG)
	include(cmake/PNG.cmake)
ENDIF()

IF(ENABLE_RT)
    ADD_SUBDIRECTORY(rt)
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

set(test_methods
		test_imageio.py:TestHdfIO.test_read_image
		test_imageio.py:TestHdfIO.test_write_image
		test_imageio.py:TestHdfIO.test_read_write_hdf
		test_imageio.py:TestPNGIO.test_write_png
		)

foreach(t ${test_methods})
	add_test(NAME ${t}
			COMMAND nosetests ${EMAN_INSTALL_PREFIX}/test/rt/pyem/${t}
			)
endforeach()

file(READ ${CMAKE_SOURCE_DIR}/tests/test_imports.py py_imports_list)
STRING(REGEX REPLACE "\n" ";" py_imports_list "${py_imports_list}")
foreach(i ${py_imports_list})
	set(import_statement ${i})
	separate_arguments(import_statement)
	list(GET import_statement 1 item1)
	list(GET import_statement 3 item2)
	if(NOT WIN32 OR NOT "${item1}" STREQUAL "libpyTomoSeg2")
		add_test(NAME import-${item1}.${item2}
				 COMMAND ${PYTHON_EXECUTABLE} -c "${i}"
				 )
	endif()
endforeach()

add_test(NAME import-inners
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_imports_inner.py
		 )

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
FOREACH(TEST ${TESTS})

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

ENDFOREACH()
