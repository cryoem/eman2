CMAKE_MINIMUM_REQUIRED(VERSION 3.7)
PROJECT(EMAN
		VERSION 2.2
		LANGUAGES CXX
		)

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
message_var(CMAKE_PREFIX_PATH)

find_package(Python REQUIRED)
find_package(NumPy  REQUIRED)

# Set SP_DIR
if("$ENV{CONDA_BUILD_STATE}" STREQUAL "BUILD" )
	set(SP_DIR $ENV{SP_DIR})
else()
	if(NOT WIN32)
		set(py_sp_dir_command "import site; print(site.getsitepackages()[0])")
	else()
		set(py_sp_dir_command "import site; print(site.getsitepackages()[1])")
	endif()
	execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "${py_sp_dir_command}"
			OUTPUT_VARIABLE SP_DIR
			OUTPUT_STRIP_TRAILING_WHITESPACE
			)
	message("Python site-packages: ${SP_DIR}")
endif()

# Write relative path that will be used to compute EMAN2DIR in libpyEM/EMAN2.py
if(NOT WIN32)
	set(eman2dir_relative_path "../../../")
else()
	set(eman2dir_relative_path "../../Library/")
endif()

if(SP_DIR)
	file(WRITE ${SP_DIR}/eman2dir_relative_path_to_sp_dir ${eman2dir_relative_path})
endif()

set(CMAKE_INSTALL_RPATH ${SP_DIR})

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
	SET(EMAN_CXX_FLAGS "${EMAN_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0" CACHE INTERNAL "EMAN CXX FLAGS")
ENDIF()

SET(CMAKE_C_COMPILER ${CMAKE_CXX_COMPILER})
SET(CMAKE_ORIG_FLAGS ${CMAKE_CXX_FLAGS})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EMAN_CXX_FLAGS} ${OPT_FLAGS}")

if(WIN32)
	message_var(MSVC_VERSION)
endif()

if(MSVC_VERSION GREATER 1500)
	add_definitions(-DNOMINMAX)
endif()

MARK_AS_ADVANCED(CMAKE_C_FLAGS)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS)

IF(ENABLE_IOCACHE)
	ADD_DEFINITIONS(-DIMAGEIO_CACHE)
ENDIF()


IF(ENABLE_FFTW_PLAN_CACHING)
	ADD_DEFINITIONS(-DFFTW_PLAN_CACHING)
ENDIF()

IF(ENABLE_FFT_CACHING)
	ADD_DEFINITIONS(-DFFT_CACHING)
ENDIF()

if(NOT(    (    ENABLE_FFTW3 AND NOT ENABLE_NATIVE_FFT AND NOT ENABLE_ACML_FFT)
   		OR (NOT ENABLE_FFTW3 AND     ENABLE_NATIVE_FFT AND NOT ENABLE_ACML_FFT)
   		OR (NOT ENABLE_FFTW3 AND NOT ENABLE_NATIVE_FFT AND     ENABLE_ACML_FFT)
   		)
   )
	message(FATAL_ERROR "Exactly one of the below has to be ON:\nENABLE_FFTW3= ${ENABLE_FFTW3}\nENABLE_NATIVE_FFT= ${ENABLE_NATIVE_FFT}\nENABLE_ACML_FFT= ${ENABLE_ACML_FFT}\n")
endif()

IF(ENABLE_NATIVE_FFT)
	ADD_DEFINITIONS(-DNATIVE_FFT)
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

add_custom_target(PythonFiles
				  COMMAND ${CMAKE_COMMAND} -DCOMPONENT=PythonFiles -P cmake_install.cmake
				  )

if(NOT WIN32)
	install(CODE "execute_process(
					COMMAND ${CMAKE_COMMAND} -E create_symlink ${EMAN_PREFIX}/bin/sxgui.py ${EMAN_PREFIX}/bin/sphire
					COMMAND ${CMAKE_COMMAND} -E create_symlink ${EMAN_PREFIX}/bin/sx.py    ${EMAN_PREFIX}/bin/sparx
				 )"
			)
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

add_test(NAME imports
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_imports.py
		 )

add_test(NAME test-EMAN2DIR
		 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_EMAN2DIR.py
		 )

find_package(Nosetests)

if(NOT WIN32)
	add_test(NAME nose-tests
			COMMAND ${NOSETESTS_EXECUTABLE} -vv --exe -m "^test_*" -e "^test_image_" -e "test_main" -e "test_result" -e "test_boxing" -a \!broken
					${CMAKE_SOURCE_DIR}/rt/pyem/ --with-coverage --cover-html --cover-html-dir=covnose --cover-erase
					--cover-package=EMAN2
					--cover-package=EMAN2_cppwrap
					--cover-package=EMAN2_meta
					--cover-package=EMAN2db
					--cover-package=EMAN2jsondb
					--cover-package=pyemtbx
			)
else()
	set(test_methods-win
			test_imageio.py:TestHdfIO.test_read_image
			test_imageio.py:TestHdfIO.test_write_image
			test_imageio.py:TestHdfIO.test_read_write_hdf
			test_imageio.py:TestPNGIO.test_write_png
		)

	foreach(t ${test_methods-win})
		add_test(NAME ${t}
				COMMAND ${NOSETESTS_EXECUTABLE} -v -m "^test_*" -a \!broken
						${CMAKE_SOURCE_DIR}/rt/pyem/${t}
				)
	endforeach()
endif()

add_custom_target(test-rt
		COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release -R nose-tests
		DEPENDS PythonFiles
		)

add_custom_target(test-verbose-broken
		COMMAND ${NOSETESTS_EXECUTABLE} -v -m "^test_*" -a broken ${CMAKE_SOURCE_DIR}/rt/pyem/*
		)

if(NOT WIN32)
	add_custom_target(test-progs
					  COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release -R progs
					  DEPENDS PythonFiles
					  )
	
	add_test(NAME progs
			 COMMAND bash ${CMAKE_SOURCE_DIR}/tests/run_prog_tests.sh
			 )
endif()

add_test(NAME py-compile
		COMMAND ${PYTHON_EXECUTABLE} -m compileall -q ${CMAKE_SOURCE_DIR}
		)

add_custom_target(test-py-compile
		COMMAND ${CMAKE_CTEST_COMMAND} -V -C Release -R py-compile
		DEPENDS PythonFiles
		)
