add_subdirectory(pmconfig)

set(LIB_EM_SOURCE byteorder.cpp
    boxingtools.cpp
    emobject.cpp
    emfft.cpp
    log.cpp
    imageio.cpp
    util.cpp
    emutil.cpp
    interp.cpp
    quaternion.cpp
    transform.cpp
    symmetry.cpp
    geometry.cpp
    pointarray.cpp
    pdbreader.cpp
    mrcio.cpp
    dm3io.cpp
    dm4io.cpp
    tifio.cpp
    hdfio.cpp
    hdfio2.cpp
    jpegio.cpp
    emdata.cpp
    emdata_io.cpp
    emdata_core.cpp
    emdata_cuda.cpp
    emdata_modular.cpp
    emdata_metadata.cpp
    emdata_transform.cpp
    pifio.cpp
    v4l2io.cpp
    vtkio.cpp
    spiderio.cpp
    sspiderio.cpp
    pgmio.cpp
    fitsio.cpp
    icosio.cpp
    lstio.cpp
    lstfastio.cpp
    pngio.cpp
    salio.cpp
    amiraio.cpp
    xplorio.cpp
    gatan2io.cpp
    emio.cpp
    imagicio.cpp
    imagicio2.cpp
    df3io.cpp
    omapio.cpp
    situsio.cpp
    serio.cpp
    emcache.cpp
    ctf.cpp
    xydata.cpp
    processor.cpp
    aligner.cpp
    projector.cpp
    cmp.cpp
    averager.cpp
    reconstructor.cpp
    reconstructor_tools.cpp
    exception.cpp
    testutil.cpp
    analyzer.cpp
    randnum.cpp
    hdf_filecache.cpp
    polardata.cpp
    tomoseg.cpp
    )

set(LIB_GLEM_SOURCE 
		glutil.cpp
		marchingcubes.cpp
		emftgl.cpp)

if(ENABLE_FFTW3)
	find_package(FFTW3 REQUIRED)
	
	if(ENABLE_NFFT2)
		find_package(NFFT2 REQUIRED)
		target_link_libraries(NFFT2 FFTW3D::FFTW3D)
	endif()
endif()

if(WIN32)
    add_library(EM2   STATIC "")
    add_library(GLEM2 STATIC "")
    target_link_libraries(EM2 FFTW3::FFTW3)
	target_compile_definitions(EM2 PUBLIC _CRT_SECURE_NO_WARNINGS _SCL_SECURE_NO_WARNINGS)
else()
    add_library(EM2   SHARED "")
    add_library(GLEM2 SHARED "")
    if(ENABLE_FFTW3)
    	target_link_libraries(EM2 FFTW3::FFTW3)
    endif()
	target_link_libraries(EM2 m)
endif()

target_sources(EM2   PRIVATE ${LIB_EM_SOURCE})
target_sources(GLEM2 PRIVATE ${LIB_GLEM_SOURCE})

add_subdirectory(gorgon)
add_subdirectory(sparx)

find_package(GSL REQUIRED)
if(WIN32)
	set_target_properties(GSL::gsl
						  PROPERTIES
						  INTERFACE_COMPILE_DEFINITIONS GSL_DLL
						  )
endif()

target_link_libraries(EM2 GSL::gsl GSL::gslcblas Freetype::Freetype)
target_link_libraries(GLEM2 EM2)

option(ENABLE_OPENGL "enable EMAN2 OpenGL usage" ON)
if(ENABLE_OPENGL)
	include(${CMAKE_SOURCE_DIR}/cmake/OpenGL.cmake)
	target_link_libraries(EM2 OpenGL::OpenGL)
endif()

IF(ENABLE_FTGL)
	TARGET_LINK_LIBRARIES(GLEM2 FTGL::FTGL)
ENDIF()

option(ENABLE_TIFF "enable TIFF support" ON)
if(ENABLE_TIFF)
	include(${CMAKE_SOURCE_DIR}/cmake/TIFF.cmake)
    target_link_libraries(EM2 TIFF::TIFF JPEG::JPEG)
endif()

option(ENABLE_JPEG "enable JPEG support" ON)
if(ENABLE_JPEG)
	include(${CMAKE_SOURCE_DIR}/cmake/JPEG.cmake)
    target_link_libraries(EM2 JPEG::JPEG)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/HDF5.cmake)
target_link_libraries(EM2 HDF5::HDF5)
if(WIN32)
	find_package(SZLIB REQUIRED)
	target_link_libraries(EM2 SZLIB::SZLIB)
endif()

option(ENABLE_PNG "enable PNG support" ON)
if(ENABLE_PNG)
	include(${CMAKE_SOURCE_DIR}/cmake/PNG.cmake)
    target_link_libraries(EM2 PNG::PNG)
endif()

if(ENABLE_SPARX_CUDA)
    target_link_libraries(EM2 EM2SparxCuda FFTW3::FFTW3)
endif()

if(ENABLE_EMAN_CUDA)
    target_link_libraries(EM2 EM2Cuda FFTW3::FFTW3)
endif()

if(ENABLE_NFFT2)
	set(LIB_EM_SOURCE ${LIB_EM_SOURCE} pointarray.cpp)
	target_link_libraries(EM2 NFFT2::NFFT2 FFTW3D::FFTW3D)
endif()

option(ENABLE_DJBFFT "enable djbfft support" OFF)
mark_as_advanced(ENABLE_DJBFFT)
if(ENABLE_DJBFFT)
	option(ENABLE_NFFT2 "enable nfft2 support" OFF)
	
	find_package(DJBFFT REQUIRED)
	target_link_libraries(EM2 DJBFFT::DJBFFT)
endif()

option(ENABLE_OPTPP "enable Opt++ support" OFF)
if(ENABLE_OPTPP)
	find_package(OPTPP REQUIRED)
	target_link_libraries(EM2 OPTPP::OPTPP)
endif()

if(ENABLE_ACML_FFT)
	find_package(ACML REQUIRED)
	target_link_libraries(EM2 ACML::ACML)
endif()

if(ENABLE_OPENGL)
    install(TARGETS GLEM2
        DESTINATION ${EMAN_PY_INSTALL_PREFIX})
endif()

install(TARGETS EM2
		DESTINATION ${EMAN_PY_INSTALL_PREFIX})

