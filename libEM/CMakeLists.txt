ADD_SUBDIRECTORY(pmconfig)

SET(LIB_EM_SOURCE byteorder.cpp
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
    sparx/emdata_sparx.cpp
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
    sparx/fundamentals.cpp
    sparx/fourierfilter.cpp
    sparx/rsconvolution.cpp
    polardata.cpp
    sparx/native_fft.cpp
    sparx/util_sparx.cpp
    sparx/lapackblas.cpp
    sparx/pca.cpp
    sparx/varimax.cpp
    sparx/lbfgsb.cpp
    sparx/steepest.cpp
    gorgon/grid_queue.cpp
    gorgon/grid_queue2.cpp
    gorgon/volume_data.cpp
    gorgon/volume.cpp
    gorgon/skeletonizer.cpp
    tomoseg.cpp
    )

SET(LIB_GLEM_SOURCE 
		glutil.cpp
		marchingcubes.cpp
		emftgl.cpp)

OPTION(ENABLE_FFTW3 "enable fftw 3 support (USE THIS)" ON)
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

IF(WIN32)
    ADD_LIBRARY(EM2 STATIC ${LIB_EM_SOURCE})
    ADD_LIBRARY(GLEM2 STATIC ${LIB_GLEM_SOURCE})
    TARGET_LINK_LIBRARIES(EM2 FFTW3::FFTW3)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
ELSE()
    ADD_LIBRARY(EM2 SHARED ${LIB_EM_SOURCE})
    ADD_LIBRARY(GLEM2 SHARED ${LIB_GLEM_SOURCE})
    IF(ENABLE_FFTW3)
    	TARGET_LINK_LIBRARIES(EM2 FFTW3::FFTW3)
    ENDIF()
	TARGET_LINK_LIBRARIES(EM2 m)
ENDIF()

find_package(GSL REQUIRED)
if(WIN32)
	set_target_properties(GSL::gsl
						  PROPERTIES
						  INTERFACE_COMPILE_DEFINITIONS GSL_DLL
						  )
endif()

TARGET_LINK_LIBRARIES(EM2 GSL::gsl GSL::gslcblas Freetype::Freetype)
TARGET_LINK_LIBRARIES(GLEM2 EM2)

OPTION(ENABLE_OPENGL "enable EMAN2 OpenGL usage" ON)
IF(ENABLE_OPENGL)
	include(${CMAKE_SOURCE_DIR}/cmake/OpenGL.cmake)
	target_link_libraries(EM2 OpenGL::OpenGL)
ENDIF()

IF(ENABLE_FTGL)
	TARGET_LINK_LIBRARIES(GLEM2 FTGL::FTGL)
ENDIF()

OPTION(ENABLE_TIFF "enable TIFF support" ON)
IF(ENABLE_TIFF)
	include(${CMAKE_SOURCE_DIR}/cmake/TIFF.cmake)
    TARGET_LINK_LIBRARIES(EM2 TIFF::TIFF JPEG::JPEG)
ENDIF()

OPTION(ENABLE_JPEG "enable JPEG support" ON)
IF(ENABLE_JPEG)
	include(${CMAKE_SOURCE_DIR}/cmake/JPEG.cmake)
    TARGET_LINK_LIBRARIES(EM2 JPEG::JPEG)
ENDIF()

include(${CMAKE_SOURCE_DIR}/cmake/HDF5.cmake)
TARGET_LINK_LIBRARIES(EM2 HDF5::HDF5)
IF(WIN32)
	find_package(SZLIB REQUIRED)
	TARGET_LINK_LIBRARIES(EM2 SZLIB::SZLIB)
ENDIF()

OPTION(ENABLE_PNG "enable PNG support" ON)
IF(ENABLE_PNG)
	include(${CMAKE_SOURCE_DIR}/cmake/PNG.cmake)
    TARGET_LINK_LIBRARIES(EM2 PNG::PNG)
ENDIF()

IF(ENABLE_SPARX_CUDA)
    TARGET_LINK_LIBRARIES(EM2 EM2SparxCuda FFTW3::FFTW3)
ENDIF()

IF(ENABLE_EMAN_CUDA)
    TARGET_LINK_LIBRARIES(EM2 EM2Cuda FFTW3::FFTW3)
ENDIF()

IF(ENABLE_NFFT2)
	SET(LIB_EM_SOURCE ${LIB_EM_SOURCE} pointarray.cpp)
	TARGET_LINK_LIBRARIES(EM2 NFFT2::NFFT2 FFTW3D::FFTW3D)
ENDIF()

OPTION(ENABLE_DJBFFT "enable djbfft support" OFF)
MARK_AS_ADVANCED(ENABLE_DJBFFT)
IF(ENABLE_DJBFFT)
	OPTION(ENABLE_NFFT2 "enable nfft2 support" OFF)
	
	find_package(DJBFFT REQUIRED)
	TARGET_LINK_LIBRARIES(EM2 DJBFFT::DJBFFT)
ENDIF()

OPTION(ENABLE_OPTPP "enable Opt++ support" OFF)
IF(ENABLE_OPTPP)
	find_package(OPTPP REQUIRED)
	TARGET_LINK_LIBRARIES(EM2 OPTPP::OPTPP)
ENDIF()

OPTION(ENABLE_ACML_FFT "enable AMD Core Math Library fft support" OFF)
MARK_AS_ADVANCED(ENABLE_ACML_FFT)
IF(ENABLE_ACML_FFT)
	OPTION(ENABLE_FFTW3 "enable fftw 3 support" OFF)
	OPTION(ENABLE_NATIVE_FFT "enable native fft support" OFF)
	
	find_package(ACML REQUIRED)
	TARGET_LINK_LIBRARIES(EM2 ACML::ACML)
ENDIF()

IF(ENABLE_OPENGL)
    INSTALL(TARGETS GLEM2
        DESTINATION ${EMAN_PY_INSTALL_PREFIX})
ENDIF()

INSTALL(TARGETS EM2
		DESTINATION ${EMAN_PY_INSTALL_PREFIX})

