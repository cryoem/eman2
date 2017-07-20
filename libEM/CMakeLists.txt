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

IF(WIN32)
    ADD_LIBRARY(EM2 STATIC ${LIB_EM_SOURCE})
    ADD_LIBRARY(GLEM2 STATIC ${LIB_GLEM_SOURCE})
    TARGET_LINK_LIBRARIES(EM2 FFTW3::FFTW3)
    ADD_DEFINITIONS(-DBOOST_DISABLE_ASSERTS)
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

TARGET_LINK_LIBRARIES(EM2 GSL::gsl GSL::gslcblas Freetype::Freetype)
TARGET_LINK_LIBRARIES(GLEM2 EM2)

IF(ENABLE_OPENGL)
	target_link_libraries(EM2 OpenGL::OpenGL)
ENDIF()

IF(ENABLE_FTGL)
	TARGET_LINK_LIBRARIES(GLEM2 FTGL::FTGL)
ENDIF()

IF(ENABLE_TIFF)
    TARGET_LINK_LIBRARIES(EM2 TIFF::TIFF JPEG::JPEG)
ENDIF()

IF(ENABLE_JPEG)
    TARGET_LINK_LIBRARIES(EM2 JPEG::JPEG)
ENDIF()

TARGET_LINK_LIBRARIES(EM2 HDF5::HDF5)
IF(WIN32)
	TARGET_LINK_LIBRARIES(EM2 SZLIB::SZLIB)
ENDIF()

IF(ENABLE_PNG)
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

IF(ENABLE_DJBFFT)
	TARGET_LINK_LIBRARIES(EM2 DJBFFT::DJBFFT)
ENDIF()

IF(ENABLE_OPTPP)
	TARGET_LINK_LIBRARIES(EM2 OPTPP::OPTPP)
ENDIF()

IF(ENABLE_ACML_FFT)
	TARGET_LINK_LIBRARIES(EM2 ACML::ACML)
ENDIF()

IF(ENABLE_OPENGL)
    INSTALL(TARGETS GLEM2
        DESTINATION ${EMAN_PY_INSTALL_PREFIX})
ENDIF()

INSTALL(TARGETS EM2
		DESTINATION ${EMAN_PY_INSTALL_PREFIX})

