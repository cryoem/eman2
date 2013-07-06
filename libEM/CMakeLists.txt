ADD_SUBDIRECTORY(sparx)
ADD_SUBDIRECTORY(gorgon)
ADD_SUBDIRECTORY(plugins)
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
    )

SET(LIB_GLEM_SOURCE 
		glutil.cpp
		marchingcubes.cpp
		emftgl.cpp)

IF(WIN32)
    ADD_LIBRARY(EM2 STATIC ${LIB_EM_SOURCE})
    ADD_LIBRARY(GLEM2 STATIC ${LIB_GLEM_SOURCE})
    TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} )
    TARGET_LINK_LIBRARIES(GLEM2 EM2)
    ADD_DEFINITIONS(-DBOOST_DISABLE_ASSERTS)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
ELSE(WIN32)
    ADD_LIBRARY(EM2 SHARED ${LIB_EM_SOURCE})
    ADD_LIBRARY(GLEM2 SHARED ${LIB_GLEM_SOURCE})
    IF(ENABLE_FFTW2 OR ENABLE_FFTW3)
    	TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m)
    ELSE(ENABLE_FFTW2 OR ENABLE_FFTW3)
    	TARGET_LINK_LIBRARIES(EM2 ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m)
    ENDIF(ENABLE_FFTW2 OR ENABLE_FFTW3)    
    TARGET_LINK_LIBRARIES(GLEM2 EM2 )  
#    TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m)
#	TARGET_LINK_LIBRARIES(EM2 ${ACML_LIBRARY} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m g2c)
ENDIF(WIN32)


IF(ENABLE_OPENGL)
	TARGET_LINK_LIBRARIES(GLEM2 ${GL_LIBRARY} ${GLU_LIBRARY} )
ENDIF(ENABLE_OPENGL)

IF(ENABLE_FTGL)
	INCLUDE_DIRECTORIES(${FTGL_INCLUDE_PATH})
	IF(ENABLE_STATIC_FTGL)
		IF(WIN32)
			TARGET_LINK_LIBRARIES(GLEM2 ${FTGL_LIBRARY} $(FREETYPE_LIBRARY))
		ELSE(WIN32)
			TARGET_LINK_LIBRARIES(GLEM2 ${FTGL_LIBRARY} freetype)
		ENDIF(WIN32)
	ELSE(ENABLE_STATIC_FTGL)
		TARGET_LINK_LIBRARIES(GLEM2 ${FTGL_LIBRARY})
	ENDIF(ENABLE_STATIC_FTGL)
ENDIF(ENABLE_FTGL)

IF (ENABLE_TIFF)
    TARGET_LINK_LIBRARIES(EM2 ${TIFF_LIBRARY})
    TARGET_LINK_LIBRARIES(EM2 ${JPEG_LIBRARY})
ENDIF(ENABLE_TIFF)

IF (ENABLE_JPEG)
    TARGET_LINK_LIBRARIES(EM2 ${JPEG_LIBRARY})
ENDIF(ENABLE_JPEG)

IF (ENABLE_HDF5)
    TARGET_LINK_LIBRARIES(EM2 ${HDF5_LIBRARY})
    TARGET_LINK_LIBRARIES(EM2 ${ZLIB_LIBRARY})
	IF (WIN32)
		TARGET_LINK_LIBRARIES(EM2 ${SZLIB_LIBRARY})
	ENDIF(WIN32)
ENDIF(ENABLE_HDF5)

IF (ENABLE_PNG)
    TARGET_LINK_LIBRARIES(EM2 ${PNG_LIBRARY})
    TARGET_LINK_LIBRARIES(EM2 ${ZLIB_LIBRARY})
ENDIF(ENABLE_PNG)

IF (ENABLE_SPARX_CUDA)
    TARGET_LINK_LIBRARIES(EM2 EM2SparxCuda)
    TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES})
ENDIF(ENABLE_SPARX_CUDA)

IF (ENABLE_EMAN_CUDA)
    TARGET_LINK_LIBRARIES(EM2 EM2Cuda)
    TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES})
ENDIF (ENABLE_EMAN_CUDA)

IF(ENABLE_NFFT)
	SET(LIB_EM_SOURCE ${LIB_EM_SOURCE} pointarray.cpp)
	TARGET_LINK_LIBRARIES(EM2 ${NFFT_LIBRARY} ${FFTW2D_LIBRARY})
ENDIF(ENABLE_NFFT)

IF(ENABLE_NFFT2)
	SET(LIB_EM_SOURCE ${LIB_EM_SOURCE} pointarray.cpp)
	TARGET_LINK_LIBRARIES(EM2 ${NFFT2_LIBRARY} ${FFTW3D_LIBRARY})
ENDIF(ENABLE_NFFT2)

IF(ENABLE_DJBFFT)
	TARGET_LINK_LIBRARIES(EM2 ${DJBFFT_LIBRARY})
ENDIF(ENABLE_DJBFFT)

IF(ENABLE_OPTPP)
	TARGET_LINK_LIBRARIES(EM2 ${OPTPP_LIBRARY} ${NEWMAT_LIBRARY})
ENDIF(ENABLE_OPTPP)

IF(ENABLE_ACML_FFT)
	TARGET_LINK_LIBRARIES(EM2 ${ACML_LIBRARY} ${G2C_LIBRARY})
ENDIF(ENABLE_ACML_FFT)

IF(ENABLE_OPENGL)
    INSTALL(TARGETS EM2 GLEM2 
        DESTINATION lib)
ELSE(ENABLE_OPENGL)
    INSTALL(TARGETS EM2 
        DESTINATION lib)
ENDIF(ENABLE_OPENGL)

# if(ENABLE_HDFIOCACHE)
#    TARGET_LINK_LIBRARIES(EM2 ${BOOST_THREAD_LIBRARY} ${BOOST_FILESYSTEM_LIBRARY} ${BOOST_SYSTEM_LIBRARY})
# endif(ENABLE_HDFIOCACHE)

FILE(GLOB headerfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${headerfiles}
  DESTINATION    include
)
