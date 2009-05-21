SUBDIRS(sparx)

SET(LIB_EM_SOURCE byteorder.cpp
    emftgl.cpp
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
    mrcio.cpp
    dm3io.cpp
    tifio.cpp
    hdfio.cpp
    hdfio2.cpp
    jpegio.cpp
    emdata.cpp
    emdata_io.cpp
    emdata_core.cpp
    emdata_cuda.cpp
    glutil.cpp
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
    emimio.cpp
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
    marchingcubes.cpp
    randnum.cpp
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
    )


IF(WIN32)
    ADD_LIBRARY(EM2 STATIC ${LIB_EM_SOURCE})
    TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} )
    ADD_DEFINITIONS(-DBOOST_DISABLE_ASSERTS)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
ELSE(WIN32)
    ADD_LIBRARY(EM2 SHARED ${LIB_EM_SOURCE})
    IF(ENABLE_FFTW2 OR ENABLE_FFTW3)
    	TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m)
    ELSE(ENABLE_FFTW2 OR ENABLE_FFTW3)
    	TARGET_LINK_LIBRARIES(EM2 ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m)
    ENDIF(ENABLE_FFTW2 OR ENABLE_FFTW3)
    
  
#    TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m)
#	TARGET_LINK_LIBRARIES(EM2 ${ACML_LIBRARY} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m g2c)
ENDIF(WIN32)

IF(ENABLE_OPENGL)
	TARGET_LINK_LIBRARIES(EM2 ${GL_LIBRARY} ${GLU_LIBRARY} )
ENDIF(ENABLE_OPENGL)

IF(ENABLE_FTGL)
	INCLUDE_DIRECTORIES(${FTGL_INCLUDE_PATH})
	IF(ENABLE_STATIC_FTGL)
		IF(WIN32)
			TARGET_LINK_LIBRARIES(EM2 ${FTGL_LIBRARY} $(FREETYPE_LIBRARY))
		ELSE(WIN32)
			TARGET_LINK_LIBRARIES(EM2 ${FTGL_LIBRARY} freetype)
		ENDIF(WIN32)
	ELSE(ENABLE_STATIC_FTGL)
		TARGET_LINK_LIBRARIES(EM2 ${FTGL_LIBRARY})
	ENDIF(ENABLE_STATIC_FTGL)
ENDIF(ENABLE_FTGL)

IF (ENABLE_TIFF)
    TARGET_LINK_LIBRARIES(EM2 ${TIFF_LIBRARY})
    TARGET_LINK_LIBRARIES(EM2 ${JPEG_LIBRARY})
ENDIF(ENABLE_TIFF)

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

IF (ENABLE_CUDA)
    TARGET_LINK_LIBRARIES(EM2 EM2Cuda cudart)
    TARGET_LINK_LIBRARIES(EM2 EM2SparxCuda)
ENDIF(ENABLE_CUDA)

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

#IF(ENABLE_CUDA_FFT)
#	TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES})
#ENDIF(ENABLE_CUDA_FFT)

IF(ENABLE_OPTPP)
	TARGET_LINK_LIBRARIES(EM2 ${OPTPP_LIBRARY} ${NEWMAT_LIBRARY})
ENDIF(ENABLE_OPTPP)

IF(ENABLE_ACML_FFT)
	TARGET_LINK_LIBRARIES(EM2 ${ACML_LIBRARY} ${G2C_LIBRARY})
ENDIF(ENABLE_ACML_FFT)

INSTALL_TARGETS(/lib EM2)
INSTALL_FILES(/include "^.+\\.h$")
