SET(LIB_EM_SOURCE byteorder.cpp
    emobject.cpp
    emfft.cpp
    log.cpp
    imageio.cpp
    util.cpp
    emutil.cpp
    interp.cpp
    quaternion.cpp
    transform.cpp
    geometry.cpp
    pointarray.cpp
    mrcio.cpp
    dm3io.cpp
    tifio.cpp
    hdfio.cpp
    emdata.cpp
    pifio.cpp
    vtkio.cpp
    spiderio.cpp
    sspiderio.cpp
    pgmio.cpp
    emimio.cpp
    icosio.cpp
    lstio.cpp
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
    exception.cpp
    testutil.cpp
    fundamentals.cpp
    normpadft.cpp
    periodogram.cpp
    fourierfilter.cpp
    rsconvolution.cpp
    tflm.cpp
    )

IF(WIN32)
    ADD_LIBRARY(EM2 STATIC ${LIB_EM_SOURCE})
    TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
ELSE(WIN32)
    ADD_LIBRARY(EM2 SHARED ${LIB_EM_SOURCE})
    TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m)
ENDIF(WIN32)


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
ENDIF(ENABLE_PNG)

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

INSTALL_TARGETS(/lib EM2)
INSTALL_FILES(/include "^.+\\.h$")
