SET(LIB_EM_SOURCE byteorder.cpp
    emobject.cpp
    emfft.cpp
    log.cpp
    imageio.cpp
    util.cpp
    emutil.cpp
    interp.cpp
    transform.cpp
    geometry.cpp
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
    filter.cpp
    aligner.cpp
    projector.cpp
    cmp.cpp
    averager.cpp
    reconstructor.cpp
    )

ADD_LIBRARY(EM2 SHARED ${LIB_EM_SOURCE})

TARGET_LINK_LIBRARIES(EM2 ${FFTW_LIBRARIES} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m)

IF (ENABLE_TIFF)
    TARGET_LINK_LIBRARIES(EM2 ${TIFF_LIBRARY})
ENDIF(ENABLE_TIFF)

IF (ENABLE_HDF)
    TARGET_LINK_LIBRARIES(EM2 ${HDF5_LIBRARY})
ENDIF(ENABLE_HDF)

IF (ENABLE_PNG)
    TARGET_LINK_LIBRARIES(EM2 ${PNG_LIBRARY})
    TARGET_LINK_LIBRARIES(EM2 z)
ENDIF(ENABLE_PNG)

INSTALL_TARGETS(/lib EM2)
INSTALL_FILES(/include "^.+\\.h$")
