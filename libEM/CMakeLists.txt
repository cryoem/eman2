SET(LIB_EM_SOURCE byteorder.C
		  emobject.C
		  emfft.C
		  filter.C
		  aligner.C
		  projector.C
		  cmp.C
		  averager.C
		  log.C
		  imageio.C
		  util.C
		  emutil.C
		  floatstat.C
		  interp.C
		  reconstructor.C
		  transform.C
		  geometry.C
		  mrcio.C
		  dm3io.C
		  tifio.C
		  hdfio.C
		  emdata.C
		  pifio.C
		  vtkio.C
		  spiderio.C
		  sspiderio.C
		  pgmio.C
		  emimio.C
		  icosio.C
		  lstio.C
		  pngio.C
		  salio.C
		  amiraio.C
		  xplorio.C
		  gatan2io.C
		  emio.C
		  imagicio.C
		  emcache.C
		  ctf.C
		  xydata.C
		  )
		  
ADD_LIBRARY(EM SHARED ${LIB_EM_SOURCE})

TARGET_LINK_LIBRARIES(EM ${RFFTW_LIBRARY} ${FFTW_LIBRARY} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} m)

IF (ENABLE_TIFF)
	TARGET_LINK_LIBRARIES(EM  ${TIFF_LIBRARY})
ENDIF(ENABLE_TIFF)

IF (ENABLE_HDF)
	TARGET_LINK_LIBRARIES(EM  ${HDF_LIBRARY})
ENDIF(ENABLE_HDF)

IF (ENABLE_PNG)
	TARGET_LINK_LIBRARIES(EM  ${PNG_LIBRARY})
	TARGET_LINK_LIBRARIES(EM  z)
ENDIF(ENABLE_PNG)

INSTALL_TARGETS(/lib EM)
INSTALL_FILES(/include "^.+\\.h$")
