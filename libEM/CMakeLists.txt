SET(LIB_EM_SOURCE byteorder.C
		  emobject.C
		  log.C
		  imageio.C
		  util.C
		  emutil.C
		  floatstat.C
		  transform.C
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
		  emcache.C )
		  
ADD_LIBRARY(EM SHARED ${LIB_EM_SOURCE})

IF (ENABLE_TIF)
	TARGET_LINK_LIBRARIES(EM  ${TIF_LIBRARY})
ENDIF(ENABLE_TIF)

IF (ENABLE_HDF)
	TARGET_LINK_LIBRARIES(EM  ${HDF_LIBRARY})
ENDIF(ENABLE_HDF)

IF (ENABLE_PNG)
	TARGET_LINK_LIBRARIES(EM  ${PNG_LIBRARY})
ENDIF(ENABLE_PNG)

INSTALL_TARGETS(/lib EM)
INSTALL_FILES(/include "^.+\\.h$")
