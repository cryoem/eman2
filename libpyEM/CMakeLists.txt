
IF(WIN32)	
    LINK_LIBRARIES(${BOOST_LIBRARY} ${PYTHON_LIBRARY} EM2 EMext)
ELSE(WIN32)
    SET(CMAKE_CXX_FLAGS_DEBUG "")
    SET(CMAKE_CXX_FLAGS "${CMAKE_ORIG_FLAGS} ${EMAN_CXX_FLAGS} -O0")
    SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -O0")
    ## On linux, linking agains libpython24.a fails, and libpython24.so is
    ## not available in the default cctbx configuration.
    ## 
    ## LINK_LIBRARIES(${BOOST_LIBRARY} ${PYTHON_LIBRARY} EM2 EMext)
    LINK_LIBRARIES(${BOOST_LIBRARY}  EM2 EMext)
ENDIF(WIN32)

#this definition is for boost.python 1.35.0 
ADD_DEFINITIONS(-DBOOST_PYTHON_NO_PY_SIGNATURES)

IF(ENABLE_FTGL)
	INCLUDE_DIRECTORIES(${FTGL_INCLUDE_PATH})
	IF(ENABLE_STATIC_FTGL)
		IF(WIN32)
			LINK_LIBRARIES(${FTGL_LIBRARY} $(FREETYPE_LIBRARY))
		ELSE(WIN32)
			LINK_LIBRARIES(${FTGL_LIBRARY} freetype)
		ENDIF(WIN32)
	ELSE(ENABLE_STATIC_FTGL)
		LINK_LIBRARIES(${FTGL_LIBRARY})
	ENDIF(ENABLE_STATIC_FTGL)
ENDIF(ENABLE_FTGL)

IF(APPLE)	
    LINK_LIBRARIES(${BOOST_LIBRARY} ${PYTHON_LIBRARY} EM2 EMext)
ENDIF(APPLE)

INCLUDE_DIRECTORIES(${EMAN_SOURCE_DIR}/libpyEM ${BOOST_INCLUDE_PATH} ${PYTHON_INCLUDE_PATH})

IF(CMAKE_SYSTEM MATCHES "IRIX.*")
    INCLUDE_DIRECTORIES(${BOOST_INCLUDE_PATH}/boost/compatibility/cpp_c_headers)
ENDIF(CMAKE_SYSTEM MATCHES "IRIX.*")


MACRO(ADD_PYSTE_LIB_WITH_SRC libname srcfile)
	IF(WIN32)
		ADD_LIBRARY(lib${libname} SHARED lib${libname}.cpp ${srcfile})
		#INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib lib${libname})
		INSTALL(FILES ${LIBRARY_OUTPUT_PATH}/release/lib${libname}.dll
				DESTINATION	lib
				RENAME lib${libname}.pyd)
	ELSE(WIN32)
		ADD_LIBRARY(${libname} SHARED lib${libname}.cpp ${srcfile})
		INSTALL_TARGETS(/lib ${libname})
	ENDIF(WIN32)
ENDMACRO(ADD_PYSTE_LIB_WITH_SRC)


MACRO(ADD_PYSTE_LIB libname)
	ADD_PYSTE_LIB_WITH_SRC(${libname} "")
ENDMACRO(ADD_PYSTE_LIB)

ADD_PYSTE_LIB(pyAligner2)
ADD_PYSTE_LIB(pyAverager2)
ADD_PYSTE_LIB(pyCmp2)
ADD_PYSTE_LIB(pyProcessor2)
ADD_PYSTE_LIB(pyReconstructor2) 
ADD_PYSTE_LIB(pyProjector2)
ADD_PYSTE_LIB(pyEMObject2) 
ADD_PYSTE_LIB_WITH_SRC(pyEMData2 emdata_wrapitems.cpp)
ADD_PYSTE_LIB(pyGeometry2)
ADD_PYSTE_LIB(pyTransform2)
ADD_PYSTE_LIB(pyPointArray2)
#ADD_PYSTE_LIB_WITH_SRC(pyUtils2 sparx_apmq)
ADD_PYSTE_LIB(pyUtils2)
ADD_PYSTE_LIB(pyBoxingTools2)
ADD_PYSTE_LIB(pyFundamentals2)
ADD_PYSTE_LIB_WITH_SRC(pyTypeConverter2 typeconverter.cpp)
ADD_PYSTE_LIB(pyPolarData2)
ADD_PYSTE_LIB(pyAnalyzer2)
ADD_PYSTE_LIB(pyMarchingCubes2)

INSTALL_FILES(/lib "^.+\\.py$")
INSTALL_FILES(/include "^.+\\.h$")
#INSTALL_FILES(/lib FILES gui/emimage.py gui/imagecontrol.py gui/emimagecanvas.py)
INSTALL_FILES(/lib FILES qtgui/emplot2d.py qtgui/emimage.py qtgui/emimage2d.py qtgui/emimagemx.py qtgui/emimage3d.py qtgui/emimageutil.py qtgui/valslider.py qtgui/pyshed.py qtgui/emshape.py)
INSTALL_FILES(/lib FILES qtgui/emimage3diso.py qtgui/emglobjects.py qtgui/emimage3dslice.py qtgui/emimage3dvol.py qtgui/emimage3dsym.py qtgui/emglplot.py qtgui/emimagemxrotor.py qtgui/emfloatingwidgets.py qtgui/emimagerotor.py qtgui/emrotor.py qtgui/emapplication.py qtgui/emanimationutil.py qtgui/emborderdecoration.py qtgui/embrowse.py  qtgui/emselector.py qtgui/emdatastorage.py qtgui/emform.py qtgui/emsprworkflow.py qtgui/emlights.py qtgui/emplot3d.py qtgui/emsave.py)
INSTALL_FILES(/lib FILES wxgui/wxemimage.py wxgui/wximagecontrol.py wxgui/wxValSlider.py)
INSTALL_PROGRAMS(/bin qtgui/emselector.py)
