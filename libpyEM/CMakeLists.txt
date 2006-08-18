LINK_LIBRARIES(${BOOST_LIBRARY} ${PYTHON_LIBRARY} EM2 EMext)

IF(WIN32)	
ELSE(WIN32)
    SET(CMAKE_CXX_FLAGS_DEBUG "")
    SET(CMAKE_CXX_FLAGS "${CMAKE_ORIG_FLAGS} ${EMAN_CXX_FLAGS} -O0")
    SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -O0")
ENDIF(WIN32)


INCLUDE_DIRECTORIES(${EMAN_SOURCE_DIR}/libpyEM ${BOOST_INCLUDE_PATH} ${PYTHON_INCLUDE_PATH})

IF(CMAKE_SYSTEM MATCHES "IRIX.*")
    INCLUDE_DIRECTORIES(${BOOST_INCLUDE_PATH}/boost/compatibility/cpp_c_headers)
ENDIF(CMAKE_SYSTEM MATCHES "IRIX.*")


MACRO(ADD_PYSTE_LIB_WITH_SRC libname srcfile)
	IF(WIN32)
		ADD_LIBRARY(lib${libname} SHARED lib${libname}.cpp ${srcfile})
		INSTALL_TARGETS(/lib lib${libname})
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
ADD_PYSTE_LIB(pyUtils2)
ADD_PYSTE_LIB(pyFundamentals2)
ADD_PYSTE_LIB_WITH_SRC(pyTypeConverter2 typeconverter.cpp)
ADD_PYSTE_LIB(pyPolarData2)


INSTALL_FILES(/lib "^.+\\.py$")
INSTALL_FILES(/include "^.+\\.h$")
#INSTALL_FILES(/lib FILES gui/emimage.py gui/imagecontrol.py gui/emimagecanvas.py)
INSTALL_FILES(/lib FILES qtgui/emimage2d.py qtgui/emimagemx.py qtgui/emimageutil.py qtgui/valslider.py)
INSTALL_FILES(/lib FILES wxgui/wxemimage.py wxgui/wximagecontrol.py wxgui/wxValSlider.py)

