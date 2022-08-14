function(ADD_PYSTE_LIB libname)
	add_library(${libname} SHARED lib${libname}.cpp)
	
	if(WIN32)
		set(suffix .pyd)
	else()
		set(suffix .so)
	endif()

	set_target_properties(${libname}
			PROPERTIES
			PREFIX lib
			SUFFIX ${suffix}
			LINK_LIBRARIES "EM2;Boost::python${boost_py_ver}"
			)
	
	install(TARGETS ${libname} DESTINATION ${Python3_SITELIB})
endfunction()

ADD_PYSTE_LIB(pyAligner2)
ADD_PYSTE_LIB(pyAverager2)
ADD_PYSTE_LIB(pyCmp2)
ADD_PYSTE_LIB(pyProcessor2)
ADD_PYSTE_LIB(pyReconstructor2) 
ADD_PYSTE_LIB(pyProjector2)
ADD_PYSTE_LIB(pyEMObject2) 
ADD_PYSTE_LIB(pyEMData2)
target_sources(pyEMData2 PRIVATE emdata_wrapitems.cpp)
ADD_PYSTE_LIB(pyGeometry2)
ADD_PYSTE_LIB(pyTransform2)
ADD_PYSTE_LIB(pyPointArray2)
ADD_PYSTE_LIB(pyUtils2)
ADD_PYSTE_LIB(pyBoxingTools2)
ADD_PYSTE_LIB(pyFundamentals2)
ADD_PYSTE_LIB(pyTypeConverter2)
target_sources(pyTypeConverter2 PRIVATE typeconverter.cpp)
ADD_PYSTE_LIB(pyPolarData2)
ADD_PYSTE_LIB(pyAnalyzer2)
ADD_PYSTE_LIB(pyPDBReader2)

if(NOT WIN32)
	ADD_PYSTE_LIB(pyTomoSeg2)
endif()

if(ENABLE_OPENGL)
	ADD_PYSTE_LIB(pyGLUtils2)
	ADD_PYSTE_LIB(pyMarchingCubes2)
	
	target_sources(pyMarchingCubes2 PRIVATE
				   ${CMAKE_SOURCE_DIR}/libEM/marchingcubes.cpp
				   )
	target_sources(pyGLUtils2 PRIVATE
				   ${CMAKE_SOURCE_DIR}/libEM/emftgl.cpp
				   ${CMAKE_SOURCE_DIR}/libEM/glutil.cpp
				   ${CMAKE_SOURCE_DIR}/libEM/marchingcubes.cpp
				   )
	target_link_libraries(pyGLUtils2       EMAN::OpenGL FTGL::FTGL)
	target_link_libraries(pyMarchingCubes2 EMAN::OpenGL)
endif()

file(GLOB empythonlibs "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
install(FILES ${empythonlibs}
		DESTINATION ${Python3_SITELIB}
		COMPONENT PythonFiles
		)

add_subdirectory(qtgui)

target_link_libraries(pyUtils2         Python3::NumPy Boost::numpy${boost_py_ver})
target_link_libraries(pyTypeConverter2 Python3::NumPy Boost::numpy${boost_py_ver})
