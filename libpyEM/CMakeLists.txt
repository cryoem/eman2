SET(CMAKE_CXX_FLAGS "${CMAKE_ORIG_FLAGS} -O0")

INCLUDE_DIRECTORIES(${BOOST_INCLUDE_PATH} ${PYTHON_INCLUDE_PATH})
LINK_LIBRARIES(${BOOST_LIBRARY} EM)

IF(CMAKE_SYSTEM MATCHES "IRIX.*")
	INCLUDE_DIRECTORIES(${BOOST_INCLUDE_PATH}/boost/compatibility/cpp_c_headers)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -LANG:std")
ENDIF(CMAKE_SYSTEM MATCHES "IRIX.*")


ADD_LIBRARY(pyGeometry SHARED 
	         libpyGeometry.cpp
)


ADD_LIBRARY(pyFactory SHARED 
	         libpyFactory.cpp
)


ADD_LIBRARY(pyUtils SHARED 
	         libpyUtils.cpp
)


ADD_LIBRARY(pyFilter SHARED 
	         libpyFilter.cpp
)


ADD_LIBRARY(pyEMData SHARED
		 pyem.C
		 libpyEMData.cpp
)

INSTALL_TARGETS(/lib pyGeometry pyFactory pyUtils pyFilter pyEMData)
INSTALL_FILES(/lib .py EMAN.py)
INSTALL_FILES(/include "^.+\\.h$")
