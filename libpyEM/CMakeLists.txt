SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

INCLUDE_DIRECTORIES(${EMAN_SOURCE_DIR}/libpyEM ${BOOST_INCLUDE_PATH} ${PYTHON_INCLUDE_PATH})

IF(WIN32)
    LINK_LIBRARIES(${BOOST_LIBRARY} ${PYTHON_LIBRARY} EM2 EMext)
ELSE(WIN32)
    LINK_LIBRARIES(${BOOST_LIBRARY} EM2 EMext)
ENDIF(WIN32)

IF(CMAKE_SYSTEM MATCHES "IRIX.*")
    INCLUDE_DIRECTORIES(${BOOST_INCLUDE_PATH}/boost/compatibility/cpp_c_headers)
ENDIF(CMAKE_SYSTEM MATCHES "IRIX.*")

IF(WIN32)
    ADD_LIBRARY(libpyGeometry2 SHARED libpyGeometry2.cpp)
    ADD_LIBRARY(libpyFactory2 SHARED libpyFactory2.cpp)
    ADD_LIBRARY(libpyUtils2 SHARED libpyUtils2.cpp)
    ADD_LIBRARY(libpyEMData2 SHARED libpyEMData2.cpp)
    INSTALL_TARGETS(/lib libpyGeometry2 libpyFactory2 libpyUtils2 libpyEMData2)
ELSE(WIN32)
    ADD_LIBRARY(pyGeometry2 SHARED libpyGeometry2.cpp)
    ADD_LIBRARY(pyFactory2 SHARED libpyFactory2.cpp)
    ADD_LIBRARY(pyUtils2 SHARED libpyUtils2.cpp pylist.cpp)
    ADD_LIBRARY(pyEMData2 SHARED libpyEMData2.cpp)
    INSTALL_TARGETS(/lib pyGeometry2 pyFactory2 pyUtils2 pyEMData2)
ENDIF(WIN32)


INSTALL_FILES(/lib "^.+\\.py$")
INSTALL_FILES(/include "^.+\\.h$")

