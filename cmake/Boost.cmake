set(Boost_USE_MULTITHREADED ON)
set(Boost_NO_BOOST_CMAKE ON)
set(boost_py_ver ${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})

find_package(Boost COMPONENTS python${boost_py_ver} numpy${boost_py_ver} REQUIRED)

message("Boost_VERSION: ${Boost_VERSION}")
message("Boost_LIBRARIES:   ${Boost_LIBRARIES}")
message("Boost_INCLUDE_DIR: ${Boost_INCLUDE_DIR}")

#this definition is for boost.python > 1.35.0
set_target_properties(Boost::python${boost_py_ver}
					  PROPERTIES
					  INTERFACE_COMPILE_DEFINITIONS BOOST_PYTHON_NO_PY_SIGNATURES
					  INTERFACE_LINK_LIBRARIES Python::Python
					  )
if(WIN32)
	ADD_DEFINITIONS(-DBOOST_DISABLE_ASSERTS)
endif()

IF(CMAKE_SYSTEM MATCHES "IRIX.*")
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR}/boost/compatibility/cpp_c_headers)
ENDIF()
