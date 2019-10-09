find_package(Python REQUIRED)
find_package(NumPy  REQUIRED)

find_package(Nosetests)

# Find Boost
include(${CMAKE_SOURCE_DIR}/cmake/Boost.cmake)

include(${CMAKE_SOURCE_DIR}/cmake/GSL.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/HDF5.cmake)
