if(WIN32)
	set(HDF5_LIBRARY ${EMAN_PREFIX_LIB}/libhdf5.lib CACHE PATH "")
endif()

CHECK_REQUIRED_LIB(HDF5 hdf5 hdf5.h "" "")

cmake_print_variables(HDF5_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5
		REQUIRED_VARS HDF5_LIBRARY
		)

if(HDF5_FOUND AND NOT TARGET HDF5)
	add_library(HDF5 INTERFACE)
	add_library(HDF5::HDF5 ALIAS HDF5)
	set_target_properties(HDF5 PROPERTIES 
			INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_PATH}"
			INTERFACE_LINK_LIBRARIES      "${HDF5_LIBRARY}"
			INTERFACE_COMPILE_DEFINITIONS USE_HDF5
			)
endif()
