if(WIN32)
	set(HDF5_LIBRARY ${EMAN_PREFIX_LIB}/libhdf5.lib CACHE PATH "")
endif()

CHECK_OPTIONAL_LIB(HDF5 hdf5 hdf5.h)

message_var(HDF5_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5
		REQUIRED_VARS HDF5_LIBRARY
		)
