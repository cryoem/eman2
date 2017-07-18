find_package(TIFF REQUIRED)
add_definitions(-DUSE_TIFF)
CHECK_LIB_ONLY(JPEG jpeg)

message_var(TIFF_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TIFF
		REQUIRED_VARS TIFF_LIBRARIES
		)
