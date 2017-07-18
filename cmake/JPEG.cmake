find_package(JPEG REQUIRED)
add_definitions(-DUSE_JPEG)

message_var(JPEG_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JPEG
								  REQUIRED_VARS JPEG_LIBRARIES
								  )
