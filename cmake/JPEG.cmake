find_package(JPEG REQUIRED)

message_var(JPEG_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JPEG
								  REQUIRED_VARS JPEG_LIBRARIES
								  )

if(JPEG_FOUND AND NOT TARGET JPEG)
	add_library(JPEG INTERFACE)
	add_library(JPEG::JPEG ALIAS JPEG)
	set_target_properties(JPEG PROPERTIES
						  INTERFACE_LINK_LIBRARIES      "${JPEG_LIBRARIES}"
						  INTERFACE_COMPILE_DEFINITIONS USE_JPEG
						  )
endif()
