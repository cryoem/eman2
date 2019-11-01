find_package(JPEG REQUIRED)

cmake_print_variables(JPEG_LIBRARIES)

if(JPEG_FOUND)
	set_target_properties(JPEG::JPEG PROPERTIES
						  INTERFACE_COMPILE_DEFINITIONS USE_JPEG
						  )
endif()
