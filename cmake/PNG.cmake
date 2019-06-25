find_package(PNG REQUIRED)

cmake_print_variables(PNG_LIBRARIES)

if(PNG_FOUND)
	set_target_properties(PNG::PNG
						  PROPERTIES
						  INTERFACE_COMPILE_DEFINITIONS USE_PNG
						  )
endif()
