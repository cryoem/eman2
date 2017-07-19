find_package(PNG REQUIRED)

message_var(PNG_LIBRARIES)

if(PNG_FOUND)
	set_target_properties(PNG::PNG
						  PROPERTIES
						  INTERFACE_COMPILE_DEFINITIONS USE_PNG
						  )
endif()
