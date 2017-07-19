find_package(PNG REQUIRED)

message_var(PNG_LIBRARIES)

if(PNG_FOUND AND NOT TARGET PNG)
	add_library(PNG INTERFACE)
	set_target_properties(PNG
						  PROPERTIES
						  INTERFACE_LINK_LIBRARIES      "${PNG_LIBRARIES}"
						  INTERFACE_COMPILE_DEFINITIONS USE_PNG
						  )
endif()
