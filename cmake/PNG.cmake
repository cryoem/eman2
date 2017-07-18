find_package(PNG REQUIRED)

message_var(PNG_LIBRARIES)

ADD_DEFINITIONS(-DUSE_PNG)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PNG
								  REQUIRED_VARS
								  PNG_LIBRARIES
								  )

if(PNG_FOUND AND NOT TARGET PNG)
	add_library(PNG INTERFACE)
	add_library(PNG::PNG ALIAS PNG)
	set_target_properties(FTGL
	set_target_properties(PNG
						  PROPERTIES
						  INTERFACE_LINK_LIBRARIES      "${PNG_LIBRARIES}"
						  )
endif()
