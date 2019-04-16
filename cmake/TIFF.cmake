find_package(TIFF REQUIRED)

message_var(TIFF_LIBRARIES)

if(TIFF_FOUND)
	set_target_properties(TIFF::TIFF PROPERTIES
			INTERFACE_COMPILE_DEFINITIONS USE_TIFF
			)
endif()
