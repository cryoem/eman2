CHECK_LIB_ONLY(SZLIB libszip)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SZLIB
								  REQUIRED_VARS SZLIB_LIBRARY
								  )

if(SZLIB_FOUND AND NOT TARGET SZLIB)
	add_library(SZLIB INTERFACE)
	add_library(SZLIB::SZLIB ALIAS SZLIB)
	set_target_properties(SZLIB
						  PROPERTIES
						  INTERFACE_LINK_LIBRARIES ${SZLIB_LIBRARY}
						  )
endif()
