CHECK_REQUIRED_LIB(FTGL ftgl FTGL/FTGL.h ftgl FTGL/ftgl.h)

if(WIN32)
	add_definitions(-DFTGL_LIBRARY_STATIC)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FTGL
								  REQUIRED_VARS
									  FTGL_INCLUDE_PATH
									  FTGL_LIBRARY
								  )

if(FTGL_FOUND AND NOT TARGET FTGL)
	add_library(FTGL INTERFACE)
	add_library(FTGL::FTGL ALIAS FTGL)
	set_target_properties(FTGL
						  PROPERTIES
						  INTERFACE_INCLUDE_DIRECTORIES "${FTGL_INCLUDE_PATH}"
						  INTERFACE_LINK_LIBRARIES      "${FTGL_LIBRARY}"
						  INTERFACE_COMPILE_DEFINITIONS  USE_FTGL
						  )
endif()
