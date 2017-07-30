CHECK_REQUIRED_LIB(FTGL ftgl FTGL/FTGL.h ftgl FTGL/ftgl.h)

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
	if(WIN32)
		if(ENABLE_STATIC_FTGL)
			target_compile_definitions(FTGL INTERFACE FTGL_LIBRARY_STATIC)
		else()
			target_compile_definitions(FTGL INTERFACE FTGL_LIBRARY)
		endif()
	endif()
	
	include(${CMAKE_SOURCE_DIR}/cmake/Freetype.cmake)
	target_link_libraries(FTGL INTERFACE Freetype::Freetype)
endif()
