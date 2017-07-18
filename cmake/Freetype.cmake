find_package(Freetype REQUIRED)

message("FREETYPE_LIBRARIES:    ${FREETYPE_LIBRARIES}")
message("FREETYPE_INCLUDE_DIRS: ${FREETYPE_INCLUDE_DIRS}")

if(Freetype_FOUND AND NOT TARGET Freetype)
	add_library(Freetype INTERFACE)
	add_library(Freetype::Freetype ALIAS Freetype)
	set_target_properties(Freetype
						  PROPERTIES
						  INTERFACE_INCLUDE_DIRECTORIES "${FREETYPE_INCLUDE_DIRS}"
						  INTERFACE_LINK_LIBRARIES      "${FREETYPE_LIBRARIES}"
						  )
endif()
