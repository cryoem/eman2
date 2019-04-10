# Use the system compiler to get the path to OS libraries in order to find system OpenGL
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
	message("CMAKE_LIBRARY_ARCHITECTURE: ${CMAKE_LIBRARY_ARCHITECTURE}")
	set(CMAKE_LIBRARY_ARCHITECTURE_BACKUP ${CMAKE_LIBRARY_ARCHITECTURE})
	execute_process(COMMAND which g++
					OUTPUT_VARIABLE SYSTEM_COMPILER
					OUTPUT_STRIP_TRAILING_WHITESPACE
					)
	execute_process(COMMAND ${SYSTEM_COMPILER} -print-multiarch
					OUTPUT_VARIABLE CMAKE_LIBRARY_ARCHITECTURE
					OUTPUT_STRIP_TRAILING_WHITESPACE
					)
endif()

find_package(OpenGL REQUIRED)

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
	set(CMAKE_LIBRARY_ARCHITECTURE ${CMAKE_LIBRARY_ARCHITECTURE_BACKUP})
endif()

message_var(OPENGL_INCLUDE_DIR)
message_var(OPENGL_LIBRARIES)

if(OpenGL_FOUND AND NOT TARGET OpenGL)
	add_library(OpenGL INTERFACE)
	add_library(OpenGL::OpenGL ALIAS OpenGL)
	set_target_properties(OpenGL PROPERTIES
						  INTERFACE_COMPILE_DEFINITIONS USE_OPENGL
						  )
	target_link_libraries(OpenGL INTERFACE OpenGL::GL OpenGL::GLU)
		
	# Symlink to GL. When the value is /usr/include, cmake ignores it.
	# So, this is a workaround to include OpenGL headers
	if(CMAKE_SYSTEM_NAME MATCHES "Linux")
		execute_process(
				COMMAND ${CMAKE_COMMAND} -E create_symlink ${OPENGL_INCLUDE_DIR}/GL ${CMAKE_CURRENT_BINARY_DIR}/GL
		)
	endif()
endif()
