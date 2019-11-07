# Use the system compiler to get the path to OS libraries in order to find system OpenGL
if(CMAKE_SYSTEM_NAME MATCHES "Linux" AND "$ENV{CONDA_BUILD}" STREQUAL "1")
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

set(OpenGL_GL_PREFERENCE GLVND)

find_package(OpenGL REQUIRED COMPONENTS OpenGL)

if(CMAKE_SYSTEM_NAME MATCHES "Linux" AND "$ENV{CONDA_BUILD}" STREQUAL "1")
	set(CMAKE_LIBRARY_ARCHITECTURE ${CMAKE_LIBRARY_ARCHITECTURE_BACKUP})
endif()

cmake_print_variables(OPENGL_INCLUDE_DIR)
cmake_print_variables(OPENGL_LIBRARIES)
cmake_print_variables(OPENGL_glu_LIBRARY)
cmake_print_variables(OPENGL_opengl_LIBRARY)
cmake_print_variables(OPENGL_gl_LIBRARY)
cmake_print_variables(OPENGL_FOUND)
cmake_print_variables(OPENGL_GLU_FOUND)
cmake_print_variables(OpenGL_OpenGL_FOUND)
cmake_print_variables(OPENGL_opengl_LIBRARY)
cmake_print_variables(OPENGL_gl_LIBRARY)

if(OpenGL_FOUND AND NOT TARGET OpenGL AND NOT TARGET EMAN::OpenGL)
	add_library(OpenGL INTERFACE)
	add_library(EMAN::OpenGL ALIAS OpenGL)

	is_target_exist(OpenGL)
	is_target_exist(OpenGL::OpenGL)
	is_target_exist(OpenGL::GL)
	is_target_exist(OpenGL::GLU)

	set_target_properties(OpenGL PROPERTIES
						  INTERFACE_COMPILE_DEFINITIONS USE_OPENGL
						  )
	target_link_libraries(OpenGL INTERFACE OpenGL::GL OpenGL::GLU)
		
	# Symlink to GL. When the value is /usr/include, cmake ignores it.
	# So, this is a workaround to include OpenGL headers
	if(CMAKE_SYSTEM_NAME MATCHES "Linux" AND "$ENV{CONDA_BUILD}" STREQUAL "1")
		execute_process(
				COMMAND ${CMAKE_COMMAND} -E create_symlink ${OPENGL_INCLUDE_DIR}/GL ${CMAKE_CURRENT_BINARY_DIR}/GL
		)
	endif()
endif()
