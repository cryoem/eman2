execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
		OUTPUT_VARIABLE numpy_include_dir
		OUTPUT_STRIP_TRAILING_WHITESPACE
		)

set_cache_var_to_var(NUMPY_INCLUDE_DIR numpy_include_dir)

cmake_print_variables(NUMPY_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NumPy
		REQUIRED_VARS NUMPY_INCLUDE_DIR
		)

if(NumPy_FOUND AND NOT TARGET NumPy)
	add_library(NumPy INTERFACE)
	set_target_properties(NumPy
			PROPERTIES
			INTERFACE_INCLUDE_DIRECTORIES ${NUMPY_INCLUDE_DIR}
			)
	
endif()
