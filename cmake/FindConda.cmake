# Assuming the active conda environment is on PATH, this finds the path of bin/ in the environment
find_program(CONDA_EXECUTABLE conda
		PATHS ${CONDA_ROOT}/bin ${CONDA_ROOT}/Scripts ENV PATH
		)
if(CONDA_EXECUTABLE)
	execute_process(COMMAND ${CONDA_EXECUTABLE} info --root
			OUTPUT_VARIABLE out_var
			OUTPUT_STRIP_TRAILING_WHITESPACE
			)
	execute_process(COMMAND ${CONDA_EXECUTABLE} info --envs
			)
else()
	message("\nNo conda environment found in PATH!\nPATH=$ENV{PATH}\n")
endif()

set_cache_var_to_var(CONDA_ROOT out_var)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Conda
		REQUIRED_VARS CONDA_ROOT
		)
