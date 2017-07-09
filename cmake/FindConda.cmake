# Assuming the active conda environment is on PATH, this finds the path of bin/ in the environment
find_program(CONDA_EXECUTABLE conda
		PATHS ${CONDA_ROOT}/bin ${CONDA_ROOT}/Scripts ENV PATH
		)
if(CONDA_EXECUTABLE)
	execute_process(COMMAND ${CONDA_EXECUTABLE} info --root
			OUTPUT_VARIABLE out_var
			OUTPUT_STRIP_TRAILING_WHITESPACE
			)
	if(out_var)
		set(CONDA_ROOT ${out_var}            CACHE PATH "")
		message("Found conda: ${CONDA_EXECUTABLE}")
		execute_process(COMMAND ${CONDA_EXECUTABLE} info --envs
				)
		message("Set CONDA_ROOT to one of the conda environment paths listed above.\n")
	endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Conda
		REQUIRED_VARS CONDA_ROOT
		)
