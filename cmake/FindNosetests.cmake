find_program(NOSETESTS_EXECUTABLE nosetests
		PATHS ${CONDA_PREFIX}/bin ${CONDA_PREFIX}/Scripts ENV PATH
		)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Nosetests
		REQUIRED_VARS NOSETESTS_EXECUTABLE
		)
