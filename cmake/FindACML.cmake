CHECK_REQUIRED_LIB(ACML acml acml.h "" "")
FIND_LIBRARY(G2C_LIBRARY NAMES g2c PATHS
			 /usr/lib64
			 /usr/lib
			 /usr/local/lib
			 $ENV{HOME}/lib
			 )

cmake_print_variables(ACML_INCLUDE_PATH)
cmake_print_variables(ACML_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ACML
								  REQUIRED_VARS
								  ACML_INCLUDE_PATH
								  ACML_LIBRARY
								  )

if(ACML_FOUND AND NOT TARGET ACML)
	add_library(ACML INTERFACE)
	add_library(ACML::ACML ALIAS ACML)
	set_target_properties(ACML PROPERTIES
						  INTERFACE_INCLUDE_DIRECTORIES "${ACML_INCLUDE_PATH}"
						  INTERFACE_LINK_LIBRARIES      "${ACML_LIBRARY};${G2C_LIBRARY}"
						  INTERFACE_COMPILE_DEFINITIONS USE_ACML
						  )
endif()
