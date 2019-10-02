FIND_PATH(OPTPP_INCLUDE_PATH Opt.h /usr/include /usr/include/opt++
		  /usr/local/include /usr/local/include/opt++ $ENV{HOME}/include
		  $ENV{HOME}/include/opt++ NO_DEFAULT_PATH)
FIND_PATH(NEWMAT_INCLUDE_PATH newmat.h /usr/include /usr/include/opt++
		  /usr/include/newmat /usr/local/include/newmat
		  /usr/local/include /usr/local/include/opt++ $ENV{HOME}/include
		  $ENV{HOME}/include/opt++ NO_DEFAULT_PATH)
FIND_LIBRARY(OPTPP_LIBRARY NAMES opt-linux PATHS /usr/lib /usr/local/lib)
FIND_LIBRARY(NEWMAT_LIBRARY NAMES newmat-linux PATHS /usr/lib /usr/local/lib)

cmake_print_variables(OPTPP_INCLUDE_PATH)
cmake_print_variables(OPTPP_LIBRARY)
cmake_print_variables(NEWMAT_INCLUDE_PATH)
cmake_print_variables(NEWMAT_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OPTPP
								  REQUIRED_VARS
								  	OPTPP_INCLUDE_PATH
								  	OPTPP_LIBRARY
								  	NEWMAT_INCLUDE_PATH
								  	NEWMAT_LIBRARY
								  )

if(OPTPP_FOUND AND NOT TARGET OPTPP)
	add_library(OPTPP INTERFACE)
	set_target_properties(OPTPP PROPERTIES
						  INTERFACE_INCLUDE_DIRECTORIES "${OPTPP_INCLUDE_PATH};${NEWMAT_INCLUDE_PATH"
						  INTERFACE_LINK_LIBRARIES      "${OPTPP_LIBRARY};${NEWMAT_LIBRARY}"
						  INTERFACE_COMPILE_DEFINITIONS USE_OPTPP
						  )
endif()
