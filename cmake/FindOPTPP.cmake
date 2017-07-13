FIND_PATH(OPTPP_INCLUDE_PATH Opt.h /usr/include /usr/include/opt++
		  /usr/local/include /usr/local/include/opt++ $ENV{HOME}/include
		  $ENV{HOME}/include/opt++ NO_DEFAULT_PATH)
FIND_PATH(NEWMAT_INCLUDE_PATH newmat.h /usr/include /usr/include/opt++
		  /usr/include/newmat /usr/local/include/newmat
		  /usr/local/include /usr/local/include/opt++ $ENV{HOME}/include
		  $ENV{HOME}/include/opt++ NO_DEFAULT_PATH)
FIND_LIBRARY(OPTPP_LIBRARY NAMES opt-linux PATHS /usr/lib /usr/local/lib)
FIND_LIBRARY(NEWMAT_LIBRARY NAMES newmat-linux PATHS /usr/lib /usr/local/lib)

message_var(OPTPP_INCLUDE_PATH)
message_var(OPTPP_LIBRARY)
message_var(NEWMAT_INCLUDE_PATH)
message_var(NEWMAT_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OPTPP
								  REQUIRED_VARS
								  	OPTPP_INCLUDE_PATH
								  	OPTPP_LIBRARY
								  	NEWMAT_INCLUDE_PATH
								  	NEWMAT_LIBRARY
								  )

ADD_DEFINITIONS(-DUSE_OPTPP)
INCLUDE_DIRECTORIES(${OPTPP_INCLUDE_PATH} ${NEWMAT_INCLUDE_PATH})
