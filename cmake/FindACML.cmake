CHECK_OPTIONAL_LIB(ACML acml acml.h)
FIND_LIBRARY(G2C_LIBRARY NAMES g2c PATHS
			 /usr/lib64
			 /usr/lib
			 /usr/local/lib
			 $ENV{HOME}/lib
			 )

message_var(ACML_INCLUDE_PATH)
message_var(ACML_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ACML
								  REQUIRED_VARS
								  ACML_INCLUDE_PATH
								  ACML_LIBRARY
								  )
