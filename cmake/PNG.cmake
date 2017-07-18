find_package(PNG REQUIRED)

message_var(PNG_LIBRARIES)

ADD_DEFINITIONS(-DUSE_PNG)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PNG
								  REQUIRED_VARS
								  PNG_LIBRARIES
								  )
