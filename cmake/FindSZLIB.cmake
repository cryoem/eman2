CHECK_LIB_ONLY(SZLIB libszip)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SZLIB
								  REQUIRED_VARS SZLIB_LIBRARY
								  )
