CHECK_REQUIRED_LIB(FTGL ftgl FTGL/FTGL.h ftgl FTGL/ftgl.h)
ADD_DEFINITIONS(-DUSE_FTGL)
IF(EXISTS ${FTGL_INCLUDE_PATH}/FTGL/FTGL.h AND COMMAND IF)
	ADD_DEFINITIONS(-DOLD_FTGL)
ENDIF()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FTGL
								  REQUIRED_VARS
									  FTGL_INCLUDE_PATH
									  FTGL_LIBRARY
								  )
