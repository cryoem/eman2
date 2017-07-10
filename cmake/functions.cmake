FUNCTION(message_var var)
	message("${var}: ${${var}}")
endfunction()

function(set_cache_var var val)
	if(val AND NOT ${var})
		set(${var} ${val} CACHE PATH "" FORCE)
	else()
		set(${var} ${var}-NOTFOUND CACHE PATH "")
	endif()
endfunction()

function(set_cache_var_to_var var val)
	if(${val} AND NOT ${var})
		set(${var} ${${val}} CACHE PATH "" FORCE)
	else()
		set(${var} ${var}-NOTFOUND CACHE PATH "")
	endif()
endfunction()

function(EMAN_CHECK_FUNCTION FUNCTION VARIABLE)
	CHECK_FUNCTION_EXISTS(${FUNCTION} ${VARIABLE})
	IF(${VARIABLE})
		ADD_DEFINITIONS(-D${VARIABLE})
	ENDIF()
endfunction()

function(CHECK_REQUIRED_LIB2 upper lower header lower2 header2)
	message("\n### BEGIN ### CHECK_REQUIRED_LIB2 upper lower header lower2 header2")
	message_var(${upper}_INCLUDE_PATH)
	message_var(${upper}_LIBRARY)
	
	IF(NOT ${upper}_INCLUDE_PATH)
		FIND_PATH(${upper}_INCLUDE_PATH
				NAMES ${header} ${header2}
				PATHS $ENV{${upper}DIR}/include ${EMAN_PREFIX_INC}
				NO_DEFAULT_PATH
				)
	ENDIF()
	
	IF(${upper}_INCLUDE_PATH)
		FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} ${lower2} PATHS $ENV{${upper}DIR}/lib ${EMAN_PREFIX_LIB})
		message("FIND_LIBRARY: ${upper}_LIBRARY NAMES ${lower} ${lower2} PATHS $ENV{${upper}DIR}/lib ${EMAN_PREFIX_LIB}")
	ENDIF()
	
	IF(NOT ${upper}_INCLUDE_PATH OR NOT ${upper}_LIBRARY)
		MESSAGE(SEND_ERROR "ERROR: ${upper} not found. please install ${upper} first!")
	ENDIF()
	
	message("### END ### CHECK_REQUIRED_LIB2 upper lower header lower2 header2")
	message_var(${upper}_INCLUDE_PATH)
	message_var(${upper}_LIBRARY)
endfunction()

function(CHECK_OPTIONAL_LIB upper lower header)
	FIND_PATH(${upper}_INCLUDE_PATH ${header} $ENV{${upper}DIR}/include ${EMAN_PREFIX_INC} NO_DEFAULT_PATH)
	
	IF(${upper}_INCLUDE_PATH)
		FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} PATHS $ENV{${upper}DIR}/lib ${EMAN_PREFIX_LIB})
		
		IF(${upper}_LIBRARY)
			ADD_DEFINITIONS(-DEM_${upper})
			INCLUDE_DIRECTORIES(${${upper}_INCLUDE_PATH})
		ENDIF()
	ENDIF()
	IF(NOT ${upper}_INCLUDE_PATH AND NOT ${upper}_LIBRARY)
		OPTION(ENABLE_${upper} "if enable ${upper} support" OFF)
	endif()
endfunction()

function(CHECK_LIB_ONLY upper lower)
	FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} PATHS $ENV{${upper}DIR}/lib ${EMAN_PREFIX_LIB} NO_DEFAULT_PATH)
	message(STATUS "CHECK_LIB_ONLY upper lower")
	message_var(${upper}_LIBRARY)
endfunction()
