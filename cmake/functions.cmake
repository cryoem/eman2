include(CMakePrintHelpers)

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

OPTION(DEBUG_CHECK_REQUIRED_LIB "enable debug output for function CHECK_REQUIRED_LIB" OFF)
function(CHECK_REQUIRED_LIB upper lower header lower2 header2)
	if(DEBUG_CHECK_REQUIRED_LIB)
		message("\n### BEGIN ### CHECK_REQUIRED_LIB ${upper} ${lower} ${header} ${lower2} ${header2}")
		cmake_print_variables(${upper}_INCLUDE_PATH)
		cmake_print_variables(${upper}_LIBRARY)
		
		message("Searching in ${EMAN_PREFIX_INC} for ${header} and ${header2} ...")
	endif()
		
	FIND_PATH(${upper}_INCLUDE_PATH
			NAMES ${header} ${header2}
			PATHS $ENV{${upper}DIR}/include ${EMAN_PREFIX_INC}
			NO_DEFAULT_PATH
			)
	
	IF(${upper}_INCLUDE_PATH)
		FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} ${lower2} PATHS $ENV{${upper}DIR}/lib ${EMAN_PREFIX_LIB})
		
		if(DEBUG_CHECK_REQUIRED_LIB)
			message("FIND_LIBRARY: ${upper}_LIBRARY NAMES ${lower} ${lower2} PATHS $ENV{${upper}DIR}/lib ${EMAN_PREFIX_LIB}")
		endif()
	ENDIF()
	
	IF(NOT ${upper}_INCLUDE_PATH OR NOT ${upper}_LIBRARY)
		MESSAGE(SEND_ERROR "ERROR: ${upper} not found. please install ${upper} first!")
	ENDIF()

	if(DEBUG_CHECK_REQUIRED_LIB)
		message("### END ### CHECK_REQUIRED_LIB ${upper} ${lower} ${header} ${lower2} ${header2}")
		cmake_print_variables(${upper}_INCLUDE_PATH)
		cmake_print_variables(${upper}_LIBRARY)
	endif()
endfunction()

function(CHECK_LIB_ONLY upper lower)
	FIND_LIBRARY(${upper}_LIBRARY NAMES ${lower} PATHS $ENV{${upper}DIR}/lib ${EMAN_PREFIX_LIB} NO_DEFAULT_PATH)
	message(STATUS "CHECK_LIB_ONLY upper lower")
	cmake_print_variables(${upper}_LIBRARY)
endfunction()

function(check_if_options_on_greater_than opts num)
	foreach(opt ${opts})
		list(APPEND options "${opt}:${${opt}}")
	endforeach()

	list(FILTER options EXCLUDE REGEX "OFF")
	list(LENGTH options count)
	list(JOIN opts "\n" opts)
		
	if(count GREATER ${num})
		message("\nOnly ${num} or less of the options below can be ON:\n${opts}\n")
		message("Found ${count}:")
		foreach(opt ${options})
			message(${opt})
		endforeach()
		message("\n")
		message(FATAL_ERROR "")
	endif()
endfunction()

function(check_if_options_not_equal_to opts num)
	foreach(opt ${opts})
		list(APPEND options "${opt}:${${opt}}")
	endforeach()

	list(FILTER options EXCLUDE REGEX "OFF")
	list(LENGTH options count)
	list(JOIN opts "\n" opts)
	
	if(NOT count EQUAL ${num})
		message("\nOnly ${num} of the options below can be ON:\n${opts}\n")
		message("Found ${count}:")
		foreach(opt ${options})
			message(${opt})
		endforeach()
		message("\n")
		message(FATAL_ERROR "")
	endif()
endfunction()
