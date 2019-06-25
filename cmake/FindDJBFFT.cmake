CHECK_REQUIRED_LIB(DJBFFT fftc4.h fftr4.h "" "")

cmake_print_variables(DJBFFT_INCLUDE_PATH)
cmake_print_variables(DJBFFT_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DJBFFT
								  REQUIRED_VARS DJBFFT_LIBRARY
								  )

if(DJBFFT_FOUND AND NOT TARGET DJBFFT)
	add_library(DJBFFT INTERFACE)
	add_library(DJBFFT::DJBFFT ALIAS DJBFFT)
	set_target_properties(DJBFFT PROPERTIES
						  INTERFACE_INCLUDE_DIRECTORIES "${DJBFFT_INCLUDE_PATH}"
						  INTERFACE_LINK_LIBRARIES      "${DJBFFT_LIBRARY}"
						  INTERFACE_COMPILE_DEFINITIONS USE_DJBFFT
						  )
endif()
