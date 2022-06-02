CHECK_REQUIRED_LIB(NFFT2 nfft nfft.h "" "")

cmake_print_variables(NFFT2_INCLUDE_PATH)
cmake_print_variables(NFFT2_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NFFT2
								  REQUIRED_VARS NFFT2_INCLUDE_PATH NFFT2_LIBRARY
								  )

if(NFFT2_FOUND AND NOT TARGET NFFT2)
	add_library(NFFT2 INTERFACE)
	set_target_properties(NFFT2 PROPERTIES
						  INTERFACE_INCLUDE_DIRECTORIES "${NFFT2_INCLUDE_PATH}"
						  INTERFACE_LINK_LIBRARIES      "${NFFT2_LIBRARY}"
						  INTERFACE_COMPILE_DEFINITIONS USE_NFFT2
						  )
		
	target_link_libraries(NFFT2 INTERFACE FFTW3D::FFTW3D)
endif()
