CHECK_REQUIRED_LIB(NFFT2 nfft nfft.h "" "")
target_link_libraries(NFFT2 FFTW3D)

message_var(NFFT2_INCLUDE_PATH)
message_var(NFFT2_LIBRARY)

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
endif()
