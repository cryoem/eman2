INCLUDE_DIRECTORIES(${NFFT2_INCLUDE_PATH})
CHECK_OPTIONAL_LIB(NFFT2 nfft nfft.h)
target_link_libraries(NFFT2 FFTW3D)

message_var(NFFT2_INCLUDE_PATH)
message_var(NFFT2_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NFFT2
								  REQUIRED_VARS NFFT2_INCLUDE_PATH NFFT2_LIBRARY
								  )
