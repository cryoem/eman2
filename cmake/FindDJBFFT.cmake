SET(DJBFFT_LIBRARIES ${DJBFFT_LIBRARY})
SET(DJBFFT_INCLUDE_PATH ${DJBFFT_INCLUDE_PATH})
CHECK_OPTIONAL_LIB(DJBFFT fftc4.h fftr4.h)

message_var(DJBFFT_INCLUDE_PATH)
message_var(DJBFFT_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DJBFFT
								  REQUIRED_VARS DJBFFT_LIBRARY
								  )
