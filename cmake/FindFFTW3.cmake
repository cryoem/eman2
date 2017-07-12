CHECK_REQUIRED_LIB(FFTW3F fftw3f fftw3.h libfftw3f-3 "")
CHECK_REQUIRED_LIB(FFTW3D fftw3 fftw3.h libfftw3-3  "")

SET(FFTW_LIBRARIES ${FFTW3F_LIBRARY} ${FFTW3D_LIBRARY})
SET(FFTW_INCLUDE_PATH ${FFTW3_INCLUDE_PATH})

if(ENABLE_CONDA)
	message("FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
	message("FFTW_INCLUDE_PATH: ${FFTW_INCLUDE_PATH}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
		REQUIRED_VARS
			FFTW3F_INCLUDE_PATH
			FFTW3F_LIBRARY
			FFTW3D_INCLUDE_PATH
			FFTW3D_LIBRARY
		)

if(FFTW3_FOUND AND NOT TARGET FFTW3 AND NOT TARGET FFTW3F AND NOT TARGET FFTW3D)
	foreach(trgt FFTW3F FFTW3D FFTW3)
		add_library(${trgt} INTERFACE)
		add_library(${trgt}::${trgt} ALIAS ${trgt})
	endforeach()
	foreach(trgt FFTW3F FFTW3D)
		set_target_properties(${trgt} PROPERTIES
							  INTERFACE_INCLUDE_DIRECTORIES "${${trgt}_INCLUDE_PATH}"
							  INTERFACE_LINK_LIBRARIES      "${${trgt}_LIBRARY}"
							  )
	endforeach()
		
	set_target_properties(FFTW3 PROPERTIES
						  INTERFACE_COMPILE_DEFINITIONS USE_FFTW3
						  )
	target_link_libraries(FFTW3 INTERFACE FFTW3F::FFTW3F FFTW3D::FFTW3D)
endif()
