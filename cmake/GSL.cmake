find_package(GSL REQUIRED)

if(WIN32)
	set_target_properties(GSL::gsl
						  PROPERTIES
						  INTERFACE_COMPILE_DEFINITIONS GSL_DLL
						  )
endif()
