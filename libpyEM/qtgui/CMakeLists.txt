add_custom_target(QtGuiFiles
				  COMMAND ${CMAKE_COMMAND} -DCOMPONENT=QtGuiFiles -P cmake_install.cmake
				  )

FILE(GLOB qtguifiles "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
FILE(GLOB exclusionlist "em3Dhelloworld.py")

# MESSAGE("gui files: ${qtguifiles}")

LIST(REMOVE_ITEM qtguifiles ${exclusionlist})

# MESSAGE("reduced gui files: ${qtguifiles}")

INSTALL(FILES ${qtguifiles}
		DESTINATION  ${SP_DIR}
		COMPONENT QtGuiFiles
		)
