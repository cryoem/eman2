FILE(GLOB qtguifiles "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

# MESSAGE("gui files: ${qtguifiles}")

INSTALL(FILES ${qtguifiles}
		DESTINATION  ${Python3_SITELIB}/eman2_gui
		COMPONENT PythonFiles
		)
