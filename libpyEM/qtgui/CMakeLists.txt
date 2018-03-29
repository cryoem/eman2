FILE(GLOB qtguifiles "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

# MESSAGE("gui files: ${qtguifiles}")

INSTALL(FILES ${qtguifiles}
		DESTINATION  ${SP_DIR}/eman2_gui
		COMPONENT PythonFiles
		)
