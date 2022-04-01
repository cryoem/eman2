ADD_CUSTOM_TARGET(mkdoc-programs-help
	DEPENDS	${CMAKE_CURRENT_SOURCE_DIR}/generate_program_options_html.py
)

ADD_CUSTOM_COMMAND(
	TARGET mkdoc-programs-help
	COMMAND ./generate_program_options_html.py
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)


INSTALL(DIRECTORY programs_help_html
		DESTINATION    doc/doxygen_html
		)
