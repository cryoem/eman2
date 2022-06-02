ADD_CUSTOM_TARGET(mkdoc-modular-class
	DEPENDS	${CMAKE_CURRENT_SOURCE_DIR}/docgen.py
)

ADD_CUSTOM_COMMAND(
	TARGET mkdoc-modular-class
	COMMAND ./docgen.py
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

INSTALL(DIRECTORY modular_class_html
		DESTINATION    doc/doxygen_html
		)
