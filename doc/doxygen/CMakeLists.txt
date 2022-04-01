ADD_CUSTOM_TARGET(mkdoc-doxy
	DEPENDS	${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile ${CMAKE_CURRENT_SOURCE_DIR}/makedoc.sh
	DEPENDS mkdoc-modular-class mkdoc-programs-help
)

ADD_CUSTOM_COMMAND(
	TARGET mkdoc-doxy
	COMMAND ./makedoc.sh
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

INSTALL(DIRECTORY doxygen_html
		DESTINATION    doc
		)
