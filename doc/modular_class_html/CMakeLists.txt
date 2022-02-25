ADD_CUSTOM_TARGET(mkdoc ALL
	DEPENDS	${CMAKE_SOURCE_DIR}/doc/modular_class_html/docgen.py
)

ADD_CUSTOM_COMMAND(
	COMMAND cp
	ARGS -f ${CMAKE_SOURCE_DIR}/doc/modular_class_html/*.py .
	TARGET mkdoc
)

ADD_CUSTOM_COMMAND(
	COMMAND ./docgen.py
	TARGET mkdoc
)

FILE(GLOB modular_html "${CMAKE_CURRENT_SOURCE_DIR}/*.html")
INSTALL(FILES ${modular_html}
  DESTINATION    doc
)
