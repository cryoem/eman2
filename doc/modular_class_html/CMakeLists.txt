ADD_CUSTOM_TARGET( mkdoc ALL
	DEPENDS	${EMAN_SOURCE_DIR}/doc/modular_class_html/docgen.py
)

ADD_CUSTOM_COMMAND(
	SOURCE  mkdoc
	COMMAND cp 
	ARGS -f ${EMAN_SOURCE_DIR}/doc/modular_class_html/*.py .
	TARGET mkdoc
)

ADD_CUSTOM_COMMAND(
	SOURCE  mkdoc
	COMMAND ./docgen.py	
	TARGET mkdoc
)

FILE(GLOB modular_html "${CMAKE_CURRENT_SOURCE_DIR}/*.html")
INSTALL(FILES ${modular_html}
  DESTINATION    doc
)

#INSTALL_FILES(/doc FILES processors.html processor_groups.html aligners.html analyzers.html averagers.html cmps.html projectors.html reconstructors.html)
