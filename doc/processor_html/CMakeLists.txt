ADD_CUSTOM_TARGET( mkdoc ALL
	DEPENDS	${EMAN_SOURCE_DIR}/doc/processor_html/docgen.py
)

ADD_CUSTOM_COMMAND(
	SOURCE  mkdoc
	COMMAND cp 
	ARGS -f ${EMAN_SOURCE_DIR}/doc/processor_html/*.py .
	TARGET mkdoc
)

ADD_CUSTOM_COMMAND(
	SOURCE  mkdoc
	COMMAND ./docgen.py	
	TARGET mkdoc
)

INSTALL_FILES(/doc FILES processors.html processor_groups.html)
