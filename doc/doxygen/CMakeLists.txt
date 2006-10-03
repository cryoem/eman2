ADD_CUSTOM_TARGET( mkdoc ALL
	DEPENDS	${EMAN_SOURCE_DIR}/doc/doxygen/Doxyfile
)

ADD_CUSTOM_COMMAND(
	SOURCE  mkdoc
	COMMAND cp 
	ARGS -f ${EMAN_SOURCE_DIR}/doc/doxygen/makedoc.sh .
	TARGET mkdoc
)

ADD_CUSTOM_COMMAND(
	SOURCE  mkdoc
	COMMAND ./makedoc.sh	
	TARGET mkdoc
)

