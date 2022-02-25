ADD_CUSTOM_TARGET(mkdoc ALL
	DEPENDS	${CMAKE_SOURCE_DIR}/doc/doxygen/Doxyfile
)

ADD_CUSTOM_COMMAND(
	COMMAND cp
	ARGS -f ${CMAKE_SOURCE_DIR}/doc/doxygen/makedoc.sh .
	TARGET mkdoc
)

ADD_CUSTOM_COMMAND(
	COMMAND ./makedoc.sh
	TARGET mkdoc
)
