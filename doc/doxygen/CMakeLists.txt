ADD_CUSTOM_TARGET(mkdoc-doxy
	DEPENDS	${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
	DEPENDS mkdoc-modular-class
)

ADD_CUSTOM_COMMAND(
	TARGET mkdoc-doxy
	COMMAND ./makedoc.sh
)
