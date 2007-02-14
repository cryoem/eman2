IF (ENABLE_AUTODOC)
	SUBDIRS(latex modular_class_html doxygen)
ENDIF (ENABLE_AUTODOC)

MAKE_DIRECTORY(${EMAN_INSTALL_PREFIX}/doc)
SET(DOC_PATH ${EMAN_INSTALL_PREFIX}/doc)

INSTALL_FILES(/doc "^.+\\.(html)$")
INSTALL_FILES(/doc FILES WindowsVisualStudioDotNet2003CompilationNotes.txt INSTALL EMAN2.0.big.jpg)
