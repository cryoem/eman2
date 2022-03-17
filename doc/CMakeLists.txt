IF(ENABLE_AUTODOC)
    ADD_SUBDIRECTORY(modular_class_html)
    ADD_SUBDIRECTORY(doxygen)

    add_custom_target(mkdoc
            DEPENDS mkdoc-doxy
            )

    FILE(GLOB static_html "${CMAKE_CURRENT_SOURCE_DIR}/*.html")
    INSTALL(FILES ${static_html}
      DESTINATION    doc/doxygen_html/
    )
ENDIF()

INSTALL(FILES EMAN2.0.big.jpg
    DESTINATION examples)
