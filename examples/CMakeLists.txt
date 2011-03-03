FILE(GLOB example_programs "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
INSTALL(PROGRAMS ${example_programs}
  DESTINATION    examples
)

INSTALL(FILES 00README DESTINATION examples)