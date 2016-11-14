FILE(GLOB tbxprograms "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
INSTALL(FILES ${tbxprograms}
  DESTINATION    lib/pyemtbx
)