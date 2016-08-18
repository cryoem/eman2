FILE(GLOB testprograms "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
INSTALL(PROGRAMS ${testprograms}
  DESTINATION    test/rt
)
