FILE(GLOB libpyfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

INSTALL(FILES ${libpyfiles}
  DESTINATION    lib
)