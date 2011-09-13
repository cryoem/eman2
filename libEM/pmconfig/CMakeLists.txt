FILE(GLOB pmfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.*")
INSTALL(FILES ${pmfiles}
  DESTINATION    lib/pmconfig
)
