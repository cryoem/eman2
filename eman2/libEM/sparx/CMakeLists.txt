FILE(GLOB headerfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${headerfiles}
  DESTINATION    include
)
