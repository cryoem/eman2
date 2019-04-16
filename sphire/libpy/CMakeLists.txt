ADD_SUBDIRECTORY(prior_calculation)

FILE(GLOB libpyfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

INSTALL(FILES ${libpyfiles}
  DESTINATION    ${SP_DIR}
)
