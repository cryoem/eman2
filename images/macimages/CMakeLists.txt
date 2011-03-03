FILE(GLOB pngfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.png")
INSTALL(FILES ${pngfiles} 
  DESTINATION    images/macimages
)
