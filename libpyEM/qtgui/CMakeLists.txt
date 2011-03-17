FILE(GLOB qtguifiles "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
FILE(GLOB exclusionlist "em3Dhelloworld.py")
LIST(REMOVE_ITEM qtguifiles ${exclusionlist})
INSTALL(FILES ${qtguifiles} DESTINATION lib)