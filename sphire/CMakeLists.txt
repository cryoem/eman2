find_program(PYTHON "python")
# ADD_SUBDIRECTORY(bin)
# ADD_SUBDIRECTORY(libpy)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/build/timestamp
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND ${PYTHON} setup.py install
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/sphire/__init__.py
)

add_custom_target(target ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/build/timestamp)
