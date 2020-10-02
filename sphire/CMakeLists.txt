# Install sphire via pip
INSTALL(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR} --no-compile)")
