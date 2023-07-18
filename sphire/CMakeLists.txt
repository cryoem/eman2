# Install sphire via pip
add_custom_target(install-sphire
        COMMAND ${Python3_EXECUTABLE} -m pip install ${CMAKE_CURRENT_SOURCE_DIR} --no-compile
        )
