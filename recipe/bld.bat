@echo on

SET builddir=%SRC_DIR%\..\build_eman
if errorlevel 1 exit 1

mkdir %builddir% && cd %builddir%
if errorlevel 1 exit 1

set CL=/MP

cmake --version
cmake "%SRC_DIR%" -G "%CMAKE_GENERATOR%"       ^
                    -DPython3_EXECUTABLE="%PYTHON%" ^
                    -DCMAKE_BUILD_TYPE=Release ^
                    -DENABLE_WARNINGS=OFF
if errorlevel 1 exit 1

cmake --build "%builddir%" --config Release --target install
if errorlevel 1 exit 1
