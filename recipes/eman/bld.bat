SET builddir=%SRC_DIR%\..\build_eman
if errorlevel 1 exit 1

mkdir %builddir% && cd %builddir%
if errorlevel 1 exit 1

set CL=/MP

cmake --version
cmake "%SRC_DIR%" -G "Visual Studio 9 2008 Win64" ^
                    -DCMAKE_BUILD_TYPE=Release    ^
                    -DENABLE_WARNINGS=OFF
if errorlevel 1 exit 1

cmake --build "%builddir%" --config Release --target install
if errorlevel 1 exit 1
