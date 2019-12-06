SET builddir=%SRC_DIR%\..\build_eman
if errorlevel 1 exit 1

mkdir %builddir% && cd %builddir%
if errorlevel 1 exit 1

set CL=/MP

cmake --version
cmake "%SRC_DIR%" -G "NMake Makefiles" ^
                    -DCMAKE_BUILD_TYPE=Release    ^
                    -DENABLE_WARNINGS=OFF
if errorlevel 1 exit 1

nmake
if errorlevel 1 exit 1

nmake install
if errorlevel 1 exit 1
