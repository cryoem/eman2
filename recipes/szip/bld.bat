SET builddir=%SRC_DIR%\..\..\build_szip
if errorlevel 1 exit 1

mkdir %builddir% && cd %builddir%
if errorlevel 1 exit 1

set CL=/MP

cmake "%SRC_DIR%" -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release ^
                                       -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX%
if errorlevel 1 exit 1

cmake --build "%builddir%" --config Release --target install
if errorlevel 1 exit 1
