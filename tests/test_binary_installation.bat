@echo on

set "dir=%1"
set fname=%2
set "installer_file=%dir%\%fname%"
set "binary_loc=%dir%\eman2-binary-test"

rmdir /q /s %binary_loc%
start /wait "" %installer_file% /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S /D=%binary_loc%
if errorlevel 1 exit 1
