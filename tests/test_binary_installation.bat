@echo on

set "installer_file=%1"
set "installation_loc=%2"

rmdir /q /s "%installation_loc%"

start /wait "" "%installer_file%" /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S /D=%installation_loc%
if errorlevel 1 exit 1

call "%installation_loc%\Scripts\activate.bat"
if errorlevel 1 exit 1

conda.exe info -a
conda.exe list
conda.exe list --explicit

call tests\run_tests.bat
if errorlevel 1 exit 1
