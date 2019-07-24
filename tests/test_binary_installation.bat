@echo on

set "BASH_EXE=C:\Program Files\Git\bin\bash.exe"

set "installer_file=%1"
set "installation_loc=%2"

rmdir /q /s %installation_loc%
start /wait "" %installer_file% /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S /D=%installation_loc%
if errorlevel 1 exit 1

call %installation_loc%\Scripts\activate.bat
if errorlevel 1 exit 1

call tests\run_tests.bat
if errorlevel 1 exit 1

"%BASH_EXE%" -c "bash tests/run_prog_tests.sh"
if errorlevel 1 exit 1
