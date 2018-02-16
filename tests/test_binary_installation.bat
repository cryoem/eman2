@echo on

set "BASH_EXE=C:\Program Files\Git\bin\bash.exe"

set "dir=%1"
set fname=%2
set "installer_file=%dir%\%fname%"
set "binary_loc=%dir%\eman2-binary-test"

rmdir /q /s %binary_loc%
start /wait "" %installer_file% /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S /D=%binary_loc%
if errorlevel 1 exit 1

call %binary_loc%\Scripts\activate.bat
if errorlevel 1 exit 1

call tests\run_tests.bat
if errorlevel 1 exit 1

"%BASH_EXE%" -c "bash tests/run_prog_tests.sh"
if errorlevel 1 exit 1
