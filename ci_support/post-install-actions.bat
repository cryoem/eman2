call %PREFIX%\Scripts\activate.bat
if errorlevel 1 exit 1

conda.exe config --env --set auto_update_conda False
if errorlevel 1 exit 1

mkdir %PREFIX%\install_logs
if errorlevel 1 exit 1

conda.exe info -a          > %PREFIX%\install_logs\info_log.txt 2>&1
if errorlevel 1 exit 1
conda.exe list             > %PREFIX%\install_logs\list_log.txt 2>&1
if errorlevel 1 exit 1
conda.exe list --explicit >> %PREFIX%\install_logs\list_log.txt 2>&1
if errorlevel 1 exit 1


if exist site-packages\nul (
    xcopy /e /y site-packages\* Lib\site-packages\
    rmdir /s /q site-packages
) > %PREFIX%\install_logs\install_log.txt 2>&1

if not defined EMAN_INSTALL_DONT_UPDATE_DEPS (
    if not "%EMAN_INSTALL_DONT_UPDATE_DEPS%"=="0" (
        conda.exe install -v eman-deps=29.1 -c cryoem -c defaults -c conda-forge -y >> %PREFIX%\install_logs\install_log.txt 2>&1
    ) else (
        echo "WARNING: Skipping installation of dependencies per user request..."
    )
) else (
    echo "WARNING: Skipping installation of dependencies per user request..."
)
if errorlevel 1 exit 1
