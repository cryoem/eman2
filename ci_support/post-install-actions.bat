call %PREFIX%\Scripts\activate.bat
if errorlevel 1 exit 1

mkdir %PREFIX%\install_logs
if errorlevel 1 exit 1

conda.exe info -a          > %PREFIX%\install_logs\info_log.txt 2>&1
if errorlevel 1 exit 1
conda.exe list             > %PREFIX%\install_logs\list_log.txt 2>&1
if errorlevel 1 exit 1
conda.exe list --explicit >> %PREFIX%\install_logs\list_log.txt 2>&1
if errorlevel 1 exit 1


conda.exe install --force-reinstall conda=4.6.14 conda-build=3.17.8 pytz backports backports.functools_lru_cache filelock tqdm -y > %PREFIX%\install_logs\install_log.txt 2>&1
if errorlevel 1 exit 1
conda.exe install -v eman-deps=14.1 -c cryoem -c defaults -c conda-forge -y > %PREFIX%\install_logs\install_log.txt 2>&1
if errorlevel 1 exit 1
