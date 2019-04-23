call %PREFIX%\Scripts\activate.bat

mkdir %PREFIX%\install_logs

conda.exe info -a > %PREFIX%\install_logs\info_log.txt 2>&1
conda.exe list    > %PREFIX%\install_logs\list_log.txt 2>&1

conda.exe install -v eman-deps=14.1 -c cryoem -c defaults -c conda-forge -y > %PREFIX%\install_logs\install_log.txt 2>&1
