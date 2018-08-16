call %PREFIX%\Scripts\activate.bat

conda.exe info -a
conda.exe list

conda.exe install eman-deps=12.2 -c cryoem -c defaults -c conda-forge -y
