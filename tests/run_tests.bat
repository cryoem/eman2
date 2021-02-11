conda.exe list

:: 1. Run e2version.py and e2speedtest.py
e2version.py
e2speedtest.py

:: 2. Test imports
python tests\test_imports.py"

:: 3. Existence tests for data files like images, font files, JSON
python tests\test_EMAN2DIR.py
