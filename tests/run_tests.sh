#!/usr/bin/env bash

set -xe

conda list

MYDIR="$(cd "$(dirname "$0")"; pwd -P)"

python -m compileall -q -x .git -x sparx -x sphire .

# 1. Run e2version.py and e2speedtest.py
e2version.py
e2speedtest.py

# 2. Test git-hash consistency
if [ ${CONDA_BUILD:-0} -ne 1 ];then
    bash "${MYDIR}/test_git_hash.sh"
fi

# 3. Test imports
python "${MYDIR}/test_imports.py"

# 4. Existence tests for data files like images, font files, JSON
python "${MYDIR}/test_EMAN2DIR.py"

# 5. Unit tests
nosetests -vv --exe -m "^test_*" \
                    -e "^test_image_" \
                    -e "test_main" \
                    -e "test_result" \
                    -e "test_boxing" \
                    -a \!broken \
                    rt/pyem/

# 6. Test openmpi
if [ $(whoami) != "root" ];then
    mpirun --oversubscribe -n 4 $(which python) ${MYDIR}/../examples/mpi_test.py
fi

# 7. Run e2*.py -h
python "${MYDIR}/run_prog_tests.py"
