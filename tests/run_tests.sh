#!/usr/bin/env bash

set -xe

conda list

MYDIR="$(cd "$(dirname "$0")"; pwd -P)"

# 1. Run e2version.py and e2speedtest.py
e2version.py
e2speedtest.py
time python -c "from EMAN2 import *; Util.version()"

# 2. Test git-hash consistency
if [ ${CONDA_BUILD:-0} -ne 1 ];then
    bash "${MYDIR}/test_git_hash.sh"
fi

# 3. Existence tests for data files like images, font files, JSON
python "${MYDIR}/test_EMAN2DIR.py"

# 4. Unit tests
nosetests -vv --exe -m "^test_*" \
                    -e "^test_image_" \
                    -e "test_main" \
                    -e "test_result" \
                    -e "test_boxing" \
                    -a \!broken \
                    rt/pyem/

# 5. Test openmpi
if [ $(whoami) != "root" ] && [ "$(uname -s)" != "Darwin" ];then
    mpirun -n 4 $(which python) ${MYDIR}/../examples/mpi_test.py
fi

# 6. Run e2*.py -h
bash "${MYDIR}/run_prog_tests.sh"
