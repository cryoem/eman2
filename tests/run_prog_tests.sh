#!/usr/bin/env bash

# Run e2 programs by running the commands with -h

set -e

if [ -z ${CONDA_PREFIX+x} ];then
    source activate root
fi

# Gather programs from CONDA_PREFIX
if [ -d "${CONDA_PREFIX}"/bin ];then
    PROGS_DIR="${CONDA_PREFIX}"/bin
else
    PROGS_DIR="${CONDA_PREFIX}"/Library/bin
fi

progs=$(find "${PROGS_DIR}" -name 'e2*.py' | xargs -n 1 basename)

if [ -z ${progs+x} ];then
    echo "Cannot find any e2 programs in ${PROGS_DIR}"
    echo "CONDA_PREFIX: ${CONDA_PREFIX-:<not set>}"
fi

# Remove programs listed in "programs_no_test.txt"
MYDIR="$(cd "$(dirname "$0")"; pwd -P)"
progs_exclude=$(cat "${MYDIR}"/programs_no_test.txt | awk '{print $1}')

echo; echo "Removing programs from test list..."
for f in ${progs_exclude[@]};do
    echo "... $f"
    progs=( ${progs[@]/$f} )
done
echo

set +e

failed_progs=()
for prog in ${progs[@]};do
    echo "Running: $prog -h"
    $prog -h > /dev/null
    if [ $? -ne 0 ];then
        failed_progs+=($prog)
    fi
done

echo
echo "Total failed programs: ${#failed_progs[@]} / ${#progs[@]}"
for prog in ${failed_progs[@]};do
    echo ${prog}
done

if [ ${#failed_progs[@]} -ne 0 ];then
    exit 1
fi
