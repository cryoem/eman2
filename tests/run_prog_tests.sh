#!/usr/bin/env bash

# Run e2 programs by running the commands with -h

set -e

MYDIR="$(cd "$(dirname "$0")"; pwd -P)"
PROGS_DIR="${MYDIR}/../programs"

progs=$(find "${PROGS_DIR}" -name 'e2*.py' -print0 | while read -d '' prog; do basename "$prog"; done)
progs_exclude=$(cat "${MYDIR}"/programs_no_test.txt | awk '{print $1}')

echo -e "\nRemoving programs from test list..."
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

echo -e "\nTotal failed programs: ${#failed_progs[@]} / ${#progs[@]}"
for prog in ${failed_progs[@]};do
    echo ${prog}
done

if [ ${#failed_progs[@]} -ne 0 ];then
    exit 1
fi
