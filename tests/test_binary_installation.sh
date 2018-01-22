#!/usr/bin/env bash

set -xe

MYDIR="$(cd "$(dirname "$0")" && pwd -P)"

for f in ${@};do
    dir=$(cd $(dirname $f) && pwd -P)
    fbase=$(basename $f)
    fbase=${fbase%\.*}
    fbase=${fbase//\./-}
    conda_loc=${dir}/${fbase}
    
    echo "... $fbase ..."
    
    bash $f -bf -p ${conda_loc}
    source ${conda_loc}/bin/activate root
    conda info -a
    conda list
    bash "${MYDIR}/run_tests_from_binary.sh"
    source deactivate
done
