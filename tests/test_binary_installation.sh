#!/usr/bin/env bash

set -xe

dir="$1"
fname=$2
installer_file="${dir}/${fname}"
binary_loc="${dir}/eman2-binary-test"
MYDIR="$(cd "$(dirname "$0")" && pwd -P)"

rm -rf ${binary_loc}
bash "${installer_file}" -bp ${binary_loc}
source ${binary_loc}/bin/activate root
conda info -a
conda list
conda list --explicit
bash "${MYDIR}/run_tests_from_binary.sh"
