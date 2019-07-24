#!/usr/bin/env bash

set -xe

dir="$1"
fname=$2
installer_file="${dir}/${fname}"
installation_loc="${dir}/eman2-binary-test"
MYDIR="$(cd "$(dirname "$0")" && pwd -P)"

rm -rf ${installation_loc}
bash "${installer_file}" -bp ${installation_loc}
source ${installation_loc}/bin/activate root
conda info -a
conda list
conda list --explicit
bash "${MYDIR}/run_tests_from_binary.sh"
