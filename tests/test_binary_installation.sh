#!/usr/bin/env bash

set -xe

dir="$1"
fname=$2
installer_file="${dir}/${fname}"
binary_loc="${dir}/eman2-binary-test"

rm -rf ${binary_loc}
bash "${installer_file}" -bp ${binary_loc}
