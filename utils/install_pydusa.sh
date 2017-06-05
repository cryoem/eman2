#!/usr/bin/env bash

# Builds and installs Pydusa against NumPy v1.8

set -xe

MYDIR=$(cd $(dirname $0) && pwd -P)

bash "${MYDIR}/build_pydusa_numpy.sh"   1.8
bash "${MYDIR}/install_pydusa_numpy.sh" 1.8
