#!/usr/bin/env bash

source ci_support/setup_conda.sh

export CPU_COUNT=2

conda install conda-build=2 -c defaults --yes --quiet

if [ "$(uname -s)" != "Darwin" ];then
    conda build recipes/eman -c cryoem -c defaults -c conda-forge --numpy 1.8
else
    export EMAN_TEST_SKIP=1
    conda build recipes/eman -c cryoem -c defaults -c conda-forge
fi
