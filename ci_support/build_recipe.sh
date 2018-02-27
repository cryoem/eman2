#!/usr/bin/env bash

if [ ! -z ${CI} ];then
    source ci_support/setup_conda.sh

    conda install conda-build=2 -c defaults --yes --quiet
fi

export CPU_COUNT=2
export NUMPY_VERSION=1.9

if [ "$(uname -s)" == "Darwin" ];then
    export EMAN_TEST_SKIP=1
fi

conda build recipes/eman -c cryoem -c defaults -c conda-forge
