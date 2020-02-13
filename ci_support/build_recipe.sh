#!/usr/bin/env bash

set -xe

MYDIR="$(cd "$(dirname "$0")"; pwd -P)"

if [ -n "${TRAVIS}" ];then
    source ci_support/setup_conda.sh
fi

if [ -n "${CIRCLECI}" ];then
    . $HOME/miniconda/etc/profile.d/conda.sh
    conda activate eman
fi

python -m compileall -q .

if [ -n "$JENKINS_HOME" ];then
    export CPU_COUNT=4
else
    export CPU_COUNT=2
fi

conda info -a
conda list
conda list --explicit
conda render recipes/eman
conda build purge-all

if [ $AGENT_OS_NAME == "win" ];then
    CONDA_BUILD_TEST="--no-test"
else
    CONDA_BUILD_TEST=""
fi

conda build recipes/eman -c cryoem -c defaults -c conda-forge ${CONDA_BUILD_TEST}
