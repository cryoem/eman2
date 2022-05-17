#!/usr/bin/env bash

# -*- mode: jinja-shell -*-

set -xeuo pipefail
export REPO_ROOT="/home/conda/repo_root"


export PYTHONUNBUFFERED=1

cat >~/.condarc <<CONDARC

conda-build:
 root-dir: ${REPO_ROOT}/../build_artifacts

CONDARC


mamba install --update-specs --yes --quiet --channel conda-forge \
    conda-build pip boa conda-forge-ci-setup=3
mamba update --update-specs --yes --quiet --channel conda-forge \
    conda-build pip boa conda-forge-ci-setup=3

/usr/bin/sudo -n yum install -y mesa-libGL


set +e
source ${REPO_ROOT}/ci_support/set_env_vars.sh
conda mambabuild "${REPO_ROOT}/recipe" -c cryoem -c conda-forge
