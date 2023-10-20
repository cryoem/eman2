#!/usr/bin/env bash

# -*- mode: jinja-shell -*-

set -xe

MINIFORGE_HOME=${HOME}/miniforge3

MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download"
MINIFORGE_FILE="Mambaforge-MacOSX-$(uname -m).sh"
curl -L -O "${MINIFORGE_URL}/${MINIFORGE_FILE}"
rm -rf ${MINIFORGE_HOME}
bash $MINIFORGE_FILE -b -p ${MINIFORGE_HOME}

source ${MINIFORGE_HOME}/etc/profile.d/conda.sh
conda activate base

mamba install --update-specs --quiet --yes --channel conda-forge \
    conda-build pip boa conda-forge-ci-setup=3 constructor
mamba update --update-specs --yes --quiet --channel conda-forge \
    conda-build pip boa conda-forge-ci-setup=3 constructor


/usr/bin/sudo mangle_homebrew
/usr/bin/sudo -k

source ./ci_support/set_env_vars.sh
conda mambabuild ./recipe -c cryoem -c conda-forge
constructor ./ci_support/ -v --debug
