#!/usr/bin/env bash

if [ $# -ne 1 ];then
    printf "\e\033[35m\n  Usage: $(basename ${0})   %s   \033[0m\n\n" "eman-recipe-dir" >&2
    exit 64
fi

set -xe

EMAN_RECIPE_DIR=$1

export PYTHONUNBUFFERED=1

# Build eman recipe
conda info -a
conda render ${EMAN_RECIPE_DIR}
conda build purge-all
conda build ${EMAN_RECIPE_DIR} -c cryoem -c defaults -c conda-forge --quiet
