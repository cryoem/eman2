#!/usr/bin/env bash

if [ $# -ne 2 ];then
    printf "\e\033[35m\n  Usage: $(basename ${0})   %s   %s\033[0m\n\n" "output-dir" "construct.yaml-dir" >&2
    exit 64
fi

set -xe

OUTPUT_DIR=$1
CONSTRUCT_YAML_DIR=$2

export PYTHONUNBUFFERED=1

# Package eman
mkdir -p ${OUTPUT_DIR} && cd ${OUTPUT_DIR}

CONSTRUCT_YAML="${CONSTRUCT_YAML_DIR}/construct.yaml"
CONDA_ROOT_DIR="$(conda info --root)"
CONDA_ROOT_DIR=${CONDA_ROOT_DIR//\\/\\\\}
sed -i.bak "s~place_holder_conda_prefix~"${CONDA_ROOT_DIR}"~" ${CONSTRUCT_YAML}
cat ${CONSTRUCT_YAML}
constructor --version
constructor ${CONSTRUCT_YAML_DIR} -v --cache-dir=${HOME_DIR}/.conda/constructor
mv ${CONSTRUCT_YAML}.bak ${CONSTRUCT_YAML}
