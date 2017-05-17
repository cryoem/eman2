#!/usr/bin/env bash

# Builds numpy version(s), if not specified build all versions, 1.5 - 1.12

set -xe

source activate root

RECIPES_DIR=$(cd $(dirname $0)/../recipes && pwd -P)

if [ $# -lt 1 ];then
    numpy_versions=( 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 )
else
    numpy_versions=${@}
fi

for v in ${numpy_versions[@]};do
    conda build "${RECIPES_DIR}/pydusa" --numpy ${v}
done
