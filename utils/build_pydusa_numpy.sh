#!/usr/bin/env bash

# Builds numpy version(s), if not specified build all versions, 1.5 - 1.12

set -xe

source activate root

RECIPES_DIR=$(cd $(dirname $0)/../recipes && pwd -P)

numpy_versions=( 5 6 7 8 9 10 11 12 )
for v in ${numpy_versions[@]};do
    conda build "${RECIPES_DIR}/pydusa" --numpy 1.${v}
done
