#!/usr/bin/env bash

# Builds numpy version(s), if not specified build all versions, 1.5 - 1.12

# unset args before using, it is set in activate script above
unset args
for elem in ${@};do
    regex=-*
    if [[ $elem == $regex ]];then
        opts=( ${opts[@]} $elem )
    else
        args=( ${args[@]} $elem )
    fi
done

set -xe

RECIPES_DIR=$(cd $(dirname $0)/../recipes && pwd -P)

if [ ${#args} -lt 1 ];then
    numpy_versions=( 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 )
else
    numpy_versions=${args[@]}
fi

for v in ${numpy_versions[@]};do
    conda build "${RECIPES_DIR}/pydusa" --numpy ${v} ${opts[@]}
done
