#!/usr/bin/env bash

# Builds numpy version(s), if not specified build all versions, 1.5 - 1.12

function print_usage(){
    printf "\e\033[35m\n  Usage: $(basename ${0}) %s %s %s\n\n\033[0m" "{1.5|...|1.12}" "[branch-to-checkout]" "[--conda-build-arg1] [--conda-build-arg2] ..." >&2
    exit 64
}

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

number_of_words=$(wc -w <<< "${args[@]}")
if [ ${number_of_words} -le 0 ];then
    numpy_versions=( 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 )
elif [ ${number_of_words} -le 1 ];then
    numpy_versions=${args[0]}
else
    numpy_versions=${args[0]}
    export GIT_PYDUSA_VERSION=${args[1]}
fi


for v in ${numpy_versions[@]};do
    conda build "${RECIPES_DIR}/pydusa" --numpy ${v} ${opts[@]}
done
