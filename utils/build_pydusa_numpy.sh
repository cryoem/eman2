#!/usr/bin/env bash

# Builds user-specified numpy version from optionally-specified branch

function print_usage(){
    printf "\e\033[35m\n  Usage: $(basename ${0}) %s %s %s\n\n\033[0m" "{1.5|...|1.13}" "[branch-to-checkout]" "[--conda-build-opt1] [--conda-build-opt2] ..." >&2
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

case ${#args[@]} in
    1|2) numpy_versions=${args[0]}
        
        export GIT_PYDUSA_BRANCH=${args[1]:-"v20180115"}
        ;;

    *) print_usage
       ;;
esac

for v in ${numpy_versions[@]};do
    conda build "${RECIPES_DIR}/pydusa" --numpy ${v} ${opts[@]}
done
