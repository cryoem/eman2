#!/usr/bin/env bash

# Builds EMAN against user-specified NumPy version(s)

if [ $# -lt 1 ];then
    echo
    echo -e '\033[35m'"  Usage: $(basename ${0}) [NumPy version(s)]"'\033[0m'
    echo
    exit 1
fi

source activate root

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

src_dir="${CONDA_PREFIX}/github.com/cryoem/eman2"
if [ ! -d "${src_dir}/.git" ];then
    mkdir -p "${src_dir}" && cd "${src_dir}/.."
    git clone https://github.com/cryoem/eman2.git --branch master
    cd -
else
    cd "${src_dir}"
    git checkout master -f
    git pull --rebase
    cd -
fi

for v in ${args[@]};do
    conda build "${src_dir}/recipes/eman" --numpy ${v} -c cryoem -c defaults -c conda-forge ${opts[@]}
done
