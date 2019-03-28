#!/bin/bash

this_dir=$(dirname ${0})
bin_dir=$(realpath ${this_dir}/../bin)
echo ${this_dir}
echo ${bin_dir}
for file_name in ${bin_dir}/*.py
do
    sed -i.bkp "s|#!/usr/bin/env python|#!${bin_dir}/python|g" ${file_name}
    #sed "s|#!/usr/bin/env python|#!${bin_dir}/python|g" ${file_name}
done
