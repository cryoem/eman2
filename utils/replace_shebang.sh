#!/bin/bash

this_dir=$(dirname ${0})
bin_dir=$(realpath ${this_dir}/../bin)
echo 'Script directory:' ${this_dir}
echo 'Bin directory:' ${bin_dir}
for file_name in ${bin_dir}/*.py
do
    echo 'Change shebang in file:' ${file_name}
    sed -i.bkp "s|^#![ \t]*/usr/bin/env[ \t]*python.*$|#!${bin_dir}/python|g" ${file_name}
    chmod +x ${file_name}
    chmod -x ${file_name}.bkp
done
