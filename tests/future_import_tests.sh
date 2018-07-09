#!/usr/bin/env bash

future_imports=( 
                print_function
               )

py_files=$( find . -not -empty -name '*.py' )

failed_cases=()
for f in ${py_files[@]};do
    for imp in ${future_imports[@]};do
        if ! grep -q "from __future__ import ${imp}" ${f};then
            failed_cases+=("$f")
            echo "  ${f} is missing \"from __future__ import ${imp}"\"
        fi
    done
done

echo

if [ ${#failed_cases[@]} -ne 0 ];then
    exit 1
fi
