#!/bin/sh

# Generate Doxygen documentation under eman2/doc
# Usage: makedoc.sh

which doxygen 2>/dev/null 1>/dev/null
if test ! $? = 0; then
    echo
    echo "Error: 'doxygen' is not found. Please install 'doxgen' first"
    echo
    exit 1
fi

echo -n "Start to generate Doxygen documentation. Be patient ... "
doxygen
echo "Done"

set -xe

rm -rf doxygen_html
mv -f html doxygen_html
