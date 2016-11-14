#!/bin/sh

# to generate Doxygen documentation under eman2/doc
# usage: makedoc.sh

which doxygen 2>/dev/null 1>/dev/null
if test ! $? = 0; then
    echo
    echo "Error: 'doxygen' is not found. Please install 'doxgen' first"
    echo
    exit 1
fi

echo -n "Start to generate Doxygen documentation. Be patient ... "
cd ${EMAN2DIR}/src/eman2
doxygen  doc/doxygen/Doxyfile
echo "Done"

rm -rf ${EMAN2DIR}/doc/doxygen_html
mv -f ${EMAN2DIR}/src/eman2/doc/html ${EMAN2DIR}/doc/doxygen_html

#echo "Documentation is at $PWD/doc/html/index.html"
