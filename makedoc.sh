#!/bin/sh

which doxygen 2>/dev/null 1>/dev/null
if test ! $? = 0; then
    echo
    echo "Error: doxygen is not found. Please install doxgen first"
    echo
    exit 1
fi

doxygen  doc/Doxyfile
