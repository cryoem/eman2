#!/bin/sh

# to purify C++ code using GNU indent.
# usage: ./purifycode.sh

which indent 2>/dev/null 1>/dev/null
if test ! $? = 0; then
    echo
    echo "Error: 'indent' is not found. Please install 'indent' first"
    echo
    exit 1
fi

DIRS="plugins"
INDENT_F="-br -nce -i4  -npcs -nprs  -npsl -l100 -ts4"

for dir1 in $DIRS; do
	cd $dir1
	files=`/bin/ls *.h *.cpp`

	for file1 in $files; do
		indent $INDENT_F $file1
		sed -i s/"const const"/"const"/g $file1
	done
	cd ..
done


