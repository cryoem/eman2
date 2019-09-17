#!/usr/bin/env bash

set -x

case $# in
    1) suffix=""
       dir="$1"
       ;;

    2) suffix="$1"
       dir="$2"
       ;;

    *) echo "Need 1 or 2 arguments!"
       ;;
esac

cmd "/C ${WORKSPACE//\//\\}\\tests\\test_binary_installation.bat ${INSTALLERS_DIR//\//\\}\eman2${suffix}.win.exe ${INSTALLERS_DIR//\//\\}\\${dir} "
