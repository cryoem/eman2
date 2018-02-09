#!/usr/bin/env bash

set -x

cmd "/C ${WORKSPACE//\//\\}\\tests\\test_binary_installation.bat ${INSTALLERS_DIR//\//\\} eman2.win.exe "
