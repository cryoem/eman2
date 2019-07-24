#!/usr/bin/env bash

set -x

INSTALLERS_DIR="${HOME_DIR}/workspace/${JOB_NAME}-installers"

cmd "/C ${WORKSPACE//\//\\}\\tests\\test_binary_installation.bat ${INSTALLERS_DIR//\//\\}\eman2.win.exe ${INSTALLERS_DIR//\//\\}\eman2-binary-test "
