#!/usr/bin/env bash

source ci_support/pre_build.sh

cmake $src_dir
make
make install

source $src_dir/ci_support/post_build.sh
