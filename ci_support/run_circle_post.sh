#!/usr/bin/env bash

# Run tests
export PATH="$HOME/miniconda2/bin/:$HOME/EMAN2/bin:$PATH"
export PYTHONPATH="$HOME/EMAN2/lib:$HOME/EMAN2/bin:$PYTHONPATH"

e2version.py
e2speedtest.py
