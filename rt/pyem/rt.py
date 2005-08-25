#!/usr/bin/env python

"""regression test: run all the python unittest under current directory."""

import os
import glob

for name in glob.glob('test_*.py'):
    os.system('python %s'%name)
    
    