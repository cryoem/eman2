#!/usr/bin/env python
from __future__ import print_function

import os

sxreal = os.path.join(os.path.abspath(os.path.dirname(__file__)), "sx_real.py")
ipython = os.path.join(os.path.abspath(os.path.dirname(__file__)), "ipython")
os.execlp(ipython,ipython,"-i",sxreal)
