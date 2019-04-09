#!/usr/bin/env python
from __future__ import print_function

import os
import platform

spreal = os.path.join(os.path.abspath(os.path.dirname(__file__)), "sp_real.py")
ipython = os.path.join(os.path.abspath(os.path.dirname(__file__)), "ipython")
try:
	if platform.system()=="Linux" and os.getenv("DISPLAY")==None: raise Exception
	os.execlp(ipython,"ipython","-i","--gui=qt5",spreal)
except:
	print("Warning: No DISPLAY available, running in non-GUI mode.")
	os.execlp(ipython,"ipython","-i",spreal)
