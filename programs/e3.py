#!/usr/bin/env python
import os
import platform
import sys

# Change default behavior to execute jupyter-lab
if len(sys.argv)>1 and sys.argv[1]=="--ipython":
	e2real = os.path.join(os.path.abspath(os.path.dirname(__file__)), "e2_real.py")
	try:
		if platform.system()=="Linux" and os.getenv("DISPLAY")==None: raise Exception
		os.execlp("ipython","ipython","-i","--gui=qt5",e2real)
		sys.exit(0)
	except:
		print("Warning: No DISPLAY available, running in non-GUI mode.")
		os.execlp("ipython","ipython","-i",e2real)


os.system("jupyter-lab")
