#!/usr/bin/env python

import sys
import os
from optparse import OptionParser
from EMAN2 import *
from emimage import *
import pyshed
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
	


if __name__ == "__main__":
	app = get_app()
	window = pyshed.Shell()
	window.show()
	
	sys.exit(app.exec_())
