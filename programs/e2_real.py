#! /usr/bin/env python
#
# Do not run this script directly. You could run this as 
# ipython --gui=qt -i e2_real.py
# otherwise it is normally run via the e2.py script
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

import EMAN2
import os
from EMAN2 import *

failed=False
try:
	if get_platform()=="Linux" and os.getenv("DISPLAY")==None: raise Exception

	from PyQt4 import QtCore, QtGui, QtOpenGL
	from emapplication import EMApp
	import IPython.lib.inputhook


	app=EMApp()
	IPython.lib.inputhook.enable_qt4(app)

	from emimage import image_update

	def ipy_on_timer():
		image_update()

	ipytimer = QtCore.QTimer()
	ipytimer.timeout.connect(ipy_on_timer)
	ipytimer.start(200)
	
	EMAN2.GUIMode=True
	EMAN2.app=app
except:
	failed=True


print """Welcome to the interactive EMAN2 Python interface, provided by ipython

NOTE: that you should NOT be running this program if your intent is to run other
  EMAN2 programs. To do that, quit e2.py and just type the name of the EMAN2 program
  directly at the system command line. This interface is for people who know some
  Python programming."""
  
if failed : print "\nWarning: Failed to initialize Qt mode. Running in fallback mode with non-interactive graphics."
