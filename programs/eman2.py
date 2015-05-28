#!/usr/bin/env python
#
# Author: Steve Ludtke 5/7/15 (sludtke@bcm.edu)
# Copyright (c) 2000-2015 Baylor College of Medicine

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

helpstring="""Hello and welcome to EMAN2 !

EMAN2 is NOT a single central graphical program, but includes a large set of programs both graphical and command-line. 

If you are just getting started with EMAN2, here are some tips:

1) If you are just interested in using the image browser and other display programs, simply run e2display.py
2) If you want to use EMAN2 for single particle reconstruction, make an empty folder, 'cd' into that folder, then run e2projectmanager.py
3) If you want to use EMAN2 for file format conversions, simple image processing, etc. here are some programs to run to get started:
	* e2proc2d.py --help
	* e2proc3d.py --help
	* e2help.py processors
	* e2help.py processors -v 2
	* e2display.py

4) Please consider going through a simple tutorial available at: http://blake.bcm.edu/emanwiki/EMAN2/Tutorials
"""

from EMAN2 import *
from e2version import EMANVERSION,CVSDATESTAMP
from emapplication import EMApp
import os, json, re, glob, signal
import subprocess

helpstring+="\n\nYou are currently running %s (%s)"%(EMANVERSION,CVSDATESTAMP[6:-2])

try:
	if os.getenv("DISPLAY")==None : raise Exception
	from PyQt4 import QtCore, QtGui
	from PyQt4.QtCore import Qt
except:
	print helpstring
	raw_input("Please press <enter> to exit")
	exit()	

app = EMApp()
QtGui.QMessageBox.warning(None,"Welcome to EMAN2!","<b><big>"+helpstring.replace("\n","<br>")+"</big></b>")

