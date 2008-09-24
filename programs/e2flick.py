#!/usr/bin/env python

#
# Author:David Woolford 09/23/2008 (woolford@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#


#from EMAN2 import *
#from emimage import EMImage
#from emimageutil import EMParentWin
#import sys
#from optparse import OptionParser
#from PyQt4 import QtCore, QtGui, QtOpenGL
#from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU,GLUT

from emimagemxrotor import *
from optparse import OptionParser
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image file> ...
	
	An interactive interface for viewing and cleaning large sets of images.

	Uses a local database to remember which images have previously been deleted.

	Execute e2flick.py with no function arguments to see an example.
	
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	(options, args) = parser.parse_args()

	logid=E2init(sys.argv)
	app = QtGui.QApplication(sys.argv)
	window = EMImageMXRotor()
	if len(sys.argv)==1 : 
		data = []
		for i in range(500):
			#if i == 0: idx = 0
			#else: idx = i%64
			#e = EMData(64,64)
			#e.set_size(64,64,1)
			#e.to_zero()
			#e.add(sin( (i/10.0) % (pi/2)))
			data.append(test_image(Util.get_irand(0,3)))
			#data.append(e)
		
		
		window.setData(data)
	else :
		window.set_image_file_name(sys.argv[1])
	window2=EMParentWin(window)
	window2.resize(*window.get_optimal_size())
	window2.show()
	sys.exit(app.exec_())
	
	E2end(logid)
	
# If executed as a program
if __name__ == '__main__':
	main()	

