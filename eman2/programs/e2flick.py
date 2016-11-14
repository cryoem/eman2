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
import sys
import os
#from optparse import OptionParser
#from PyQt4 import QtCore, QtGui, QtOpenGL
#from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU,GLUT

#from emimagemxrotor import * #emimagemxrotor was deprecated and removed from CVS
#from emapplication import EMStandAloneApplication #EMStandAloneApplication was deprecated and removed from CVS
from optparse import OptionParser

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image file> ...
	
	An interactive interface for viewing and cleaning large sets of images.

	Uses a local database to remember which images have previously been deleted.

	Execute e2flick.py with no function arguments to see an example.
	
"""

	print "WARNING: This program is currently broken. We intend to ressurect it in future."
	sys.exit(1)

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	logid=E2init(sys.argv,options.ppid)
	
	
	app = EMStandAloneApplication()
	window = EMImageMXRotorModule(application=app)
	
	if len(sys.argv)==1 : 
		data = []
		for i in range(0,500):
			e = test_image(Util.get_irand(0,9))
			if ( Util.get_irand(0,4) == 0):	e.set_attr("excluded",True)
			data.append(e)
			
		window.set_data(data)
	else :
		window.set_image_file_name(sys.argv[1])
		
	app.show()
	app.execute()
	
	E2end(logid)
	
# If executed as a program
if __name__ == '__main__':
	main()	

