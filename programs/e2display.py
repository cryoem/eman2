#!/usr/bin/env python

#
# Author: Steven Ludtke 11/09/2006 (sludtke@bcm.edu)
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

from EMAN2 import *
from emimage import EMImage
import sys
from optparse import OptionParser
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU,GLUT
#from valslider import ValSlider
#from math import *
#import numpy
#from emimageutil import ImgHistogram
#from weakref import WeakKeyDictionary

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image file> ...
	
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
#	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: circle, ref, grid",default=[])
#	parser.add_option("--threshold","-T",type="float",help="Threshold for keeping particles. 0-4, 0 excludes all, 4 keeps all.",default=2.0)
#	parser.add_option("--maxbad","-M",type="int",help="Maximumum number of unassigned helices",default=2)
#	parser.add_option("--minhelix","-H",type="int",help="Minimum residues in a helix",default=6)
#	parser.add_option("--apix","-P",type="float",help="A/Pixel",default=1.0)
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	logid=E2init(sys.argv)
	display(args)
	E2end(logid)

def display(img):
	if isinstance(img,str) : img=[img]
	app = QtGui.QApplication(sys.argv)
	
	win=[]
	for f in img:
		a=EMData.read_images(f)
		if len(a)==1 : a=a[0]
		w=EMImage(a)
		w.setWindowTitle("EMImage (%s)"%f)
		w.show()
		win.append(w)
	
	sys.exit(app.exec_())

# If executed as a program
if __name__ == '__main__':
	main()	
