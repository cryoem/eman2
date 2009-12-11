#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/28/2008 (sludtke@bcm.edu)
# Copyright (c) 2000-2008 Baylor College of Medicine
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

from EMAN2 import *
import sys
from math import *
from optparse import OptionParser
import os

def main(argv,app=None) :

	ome = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input stack/image> ...
	
Various CTF-related operations on images."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--size",type="string",help="nx,ny",default="512,256")
	parser.add_option("--origin",type="string",help="x,y",default="256,256")
	parser.add_option("--xrange",type="string",help="x0,x1,dx",default="0,1,1")
	parser.add_option("--yrange",type="string",help="y0,y1,dy",default="0,1,1")
#	parser.add_option("--zrange",type="string",help="z0,zx,z1",default="0,1,1")
	parser.add_option("--wavelen",type="float",help="Wavelength in pixels",default=16.0)
	parser.add_option("--phase",type="float",help="Phase shift",default=0)
	parser.add_option("--plane",action="store_true",help="Make a plane wave rather than a circular wave",default=False)

	(options, args) = parser.parse_args()
	
	size=eval(options.size)
	origin=eval(options.origin)
	
	for i in range(8):
		nsrc=0
		sm=EMData(*size)
		t=EMData(*size)
		for x in range(*eval(options.xrange)):
			for y in range(*eval(options.yrange)):
				if options.plane:
					t.process_inplace("testimage.sinewave",{"wavelength":options.wavelen,"phase":options.phase+2*pi*x/(options.wavelen)-0.25*pi*i,"axis":"x"})
				else :
					t.process_inplace("testimage.sphericalwave",{"wavelength":options.wavelen,"phase":options.phase+2*pi*x/(options.wavelen)-0.25*pi*i,"x":x+origin[0],"y":y+origin[1],"z":0})
				sm+=t
				nsrc+=1
		if options.plane :
			sm["render_min"]=sm["minimum"]
			sm["render_max"]=sm["maximum"]
		else :
			sm["render_min"]=-10.0*nsrc/256.0
			sm["render_max"]=10.0*nsrc/256.0
		sm.write_image("img.%02d.png"%i)
		
	os.system("ffmpeg -i img.%02d.png out.mov")
#	os.system("animate -delay 2 img.?.png")

if __name__ == "__main__":
   	main(sys.argv)
