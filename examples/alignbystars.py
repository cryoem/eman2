#!/bin/env python
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

from EMAN2 import *
from emimage import *
import time
from optparse import OptionParser

def findstars(img):
	"""This will find localized peaks to subpixel accuracy and return a list of
	x,y,peak,rad_gyr"""
	imc=img.copy()
	thr=imc.get_attr("mean")+imc.get_attr("sigma")*8.0
	imc.process_inplace("filter.lowpass.gauss",{"sigma":0.3})
	imc.process_inplace("eman1.mask.onlypeaks",{"npeaks":0})
	peaks=imc.calc_highest_locations(thr)
	
	ret=[]
	for i in peaks:
		c=img.get_clip(Region(i.x-5,i.y-5,11,11))
		cg=c.cog()[:3]
		ret.append((i.x+cg[0],i.y+cg[1],i.value,cg[2]))	# x,y,peak,rad
	
	return ret

def centerofstars(a):
	"""takes a list of x,y,peak,rad_gyr and returns a 'center of mass'"""
	
	x=0
	y=0
	s=0
	for i in a:
		x+=i[0]
		y+=i[1]
		s+=1
	
	return (x/s,y/s)

def alignstars(a,b):
	"""This will take two lists of x,y,peak,rad_gyr and align them in 2-d"""
	
	print centerofstars(a),centerofstars(b)


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image 1> <image 2> ...
	
Finds isolated spots in the image and uses them as a basis for alignment"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--csym", action="store_true", help="Restrict axes for c/d symmetric objects", default=False)
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("At least 2 inputs required")

	pats=[]
	for i in args:
		img=EMData(i,0)
		pats.append(findstars(img))
	
	print alignstars(pats[0],pats[1])
	print alignstars(pats[1],pats[2])
	print alignstars(pats[2],pats[3])
	
if __name__ == "__main__":
    main()
