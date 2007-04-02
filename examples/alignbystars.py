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
import sys
from optparse import OptionParser
from pprint import pprint

def findstars(img,scale=1.0):
	"""This will find localized peaks to subpixel accuracy and return a list of
	x,y,peak,rad_gyr"""
	imc=img.copy()
	thr=imc.get_attr("mean")+imc.get_attr("sigma")*3.0
	imc.process_inplace("filter.lowpass.gauss",{"sigma":0.3})
	imc.process_inplace("eman1.mask.onlypeaks",{"npeaks":0})
	peaks=imc.calc_highest_locations(thr)
	
	ret=[]
	for i in peaks:
		c=img.get_clip(Region(i.x-5,i.y-5,11,11))
		
		# this should remove 'hot' ccd pixels
		if c.get_value_at(5,5)>c.get_value_at(4,5)+c.get_value_at(6,5)+c.get_value_at(5,4)+c.get_value_at(5,6) : continue
		
		# find the center of mass of each peak
		cg=c.cog()[:3]
		ret.append(((i.x+cg[0]-img.get_xsize()/2)*scale,(i.y+cg[1]-img.get_ysize()/2)*scale,i.value,cg[2]))	# x,y,peak,rad
	
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

def l2pa(a):
	"""Convert a list, as output by findstars, into a PointArray"""
	r=PointArray(len(a))
	for i,j in enumerate(a):
		r.set_vector_at(i,Vec3f(j[0],j[1],0),j[2])
	
	return r

def alignstars(a,b):
	"""This will take two lists of x,y,peak,rad_gyr and align them in 2-d"""
	a=l2pa(a)
	b=l2pa(b)
#	print a.align_trans_2d(b)
	print a.align_2d(b)
	print a.align_trans_2d(b,1,0,0)
#	print centerofstars(a),centerofstars(b)


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image 1> <image 2> ...
	
Finds isolated spots in the image and uses them as a basis for alignment"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--scale", action="store_true", help="Restrict axes for c/d symmetric objects", default=False)
	parser.add_option("--scale", "-S", type="float", help="scale factor", default=1.0)
	parser.add_option("--out", "-o", type="string", help="Output filespec", default="out.mrc")
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("At least 2 inputs required")
	scale=options.scale

	pats=[]
	for j,i in enumerate(args):
		img=EMData(i,0)
		pats.append(l2pa(findstars(img,scale)))
		print " %04d\r"%j,
		sys.stdout.flush()
		
	avg=EMData(args[0],0)
	avg.process_inplace("eman1.normalize.edgemean",{})
	nx=avg.get_xsize()
	ny=avg.get_ysize()
	if scale!=1.0 : 
		avg=avg.get_clip(Region(-nx*(scale-1)/2,-ny*(scale-1)/2,nx*scale,ny*scale))
		avg.scale(1.0/scale)
	avgn=avg.copy()
	avgn.to_one()
		
	for i in range(1,len(pats)):
		xf=pats[i].align_2d(pats[0])
		print "%d. %6.2f\t%6.2f"%(i,xf.get_pretrans().at(0),xf.get_pretrans().at(1))
		img=EMData(args[i],0)
		img.process_inplace("eman1.normalize.edgemean",{})
		if scale!=1.0 : 
			img=img.get_clip(Region(-nx*(scale-1)/2,-ny*(scale-1)/2,nx*scale,ny*scale))
			img.scale(1.0/scale)
		imgn=img.copy()
		imgn.to_one()
		img.rotate_translate(xf)
		imgn.rotate_translate(xf)
		avg+=img
		avgn+=imgn
	
	avg/=avgn
	avg.write_image(options.out)
	#print str(xf)
	#print xf.get_scale()
	#print xf.get_center()
	#print xf.get_posttrans()
	#print xf.get_rotation(EULER_EMAN)

	

	
if __name__ == "__main__":
    main()
