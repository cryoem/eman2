#!/usr/bin/env python

#
# Author: Steven Ludtke, 01/03/07 (sludtke@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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

# e2stacktoanim.py  01/03/07  Steven Ludtke
# This program will convert a 2D image stack into a GIF animation


from EMAN2 import *
from math import *
import os
import sys

def stacktoanim(stack,outpath,ntk):
	"""Takes an input list of images and and output pathname. Converts each image to a standard 2D image
	format, then produces a GIF animation. Requires functional ImageMagick installation."""
	for i in range(ntk+1):
		im=stack[i]
		im.set_attr("render_min",im.get_attr("mean")-im.get_attr("sigma")*3.0)
		im.set_attr("render_max",im.get_attr("mean")+im.get_attr("sigma")*3.0)
		im.write_image("tmp_img-%03d.pgm"%i)
		print "%d. %1.3f - %1.3f"%(i,im.get_attr("render_min"),im.get_attr("render_max"))
	os.system("convert -delay 10 tmp_img-*.pgm %s "%outpath)
	
	for i in range(ntk+1):
		os.unlink("tmp_img-%03d.pgm"%i)

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] input_stack.hed output.(gif:pnm)
	
	Converts a 2D image stack into a GIF/PNM animation using ImageMagick (which must be installed for this to work)"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--scale", "-S", type=float, help="Scale factor",default=1.0)
	parser.add_argument("--pingpong",action="store_true",default=False,help="Cycle through the sequence forwards then backwards")
	parser.add_argument("--last","-M", type=int, help="Number of last image to use",default=0)
#	parser.add_argument("--mode",type=str,help="centering mode 'modeshift', 'censym' or 'region,<x>,<y>,<clipsize>,<alisize>",default="censym")
#	parser.add_argument("--twopass",action="store_true",default=False,help="Skip automatic tilt axis location, use fixed angle from x")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	
	a=EMData.read_images(args[0])
	if options.last>0 and options.last<len(a): ntk=options.last
	else : ntk=len(a)-1
	
	# rescale images if requested
	if options.scale!=1.0 :
		olds=(a[0].get_xsize(),a[0].get_ysize())
		news=(int(olds[0]*options.scale),int(olds[1]*options.scale))
		for i in range(ntk+1):
			if options.scale<1.0: a[i].scale(options.scale)
			a[i]=a[i].get_clip(Region((olds[0]-news[0])/2.0,(olds[1]-news[1])/2.0,news[0],news[1]))
			if options.scale>1.0: a[i].scale(options.scale)

	if options.pingpong :
		for i in range(ntk-1,-1,-1):
			a.append(a[i])
			
	stacktoanim(a,args[1],ntk)

if __name__ == "__main__":
    main()
