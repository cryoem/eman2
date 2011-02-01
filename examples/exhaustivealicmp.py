#!/usr/bin/env python

#
# Author: Steve Ludtke, 2/1/2011 (stevel@bcm.edu)
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

# This program will rotate/translate an input image and compare to a reference
# using a selected comparator


from EMAN2 import *
from optparse import OptionParser
from math import *
from copy import deepcopy
import os
import sys


def main():
        progname = os.path.basename(sys.argv[0])
        usage = """%prog <file1> <N1> <file2> <N2> <outfile> [options]
	
	Using the selected comparator this program will rotate/translate an image exhaustively over a range
	and write the similarity value to an output volume. (technically it writes the similarity * -1.0)

        """

        parser = OptionParser(usage=usage,version=EMANVERSION)

        parser.add_option("--xy0",type="float",help="How far to shift x/y from 0,0. Default = 20.0",default=20.0)
        parser.add_option("--dxy",type="float",help="Step (in pixels) for x/y translation. Default 1.0",default=1.0)
	parser.add_option("--dalpha",type="float",help="Angular step (in degrees). Default=3.0",default=5.0)
	parser.add_option("--cmp",type="string",help="Comparator to use. Default=ccc",default="ccc")

        (options, args) = parser.parse_args()

	cmpopt=parsemodopt(options.cmp)
	nxy=int(options.xy0/options.dxy)*2+1
	nz=int(360.0/options.dalpha)-1

	im1=EMData(args[0],int(args[1]))
	im2=EMData(args[2],int(args[3]))

	ali=im2.align("rotate_translate_flip",im1,{},cmpopt[0],cmpopt[1])
	a2=ali["xform.align2d"]
	print a2.inverse()
	if a2.get_mirror():
		im2.process_inplace("xform.flip",{"axis":"x"})
		ali=im2.align("rotate_translate_flip",im1,{},cmpopt[0],cmpopt[1])
		a2=ali["xform.align2d"]
		print a2.inverse()
	
	ali=im2.align("refine",im1,{"xform.align2d":a2},cmpopt[0],cmpopt[1])
	a2=ali["xform.align2d"]
	print a2.inverse(),"\n\n"

	output=EMData(nxy,nxy,nz)
	output.to_zero()

	best=(-1.0e8,0,0,0)
	alpha=0.0
	k=0
	while alpha<360.0:
		y=-options.xy0
		j=0
		print "\r alpha: ",alpha,"       ",
		sys.stdout.flush()
		while y<options.xy0:
			x=-options.xy0
			i=0
			while x<options.xy0:
				im2a=im2.process("xform",{"transform":Transform({"type":"2d","tx":x,"ty":y,"alpha":alpha}).inverse()})
				val=-im1.cmp(cmpopt[0],im2a,cmpopt[1])
				output[i,j,k]=val
				if val>best[0] : best=(val,x,y,alpha)
				x+=options.dxy
				i+=1
			y+=options.dxy
			j+=1
		alpha+=options.dalpha
		k+=1

	output.write_image(args[4],0)
	print "Best :",best

	best2=(-1.0e8,0,0,0)
	output=EMData(41,41,41)
	output.to_zero()
	for k in range(-20,21):
		print "\r alpha: ",k,"       ",
		sys.stdout.flush()
		for j in range(-20,21):
			for i in range(-20,21):
				im2a=im2.process("xform",{"transform":Transform({"type":"2d","tx":best[1]+i*0.1,"ty":best[2]+j*0.1,"alpha":best[3]+k*0.2}).inverse()})
				val=-im1.cmp(cmpopt[0],im2a,cmpopt[1])
				if val>best2[0] : best2=(val,best[1]+i*0.1,best[2]+j*0.1,best[3]+k*0.2)
				output[i+20,j+20,k+20]=val
				
				
	output.write_image("o_"+args[4],0)
	print "Best :",best2



if __name__ == "__main__":
    main()


