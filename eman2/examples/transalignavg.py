#!/usr/bin/env python

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

# transalignavg 12/02/2005  Steven Ludtke
# This will read a series of images, translationally align them, average them
# together, and optionally iterate. Translational alignment only.
# transalignavg.py <infile> <dot threshold> <iterations> <background> <gamma>

from EMAN2 import *
import sys
from math import *


n=EMUtil.get_image_count(sys.argv[1])
thr=float(sys.argv[2])

if len(sys.argv)>3: iter=int(sys.argv[3])
else: iter=1

darkref=None
if len(sys.argv)>4 :
	darkref=EMData()
	darkref.read_image(sys.argv[4],0)
else: darkref=None

if len(sys.argv)>5 : gamma=float(sys.argv[5])
else : gamma=0

# Pass 1, sequential alignment to average
avg=EMData()
avg.read_image(sys.argv[1],0)
ref0=avg.process("filter.lowpass.gauss",{"cutoff_abs":.15})
nx=avg.get_xsize()
if darkref :
	darkref.process_inplace("normalize.toimage",{"to":avg})
	avg-=darkref
avg-=avg.get_edge_mean()
#if gamma : avg.process_inplace("math.pow",{"pow":gamma})
sum=1
for i in range(1,n):
	aa=EMData()
	aa.read_image(sys.argv[1],i)
	if darkref :
		darkref.process_inplace("normalize.toimage",{"to":a})
		aa-=darkref
	aa-=aa.get_edge_mean()
	
	a=aa.process("filter.lowpass.gauss",{"cutoff_abs":.1})
	
	b=a.align("rotate_translate",ref0,{"nozero":1,"maxshift":nx/2},"dot",{"negative":1})
#	b=a.align("translational",ref0,{"nozero":1,"maxshift":nx/2},"dot",{"negative":1})
#	b.set_attr("rotational",0.0)
	dot=b.cmp("dot",ref0,{"negative":0,"normalize":1})
	#bdot=0
	#bang=0
	#for angi in range(-8,8):
		#ang=angi/10.0+.05
		#a=aa.copy()
		#a.rotate(ang,0,0)
	##	if gamma : a.process_inplace("math.pow",{"pow":gamma})
	##	b=a.align("translational",ref0,{"nozero":1,"maxshift":nx/2})
		#b=a.align("translational",ref0,{"nozero":1,"maxshift":nx/2})
		#dot=b.cmp("dot",ref0,{"negative":0,"normalize":1})
		#if dot>bdot : 
			#bdot=dot
			#bang=ang
			#bdx=b.get_attr("translational.dx")
			#bdy=b.get_attr("translational.dy")
	
	a=aa.copy()
	a.rotate(b.get_attr("rotational"),0,0)
	a.translate(b.get_attr("translational.dx"),b.get_attr("translational.dy"),0)
#	b=a.align("translational",ref0,{"nozero":1,"maxshift":nx/2})
	print "%4d. %3d\t%3d\t%1.2f\t%1.4f"%(i,b.get_attr("translational.dx"),b.get_attr("translational.dy"),b.get_attr("rotational"),dot)
	if dot>thr : 
		avg=avg+a
		sum+=1
	
print "%d/%d used"%(sum,n)
avg-=avg.get_attr("minimum")
avg/=avg.get_attr("maximum")
#avg.process_inplace("math.pow",{"pow":gamma})
avg.write_image("avg.mrc")
display(avg)
