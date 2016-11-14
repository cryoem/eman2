#!/usr/bin/env python

#
# Author: Steven Ludtke, 4/17/15 (sludtke@bcm.edu)
# Copyright (c) 2000-2010 Baylor College of Medicine
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
import time
import random


out=file("fftspeed.txt","w")
# for 1-D python loop is too slow, so we have to get tricky
#for size in xrange(20,256,2):
#	img=EMData(size,10000,1)
	#img2=EMData(size,10000,1)
	#reps=1000
#
	#img.process_inplace("math.addnoise",{"noise":1.0})
	#img2.process_inplace("math.addnoise",{"noise":1.0})
#
	#t0=time.time()
	#x=img.calc_ccfx(img2,0,9999,1)
	#t1=time.time()
#
	#rslt="%d\t%1.2f\t%d"%(size,t1-t0,1)
	#out.write(rslt+"\n")
	#print rslt

for dim in xrange(2,4):
	times=[]
	for size in xrange(20,514,2):
		if dim==2 : 
			img=EMData(size,size,1)
			reps=int(200000/size)
		elif dim==3 : 
			img=EMData(size,size,size)
			reps=int(2000/(size*size))+1
	

		img.process_inplace("math.addnoise",{"noise":1.0})
	
		t0=time.time()
		for r in xrange(reps):
			x=img.do_fft()
			
		t1=time.time()

		times.append(((t1-t0)/reps,size))
		rslt="%d\t%1.3f\t%d"%(size,10000*(t1-t0)/reps,dim)
		out.write(rslt+"\n")
		out.flush()
		print rslt,"\t",reps

	times.reverse()
	times2=[times[0]]
	for i in times[1:]:
		if i<times2[-1] : times2.append(i)
	
	times2.reverse()

	for t in times2: print t[1],
	print "\n"

	out.write(str(times2)+"\n")

