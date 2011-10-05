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

# This program will find the highest density point, then build a high density path until
# one putative subunit is filled.

from EMAN2 import *
from math import *
from bisect import insort
from random import random

def main() :
	global threshold
	
        usage = """prog [options] <input> <output>


	This program attempts to extract one subunit from a volume by starting at
	a high point then iteratively expanding until a subunit has been located.
"""

	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--sym", "-S", type=str, help="Symmetry", default="c4")
	parser.add_argument("--thr", "-T", type=float, help="Isosurface threshold", default=1.0)
	parser.add_argument("--random","-R",action="store_true", help="Randomize the starting location", default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	print """WARNING: Experimental program. Contact sludtke@bcm.edu before relying on its results."""
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")

	data=EMData()
	data.read_image(args[0])
	data.process_inplace("threshold.belowtozero",{"minval":options.thr})
	
	totest=[]
	for x in range(-1,2):
		for y in range(-1,2):
			for z in range(-1,2):
				totest.append((x,y,z))
	
	# need to replace this with real symmetry code, for now just C supported
	# ths is a symmetry operator for 1-fold rotation
	if options.sym[0]=='c': csym=int(options.sym[1:])
	symop=Transform()
	symop.set_rotation(symop.EulerType.EMAN,2.0*pi/csym,0,0)
	
	# now seed the iterative process with the highest value in the map
	cur=data.calc_max_location()
	
	if (options.random) :
		cur=(int(floor(random()*data.get_xsize())),int(floor(random()*data.get_ysize())),int(floor(random()*data.get_zsize())))
		while data.get_value_at(cur[0],cur[1],cur[2])==0 :
			cur=(int(floor(random()*data.get_xsize())),int(floor(random()*data.get_ysize())),int(floor(random()*data.get_zsize())))
	
	plist=[]
	nvox=0
	while (1):
		# invert the voxel, and zero the symmetric ones
		data.set_value_at(cur[0],cur[1],cur[2],-data.get_value_at(cur[0],cur[1],cur[2]))
		c2=(cur[0]-data.get_xsize()/2,cur[1]-data.get_ysize()/2,cur[2]-data.get_zsize()/2)
		for a in range(1,csym):
			c2=(c2*symop).as_list()
			c3=(int(c2[0]+data.get_xsize()/2),int(c2[1]+data.get_ysize()/2),int(c2[2]+data.get_zsize()/2))
			for l in totest:
				data.set_value_at(c3[0]+l[0],c3[1]+l[1],c3[2]+l[2],0.0)
		
		c2=cur
		# insert any neighboring pixels into our list if they are above the threshold
		for l in totest:
			v=data.get_value_at(c2[0]+l[0],c2[1]+l[1],c2[2]+l[2])
			if v>options.thr:
				ti=(v,c2[0]+l[0],c2[1]+l[1],c2[2]+l[2])
				ii=bisect_left(plist,ti)					# find the insertion spot
				try: 
					if plist[ii]==ti : continue					# skip it if it's already there
				except: pass
				plist.insert(ii,ti)
		
		while (1):
			if len(plist)==0: break
			cur=plist.pop()[1:]
#			print cur,data.get_value_at(cur[0],cur[1],cur[2])
			if data.get_value_at(cur[0],cur[1],cur[2])>options.thr: break
			
		if len(plist)==0: break
			
		nvox+=1
		if nvox%100==0 : print nvox,len(plist)
		
	data*=-1.0
	data.process_inplace("threshold.belowtozero",{"minval":options.thr})
	data.update()
	data.write_image(args[1])
		
def ilist(d,l,x,y,z):
	global threshold
	
	v=d.get_value_at(x,y,z)
	if v>threshold : insort(l,(v,x,y,z))
	
if __name__ == "__main__":
    main()
