#!/usr/bin/env python

# This program will find the highest density point, then build a high density path until
# one putative subunit is filled.

from EMAN2 import *
from math import *
from bisect import insort
from optparse import OptionParser
from random import random

def main() :
	global threshold
	
	usage=""
	
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--sym", "-S", type="string", help="Symmetry", default="c4")
	parser.add_option("--thr", "-T", type="float", help="Isosurface threshold", default=1.0)
	parser.add_option("--random","-R",action="store_true", help="Randomize the starting location", default=False)
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")

	data=EMData()
	data.read_image(args[0])
	data.filter("threshold.belowtozero",{"minval":options.thr})
	
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
	data.filter("threshold.belowtozero",{"minval":options.thr})
	data.update()
	data.write_image(args[1])
		
def ilist(d,l,x,y,z):
	global threshold
	
	v=d.get_value_at(x,y,z)
	if v>threshold : insort(l,(v,x,y,z))
	
if __name__ == "__main__":
    main()
