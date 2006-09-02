#!/usr/bin/env python
# e2skelpath.py  09/01/2006  Steven Ludtke, Matt Baker

from EMAN2 import *
from optparse import OptionParser
from math import *
import time
import os,re
import sys
from pprint import pprint

pl=()

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <skeleton map> <dejavu file>
	
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
#	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: circle, ref, grid",default=[])
#	parser.add_option("--threshold","-T",type="float",help="Threshold for keeping particles. 0-4, 0 excludes all, 4 keeps all.",default=2.0)
# 	parser.add_option("--maxbad","-M",type="int",help="Maximumum number of unassigned helices",default=2)
# 	parser.add_option("--minhelix","-H",type="int",help="Minimum residues in a helix",default=6)
 	parser.add_option("--apix","-P",type="float",help="A/Pixel",default=1.0)
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	logid=E2init(sys.argv)
		
	skeleton=EMData()
	skeleton.read_image(args[0])

	dejavupoints=getPoints(args[1],options.apix)
	mappoints=[getNearest(i[:3],3,skeleton) for i in dejavupoints]
	mappoints+=[getNearest(i[3:],3,skeleton) for i in dejavupoints]
	
# 	for s in mappoints:
# 		print s[0],s[1],s[2],skeleton.get_value_at(*s)
	
	pairlist=[]
	for i,j in enumerate(mappoints):
		findPath(skeleton,mappoints,[j],1,pairlist,i+2)

	print pairlist
	
	for i in pairlist:
		print i[0],i[1],mappoints[i[0]],mappoints[i[1]],vecdist(mappoints[i[0]],mappoints[i[1]]),i[2]
	
def vecdist(a,b):
	return sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)	

def getPoints(dejavufile, apix):
	pattern=re.compile(r"ALPHA\s'(?P<chain>[\w]+)'\s'(?P<startres>[\w]+)'\s'(?P<stopres>[\w]+)'\s(?P<reslength>[\d]+)\s(?P<x1>[\d.,-]+)\s(?P<y1>[\d.,-]+)\s(?P<z1>[\d.,-]+)\s(?P<x2>[\d.,-]+)\s(?P<y2>[\d.,-]+)\s(?P<z2>[\d.,-]+)")
	endpoints=[]
	#print dejavufile, apix
	for line in file(dejavufile,"r").readlines():
		result = pattern.search(line)
		if result:
			coord=[int(float(result.group('x1'))/apix), int(float(result.group('y1'))/apix), int(float(result.group('z1'))/apix), int(float(result.group('x2'))/apix), int(float(result.group('y2'))/apix), int(float(result.group('z2'))/apix)]
			endpoints.append(coord)
	return(endpoints)
	

def getNearest(coord, searchrange, skeleton):
	maxdistance=10000
	bestcoord=coord
	searchrange=range(-searchrange,searchrange+1,1)
	for dx in searchrange:
		for dy in searchrange:
			for dz in searchrange:
				if skeleton.get_value_at(coord[0]+dx, coord[1]+dy, coord[2]+dz)>=0.1:
					distance=sqrt(dx**2+dy**2+dz**2)
					if distance < maxdistance:
						maxdistance=distance
						bestcoord=[coord[0]+dx, coord[1]+dy, coord[2]+dz]
	#print coord, bestcoord
	return(bestcoord)

def findPath(skeleton,mappoints,seeds,it,pairs,n):
	# Iterate over all current trace edge points
	newseeds=[]
	for s in seeds:
		# search nearest neighbors to continue trace

		skeleton.set_value_at(s[0],s[1],s[2],n)
		for z in range(s[2]-1,s[2]+2):
			for y in range(s[1]-1,s[1]+2):
				for x in range(s[0]-1,s[0]+2):
					
#					print n-2,it,x,y,z,skeleton.get_value_at(x,y,z)
					if skeleton.get_value_at(x,y,z)>0 and skeleton.get_value_at(x,y,z)!=n:
						for i,j in enumerate(mappoints):
							if i==n: continue
							if [x,y,z]==j :
								pairs.append((n-2,i,it))
								setbox(x,y,z,skeleton,n)
								continue
						newseeds.append((x,y,z))
						skeleton.set_value_at(x,y,z,n)
	if len(newseeds):
		findPath(skeleton,mappoints,newseeds,it+1,pairs,n)

def setbox(x,y,z,img,n):
	for xx in range(x-1,x+2):
		for yy in range(y-1,y+2):
			for zz in range(z-1,z+2):
				img.set_value_at(xx,yy,zz,n)
				
if __name__ == "__main__":
	main()
