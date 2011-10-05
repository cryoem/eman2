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

# e2skelpath.py  09/01/2006  Steven Ludtke, Matt Baker

from EMAN2 import *
from math import *
import time
import os,re
import sys
from pprint import pprint

pl=()

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <skeleton map> <dejavu file> <output>

	Experimental program for tracing paths in skeletonized density maps. Contact mbaker@bcm.edu for more details.
	
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

# 	parser.add_argument("--maxbad","-M",type=int,help="Maximumum number of unassigned helices",default=2)
# 	parser.add_argument("--minhelix","-H",type=int,help="Minimum residues in a helix",default=6)
 	parser.add_argument("--apix","-P",type=float,help="A/Pixel",default=1.0)
 	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
 	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	logid=E2init(sys.argv,options.ppid)

	skeleton=EMData()
	skeleton.read_image(args[0])
	
	originX=0.5*skeleton.get_xsize()
	originY=0.5*skeleton.get_ysize()
	originZ=0.5*skeleton.get_zsize()
	dejavupoints=getPoints(args[1],options.apix, originX,originY,originZ)
	mappoints=[getNearest(i[:3],3,skeleton) for i in dejavupoints]
	mappoints+=[getNearest(i[3:],3,skeleton) for i in dejavupoints]	

	for i in range(len(mappoints)/2):
		print "Helix %d:  "%(i), mappoints[i+len(mappoints)/2], mappoints[i]
		erasePairs(skeleton,mappoints[i+len(mappoints)/2],((mappoints[i],()),),i+2)		
		
		"""startpoint=mappoints[i+len(mappoints)/2]
		endpoint=mappoints[i]
		helixdistance=vecdist(startpoint,endpoint)
		midpoint=(startpoint[0]+endpoint[0]/2,startpoint[1]+endpoint[1]/2,startpoint[2]+endpoint[2]/2)	
		nearestmidpoint=getNearest(midpoint,3, skeleton)
		removeHelix(mappoints[i+len(mappoints)/2], mappoints[i], skeleton, nearestmidpoint)
		removeHelix(mappoints[i], mappoints[i+len(mappoints)/2], skeleton, nearestmidpoint)"""
		
	#skeleton.process_inplace("threshold.binary",{"value":0.5})
	skeleton.write_image("skel.noh.mrc")
	
	pairlist=[]
	for i,j in enumerate(mappoints):
		findPath(skeleton,mappoints,[((j[0],j[1],j[2],0),[j[:3]])],1,pairlist,i+2)

#	print pairlist
#	for i in pairlist:
#		print i[0],i[1],mappoints[i[0]],mappoints[i[1]],vecdist(mappoints[i[0]],mappoints[i[1]]),i[2]
	
	print "%d paths detected"%len(pairlist)
#	skeleton.write_image("zz.mrc")
	
	pts=len(dejavupoints)
	pairlist=[((i[0]%pts,i[0]/pts),(i[1]%pts,i[1]/pts),i[2],j) for i,j in pairlist]
	
	out=file(args[2],"w")
	for i in pairlist: 
		if i[0][0]!=i[1][0] : out.write("%d %d\t%d %d\t%1.3f\t%s\n"%(i[0][1],i[0][0],i[1][1],i[1][0],i[2]*options.apix,i[3]))
	out.close()
	
def vecdist(a,b):
	return sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)	

def getPoints(dejavufile, apix,originX,originY,originZ):
	"""get the endpoint of the helices from a dejavu sse list. coordinates are in angstroms; origin is center of map"""
	
	pattern=re.compile(r"ALPHA\s'(?P<chain>[\w]+)'\s'(?P<startres>[\w]+)'\s'(?P<stopres>[\w]+)'\s(?P<reslength>[\d]+)\s(?P<x1>[\d.,-]+)\s(?P<y1>[\d.,-]+)\s(?P<z1>[\d.,-]+)\s(?P<x2>[\d.,-]+)\s(?P<y2>[\d.,-]+)\s(?P<z2>[\d.,-]+)")
	endpoints=[]
	#print "in get points"
	for line in file(dejavufile,"r").readlines():
		result = pattern.search(line)
		if result:
			coord=[int(originX+(float(result.group('x1'))/apix)), int(originY+(float(result.group('y1'))/apix)), int(originZ+(float(result.group('z1'))/apix)), int(originX+(float(result.group('x2'))/apix)), int(originY+(float(result.group('y2'))/apix)), int(originZ+(float(result.group('z2'))/apix))]
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
						bestcoord=(coord[0]+dx, coord[1]+dy, coord[2]+dz)
	return(bestcoord)


"""def removeHelix(startpoint, endpoint, skeleton, maxdistance):
	searchrange=range(-1,2)
	STOPFLAG=0	
	newpoints=[startpoint,endpoint]
	nmd=maxdistance
	for dx in searchrange:
		for dy in searchrange:
			for dz in searchrange:
				sp=[startpoint[0]+dx, startpoint[1]+dy, startpoint[2]+dz]
				#ep=[endpoint[0]+dx, endpoint[1]+dy, endpoint[2]+dz]
				if sp!=startpoint or  sp!=endpoint:
					toep=vecdist(sp,endpoint)
					if toep<maxdistance and skeleton.get_value_at(sp[0],sp[1],sp[2])!=0:
						newpoints[0]=sp
						nmd=min(nmd,toep)
						skeleton.set_value_at(sp[0],sp[1],sp[2],0)
						STOPFLAG=1

				#tosp=sqrt((startpoint[0]-ep[0])**2+(startpoint[1]-ep[1])**2+(startpoint[2]-ep[2])**2)
				#if tosp<maxdistance and skeleton.get_value_at(ep[0],ep[1],ep[2])!=0:
				#	newpoints[1]=ep
				#	nmd=min(nmd,tosp)
				#	skeleton.set_value_at(ep[0],ep[1],ep[2],0)
				#	STOPFLAG=1
	maxdistance=nmd
	if STOPFLAG==1:
		removeHelix(newpoints[0],newpoints[1], skeleton, maxdistance)
		
def removeHelix(startpoint, endpoint, skeleton, nearestmidpoint):
	searchrange=range(-1,2)
	STOPFLAG=0
	maxdistance=max(vecdist(startpoint,nearestmidpoint),vecdist(endpoint,nearestmidpoint))
	for dx in searchrange:
		for dy in searchrange:
			for dz in searchrange:
				sp=[nearestmidpoint[0]+dx, nearestmidpoint[1]+dy, nearestmidpoint[2]+dz]
				if sp!=startpoint or  sp!=endpoint:
					print "yes"
					toend=vecdist(sp,endpoint)
					if toend<maxdistance and skeleton.get_value_at(sp[0],sp[1],sp[2])!=0:
							print "dude"
							nearestmidpoint=sp
							maxdistance=toend
							STOPFLAG=1
				if STOPFLAG==1:
					skeleton.set_value_at(nearestmidpoint[0],nearestmidpoint[1],nearestmidpoint[2],0)
	if STOPFLAG==1:
		removeHelix(startpoint,endpoint, skeleton, nearestmidpoint)"""					



def erasePairs(skeleton,target,seeds,n):
	"""Iteratively erases the shortest path connecting two points (inital seed and target)
	seed must be passed as ((x,y,z),()),)"""
	newseeds=[]

	print n,len(seeds),seeds[0][0], skeleton.get_value_at(seeds[0][0][0],seeds[0][0][1],seeds[0][0][2])
	for s in seeds:
		if s[0]==target :
			print "trace ",len(s[1])
			for i in s[1][5:-5]:
#				setbox(i[0],i[1],i[2],skeleton,0,1)
				print "I am setting this point to zero: %d, %d, %d "%(i[0],i[1],i[2])
				skeleton.set_value_at(i[0],i[1],i[2],0)
			return

	for ss in seeds:
		s=ss[0]
		# 2nd pass to find new seeds away from new endpoints
		for z in range(s[2]-1,s[2]+2):
			for y in range(s[1]-1,s[1]+2):
				for x in range(s[0]-1,s[0]+2):
					
					if skeleton.get_value_at(x,y,z)>0 and skeleton.get_value_at(x,y,z)!=n:
						newseeds.append(((x,y,z),ss[1]+(ss[0],)))
						print "Setting %d,%d,%d to %f"%(x,y,z,n)
						skeleton.set_value_at(x,y,z,n)
	if len(newseeds):
		erasePairs(skeleton,target,newseeds,n)


def findPath(skeleton,mappoints,seeds,it,pairs,n):
	# Iterate over all current trace edge points
	newseeds=[]
	mp2=len(mappoints)/2
	for s,p in seeds:
		# search nearest neighbors to continue trace

		skeleton.set_value_at(s[0],s[1],s[2],n)
		for z in range(s[2]-1,s[2]+2):
			for y in range(s[1]-1,s[1]+2):
				for x in range(s[0]-1,s[0]+2):
					
#					print n-2,it,x,y,z,skeleton.get_value_at(x,y,z)
					if skeleton.get_value_at(x,y,z)>0 and skeleton.get_value_at(x,y,z)!=n:
						l=sqrt((z-s[2])**2+(y-s[1])**2+(x-s[0])**2)
						for i,j in enumerate(mappoints):
							if i==n: continue
							if (x,y,z)==j :
								pairs.append(((n-2,i,s[3]+l),p+[(x,y,z)]))
								setbox(x,y,z,skeleton,n)
								continue
							
	for s,p in seeds:
		# 2nd pass to find new seeds away from new endpoints
		for z in range(s[2]-1,s[2]+2):
			for y in range(s[1]-1,s[1]+2):
				for x in range(s[0]-1,s[0]+2):
					
					if skeleton.get_value_at(x,y,z)>0 and skeleton.get_value_at(x,y,z)!=n:
						l=sqrt((z-s[2])**2+(y-s[1])**2+(x-s[0])**2)
						newseeds.append(((x,y,z,l+s[3]),p+[(x,y,z)]))
						skeleton.set_value_at(x,y,z,n)
	if len(newseeds):
#		if n==2 :
#			if it<20 : skeleton.write_image("zzzz.%02d.mrc"%it)
#			print n,it,len(newseeds)
		findPath(skeleton,mappoints,newseeds,it+1,pairs,n)

def setbox(x,y,z,img,n,rng=2):
	for xx in range(x-rng,x+rng+1):
		for yy in range(y-rng,y+rng+1):
			for zz in range(z-rng,z+rng+1):
				if img.get_value_at(xx,yy,zz) : img.set_value_at(xx,yy,zz,n)
				
if __name__ == "__main__":
	main()
