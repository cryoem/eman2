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
	usage = """%prog [options] <skeleton map> <dejavu file> <output>
	
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
	skeleton.set_attr("apix_x",1.06)
	skeleton.set_attr("apix_y",1.06)
	skeleton.set_attr("apix_z",1.06)
	skeleton.set_attr("MRC.nxstart",1)
	skeleton.set_attr("MRC.nystart",1)
	skeleton.set_attr("MRC.nzstart",1)

	dejavupoints=getPoints(args[1],options.apix)
	mappoints=[getNearest(i[:3],3,skeleton) for i in dejavupoints]
	mappoints+=[getNearest(i[3:],3,skeleton) for i in dejavupoints]
	
# 	for s in mappoints:
# 		print s[0],s[1],s[2],skeleton.get_value_at(*s)
	

	for i in range(len(mappoints)/2):
		erasePairs(skeleton,mappoints[i+len(mappoints)/2],((mappoints[i],()),),i+2)
	
	skeleton.process_inplace("eman1.threshold.binary",{"value":0.5})
#	skeleton.write_image("skel.noh.mrc")
# 	for s in mappoints:
# 		print s[0],s[1],s[2],skeleton.get_value_at(*s)
		
	pairlist=[]
	for i,j in enumerate(mappoints):
		findPath(skeleton,mappoints,[(j[0],j[1],j[2],0)],1,pairlist,i+2)

#	print pairlist
	
#	for i in pairlist:
#		print i[0],i[1],mappoints[i[0]],mappoints[i[1]],vecdist(mappoints[i[0]],mappoints[i[1]]),i[2]
	
	print "%d paths detected"%len(pairlist)
#	skeleton.write_image("zz.mrc")
	
	pts=len(dejavupoints)
	pairlist=[((i[0]%pts,i[0]/pts),(i[1]%pts,i[1]/pts),i[2]) for i in pairlist]
	
	out=file(args[2],"w")
	for i in pairlist: 
		if i[0][0]!=i[1][0] : out.write("%d %d\t%d %d\t%1.3f\n"%(i[0][1],i[0][0],i[1][1],i[1][0],i[2]*options.apix))
	out.close()
	
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
						bestcoord=(coord[0]+dx, coord[1]+dy, coord[2]+dz)
	#print coord, bestcoord
	return(bestcoord)

def zapPairs(skeleton,mp):
	"""Erases a path connecting the pairs of Helix endpoints, so no connectivities
	following paths through helices will be considered. 'mp' is the mappoints array.
	Doesn't work very well in most cases..."""
	n=len(mp)/2
	for i in range(n):
		l=sqrt((mp[i][0]-mp[i+n][0])**2+(mp[i][1]-mp[i+n][1])**2+(mp[i][2]-mp[i+n][2])**2)
		l=floor(l)
		for f in [x/l for x in range(3,int(l)-2)]:
			x=f*mp[i][0]+(1.0-f)*mp[i+n][0]
			y=f*mp[i][1]+(1.0-f)*mp[i+n][1]
			z=f*mp[i][2]+(1.0-f)*mp[i+n][2]
			setbox(int(x),int(y),int(z),skeleton,0)

def removeHelix(startpoint, endpoint, skeleton, maxdistance):
	searchrange=range(-2,3)
	STOPFLAG=0
	if maxdistance<0: 
		maxdistance=sqrt((startpoint[0]-endpoint[0])**2+(startpoint[1]-endpoint[1])**2+(startpoint[2]-endpoint[2])**2)
	
	newpoints=[startpoint,endpoint]
	nmd=maxdistance
	for dx in searchrange:
		for dy in searchrange:
			for dz in searchrange:
				sp=[startpoint[0]+dx, startpoint[1]+dy, startpoint[2]+dz]
				ep=[endpoint[0]+dx, endpoint[1]+dy, endpoint[2]+dz]

				toep=sqrt((sp[0]-endpoint[0])**2+(sp[1]-endpoint[1])**2+(sp[2]-endpoint[2])**2)
				if toep<maxdistance and skeleton.get_value_at(sp[0],sp[1],sp[2])!=0:
					newpoints[0]=sp
					nmd=min(nmd,toep)
					skeleton.set_value_at(sp[0],sp[1],sp[2],0)
					STOPFLAG=1

				tosp=sqrt((startpoint[0]-ep[0])**2+(startpoint[1]-ep[1])**2+(startpoint[2]-ep[2])**2)
				if tosp<maxdistance and skeleton.get_value_at(ep[0],ep[1],ep[2])!=0:
					newpoints[1]=ep
					nmd=min(nmd,tosp)
					skeleton.set_value_at(ep[0],ep[1],ep[2],0)
					STOPFLAG=1
	maxdistance=nmd
	if STOPFLAG==1:
		removeHelix(newpoints[0],newpoints[1], skeleton, maxdistance)


def erasePairs(skeleton,target,seeds,n):
	"""Iteratively erases the shortest path connecting two points (inital seed and target)
	seed must be passed as ((x,y,z),()),)"""
	newseeds=[]
	
#	print n,len(seeds),seeds[0][0]
	for s in seeds:
		if s[0]==target :
#			print "trace ",len(s[1])
			for i in s[1][3:-3]:
#				setbox(i[0],i[1],i[2],skeleton,0,1)
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
						skeleton.set_value_at(x,y,z,n)
	if len(newseeds):
		erasePairs(skeleton,target,newseeds,n)


def findPath(skeleton,mappoints,seeds,it,pairs,n):
	# Iterate over all current trace edge points
	newseeds=[]
	mp2=len(mappoints)/2
	for s in seeds:
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
								pairs.append((n-2,i,s[3]+l))
								setbox(x,y,z,skeleton,n)
								continue
							
	for s in seeds:
		# 2nd pass to find new seeds away from new endpoints
		for z in range(s[2]-1,s[2]+2):
			for y in range(s[1]-1,s[1]+2):
				for x in range(s[0]-1,s[0]+2):
					
					if skeleton.get_value_at(x,y,z)>0 and skeleton.get_value_at(x,y,z)!=n:
						l=sqrt((z-s[2])**2+(y-s[1])**2+(x-s[0])**2)
						newseeds.append((x,y,z,l+s[3]))
						skeleton.set_value_at(x,y,z,n)
	if len(newseeds):
		if n==2 :
#			if it<20 : skeleton.write_image("zzzz.%02d.mrc"%it)
			print n,it,len(newseeds)
		findPath(skeleton,mappoints,newseeds,it+1,pairs,n)

def setbox(x,y,z,img,n,rng=2):
	for xx in range(x-rng,x+rng+1):
		for yy in range(y-rng,y+rng+1):
			for zz in range(z-rng,z+rng+1):
				if img.get_value_at(xx,yy,zz) : img.set_value_at(xx,yy,zz,n)
				
if __name__ == "__main__":
	main()
