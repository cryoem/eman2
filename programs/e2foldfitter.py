#!/bin/env python

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
from bisect import insort


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: %prog [options] target.mrc probe.mrc
	
Locates the best 'docking' locations for a small probe in a large target map."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--apix", "-A", type="float", help="A/voxel", default=1.0)
	parser.add_option("--res", "-R", type="float", help="Resolution in A, equivalent to Gaussian lowpass with 1/e width at 1/res",default=2.8)
	parser.add_option("--box", "-B", type="string", help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	parser.add_option("--het", action="store_true", help="Include HET atoms in the map", default=False)
	parser.add_option("--chains",type="string",help="String list of chain identifiers to include, eg 'ABEFG'")
	parser.add_option("--quiet",action="store_true",default=False,help="Verbose is the default")
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	try: chains=options.chains
	except: chains=None
	
	try : infile=open(args[0],"r")
	except : parser.error("Cannot open input file")
	
	# read the target and probe
	target=EMData()
	target.read_image(args[0])
	
	probe=EMData()
	probe.read_image(args[1])
	
	tdim=(target.get_xsize(),target.get_ysize(),target.get_zsize())
	pdim=(probe.get_xsize(),probe.get_ysize(),probe.get_zsize())
	
	if (pdim[0]>tdim[0] or pdim[1]>tdim[1] or pdim[2]>tdim[2]) :
		print "Probe must fit within target"
		exit(1)
	
	# shrink both by some factor which keeps the smallest axis of the probe at least 10 pixels
	# we'll have to reread the files if we want to recover the unscaled images
	sfac=int(floor(min(pdim)/10.0))
	print "Shrink by %d"%sfac
	target.mean_shrink(sfac)
	probe.mean_shrink(sfac)
	tdim2=(target.get_xsize(),target.get_ysize(),target.get_zsize())
	pdim2=(probe.get_xsize(),probe.get_ysize(),probe.get_zsize())
	probe.filter("NormalizeEdgeMean")
	probeclip=probe.get_clip(Region((pdim2[0]-tdim2[0])/2,(pdim2[1]-tdim2[1])/2,(pdim2[2]-tdim2[2])/2,tdim2[0],tdim2[1],tdim2[2]))
	
	roughang=[(0,0),(45,0),(45,90),(45,180),(45,270),(90,0),(90,60),(90,120),(90,180),(90,240),(90,300),(135,0),(135,90),(135,180),(135,270),(180,0)]

#	Log.logger().set_level(Log.LogLevel.DEBUG_LOG)
	
	edge=max(pdim2)/2		# technically this should be max(pdim), but generally there is some padding in the probe model, and this is relatively harmless
	print "edge ",edge
	best=[]
	for a1,a2 in roughang:
		for a3 in range(0,360,45):
			prr=probeclip.copy()
			prr.rotate(a1,a2,a3)
#			prr.write_image('prr.%0d%0d%0d.mrc'%(a1,a2,a3))
			
			ccf=target.calc_ccf(prr,1)
			mean=float(ccf.get_attr("mean"))
			sig=float(ccf.get_attr("sigma"))
			ccf.filter("ZeroEdgePlane",{"x0":edge,"x1":edge,"y0":edge,"y1":edge,"z0":edge,"z1":edge})
			ccf.filter("PeakOnly",{"npeaks":0})		# only look at peak values in the CCF map
			
#			ccf.write_image('ccf.%0d%0d%0d.mrc'%(a1,a2,a3))
			vec=ccf.calc_highest_locations(mean+sig+.0000001)
			for v in vec: best.append([v.value,a1,a2,a3,v.x,v.y,v.z,0])
			
#			print a1,a2,a3,mean+sig,float(ccf.get_attr("max")),len(vec)
	
	best.sort()		# this is a list of all reasonable candidate locations
	
	print len(best)
	print best[0]
	print best[-1]
	
	for i in best:
		for j in best:
			if (i[4]-j[4])**2+(i[5]-j[5])**2+(i[6]-j[6])**2>8.8 : continue
			if j[0]>i[0] : i[7]=1
	
	best2=[]
	for i in best:
		if not i[7]: best2.append(i)

	print len(best2)
	for i in best2: print i
#	print best2[0]
#	print best2[-1]
		
if __name__ == "__main__":
    main()
