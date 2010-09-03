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

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
from Simplex import Simplex
from bisect import insort

cmp_probe=None
cmp_target=None
tdim=None
pdim=None
c2alt=0

def display(img):
	img.write_image("tmploc.mrc")
	os.system("v2 tmploc.mrc")

def compare(vec):
	"""Given an (alt,az,phi,x,y,z) vector, calculate the similarity
	of the probe to the map"""
	global cmp_probe,cmp_target
	
#	print vec,pdim
#	print "\n%6.3f %6.3f %6.3f    %5.1f %5.1f %5.1f"%(vec[0],vec[1],vec[2],vec[3],vec[4],vec[5])
	a=cmp_target.get_rotated_clip(Transform3D([vec[3]+tdim[0]/2,vec[4]+tdim[1]/2,vec[5]+tdim[2]/2],*vec[0:3]),pdim,1.0)
#	a.write_image("clip.mrc")
#	os.system("v2 clip.mrc")
	return -cmp_probe.cmp("dot",a,{})

def compare2(vec):
	"""Given an (alt,az,phi,x,y,z) vector, calculate the similarity
	of the probe to the map"""
	global cmp_probe,cmp_target,c2alt
	
#	print vec,pdim
#	print "\n%6.3f %6.3f %6.3f    %5.1f %5.1f %5.1f"%(vec[0],vec[1],vec[2],vec[3],vec[4],vec[5])
#	a=cmp_target.get_rotated_clip((tdim[0]/2,tdim[1]/2,vec[1]+tdim[2]/2),Rotation(c2alt,0,vec[0],Rotation.EulerType.EMAN),pdim,1.0)
	a=cmp_target.get_rotated_clip(Transform3D([tdim[0]/2,tdim[1]/2,vec[1]+tdim[2]/2],c2alt,0,vec[0]),pdim,1.0)
#	a.write_image("clip.mrc")
#	os.system("v2 clip.mrc")
	return -cmp_probe.cmp("dot",a,{})
	
def main():
	global tdim,pdim
	global cmp_probe,cmp_target
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] fixed.mrc input.mrc output.mrc

	This program is not currently considered stable.
	
Locates the best 'docking' locations for a small probe in a large target map."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--csym", action="store_true", help="Restrict axes for c/d symmetric objects", default=False)
	
	(options, args) = parser.parse_args()
	if len(args)<3 : parser.error("Input and output files required")
	try: chains=options.chains
	except: chains=None

	print "WARNING: This program is not currently considered stable"
		
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
#	sfac=int(floor(min(pdim)/10.0))
	sfac=int(floor(min(pdim)/16.0))
	while min(pdim)%(sfac)!=0 or (min(pdim)/sfac)%2!=0 : sfac-=1 
	print "Shrink by %d"%sfac
	target.process_inplace("math.meanshrink",{"n":sfac})
	probe.process_inplace("math.meanshrink",{"n":sfac})
	tdim2=(target.get_xsize(),target.get_ysize(),target.get_zsize())
	pdim2=(probe.get_xsize(),probe.get_ysize(),probe.get_zsize())
#	print (pdim2[0]-tdim2[0])/2,(pdim2[1]-tdim2[1])/2,(pdim2[2]-tdim2[2])/2,tdim2[0],tdim2[1],tdim2[2]
	probe.process_inplace("normalize.edgemean")
	probeclip=probe.get_clip(Region((pdim2[0]-tdim2[0])/2,(pdim2[1]-tdim2[1])/2,(pdim2[2]-tdim2[2])/2,tdim2[0],tdim2[1],tdim2[2]))
	
	if options.csym: roughang=[(0,0),(180,0)]
	else : roughang=[(0,0),(45,0),(45,90),(45,180),(45,270),(90,0),(90,60),(90,120),(90,180),(90,240),(90,300),(135,0),(135,90),(135,180),(135,270),(180,0)]

#	Log.logger().set_level(Log.LogLevel.DEBUG_LOG)
	
#	probeclip.write_image("m0.mrc")
	best=[]
#	sm=probeclip.copy_head()
#	sm.to_zero()
	
	sm=EMData()
	sm.set_size(probeclip.get_xsize(),probeclip.get_ysize(),probeclip.get_zsize())
	sm.to_one()	
	for a1,a2 in roughang:
		for a3 in range(0,360,45):
			prr=probeclip.copy()
			prr.rotate(a1,a2,a3)
#			print a1,a2,a3
#			display(prr)
#			prr.write_image('prr.%0d%0d%0d.mrc'%(a1,a2,a3))
			
#			target.write_image("m1.mrc")
#			prr.write_image("m2.mrc")
			ccf=target.calc_ccf(prr,fp_flag.CIRCULANT)
			mean=float(ccf.get_attr("mean"))
			sig=float(ccf.get_attr("sigma"))
#			ccf.write_image('ccf1.%0d%0d%0d.mrc'%(a1,a2,a3))
			ccf.process_inplace("mask.zeroedge3d",{"x0":pdim2[0]/4,"x1":pdim2[0]/4,"y0":pdim2[1]/4,"y1":pdim2[1]/4,"z0":pdim2[2]/4,"z1":pdim2[2]/4})
			ccf.process_inplace("mask.onlypeaks",{"npeaks":1})		# only look at peak values in the CCF map
			ccf.write_image('ccf2.%0d%0d%0d.mrc'%(a1,a2,a3))
			print ccf.get_attr("sigma"),
			sm.add(ccf)
			sm.write_image("sm.mrc")
			print sm.get_attr("maximum")

#			print mean,sig
#			ccf.write_image('ccf.%0d%0d%0d.mrc'%(a1,a2,a3))
			vec=ccf.calc_highest_locations(mean+sig+.0000001)
			for v in vec: best.append([v.value,a1,a2,a3,v.x-tdim2[0]/2,v.y-tdim2[1]/2,v.z-tdim2[2]/2,0])
#			print len(vec)
#			for v in vec: print v.value,v.x,v.y,v.z
#			print a1,a2,a3,mean+sig,float(ccf.get_attr("max")),len(vec)
	
	best.sort()		# this is a list of all reasonable candidate locations
	best.reverse()
#	for i in best: print i
	
	print len(best)," candidate locations"
	
	# this is designed to eliminate angular redundancies in peak location
#	print best[0]
#	print best[-1]
	
#	for i in best:
#		for j in best:
#			if (i[4]-j[4])**2+(i[5]-j[5])**2+(i[6]-j[6])**2>8.8 : continue
#			if j[0]>i[0] : i[7]=1
#	
#	best2=[]
#	for i in best:
#		if not i[7]: best2.append(i)

	# now we find peaks in the sm of all CCF calculations, and find the best angle associated with each peak
#	sm.write_image("sm.mrc")
	sm.process_inplace("mask.onlypeaks",{"npeaks":1})
	vec=sm.calc_highest_locations(mean+sig+.0000001)
	best2=[]
	for v in vec:
#		print "%5.1f  %5.1f  %5.1f"%(v.x*sfac-tdim[0]/2,v.y*sfac-tdim[1]/2,v.z*sfac-tdim[2]/2)
		for i in best:
			if i[4]+tdim2[0]/2==v.x and i[5]+tdim2[1]/2==v.y and i[6]+tdim2[2]/2==v.z :
				best2.append([i[0],i[1],i[2],i[3],i[4]*sfac,i[5]*sfac,i[6]*sfac,i[7]])
				break

	best2.sort()
	best2.reverse()
	print len(best2)," final candidtates"
	for i in best2: print i

	"""print "Preoptimize:"
	cmp_target=target
	cmp_probe=probe
	if options.csym :
		for j in range(len(best2)):
			c2alt=best2[j][1]
			sm=Simplex(compare2,[best2[j][3],best2[j][6]],[1,2.])
			bt=sm.minimize()
			b=bt[0]
			print "\n",j,"\t%5.2f,%5.2f  %5.1f"%(c2alt/degrad,b[0]/degrad,b[1])
			best2[j][3]=b[0]
			best2[j][6]=b[1]
			
	else :
		for j in range(len(best2)):
			sm=Simplex(compare,best2[j][1:7],[1,1,1,2.,2.,2.])
			bt=sm.minimize()
			b=bt[0]
			print "\n",j,"\t(%5.2f  %5.2f  %5.2f    %5.1f  %5.1f  %5.1f"%(b[0]/degrad,b[1]/degrad,b[2]/degrad,b[3],b[4],b[5])
			best2[j][1:7]=b"""
			

		
	# reread the original images
	target.read_image(args[0])
	probe.read_image(args[1])
	cmp_target=target
	cmp_probe=probe
	
#	for i in best2:
#		c=probe.get_clip(Region((pdim[0]-tdim[0])/2,(pdim[1]-tdim[1])/2,(pdim[2]-tdim[2])/2,tdim[0],tdim[1],tdim[2]))
#		c.rotate_translate(*i[1:7])
#		c.write_image("z.%02d.mrc"%best2.index(i))
	
	print "Optimize:"
	if options.csym!=False :
		for j in range(len(best2)):
			c2alt=best2[j][1]
			sm=Simplex(compare2,[best2[j][3],best2[j][6]],[1,2.])
			bt=sm.minimize()
			b=bt[0]
			print "\n",j,"\t%5.2f,%5.2f  %5.1f"%(c2alt,b[0],b[1])
#			a=cmp_target.get_rotated_clip((b[3]+tdim[0]/2,b[4]+tdim[1]/2,b[5]+tdim[2]/2),*b[0:3],EULER_EMAN,pdim,1.0)
#			a.write_image("clip.%02d.mrc"%j)
			pc=probe.get_clip(Region((pdim[0]-tdim[0])/2,(pdim[1]-tdim[1])/2,(pdim[2]-tdim[2])/2,tdim[0],tdim[1],tdim[2]))
	#		pc.write_image("finala.mrc")
			pc.rotate_translate(c2alt,0,-b[0],0,0,b[1])
			pc.write_image("final.%02d.mrc"%j)
			
		return
			
	for j in range(len(best2)):
		sm=Simplex(compare,best2[j][1:7],[1,1,1,2.,2.,2.])
		bt=sm.minimize()
		b=bt[0]
		print "\n",j,"\t(%5.2f  %5.2f  %5.2f    %5.1f  %5.1f  %5.1f"%(b[0],b[1],b[2],b[3],b[4],b[5])
		a=cmp_target.get_rotated_clip(Transform3D([b[3]+tdim[0]/2,b[4]+tdim[1]/2,b[5]+tdim[2]/2],*b[0:3]),pdim,1.0)
		a.write_image("clip.%02d.mrc"%j)
		pc=probe.get_clip(Region((pdim[0]-tdim[0])/2,(pdim[1]-tdim[1])/2,(pdim[2]-tdim[2])/2,tdim[0],tdim[1],tdim[2]))
                pc.rotate(-b[0],-b[2],-b[1])
                pc.rotate_translate(0,0,0,b[3],b[4],b[5])               # FIXME, when rotate_translate wi
		pc.write_image("final.%02d.mrc"%j)

	
#	print compare(best2[0][1:7])
	
#	print best2[0]
#	print best2[-1]
		
if __name__ == "__main__":
    main()
