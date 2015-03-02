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
from math import *
import os
import sys
from Simplex import Simplex
from bisect import insort

cmp_probe=None
cmp_target=None
tdim=None
pdim=None
tdim2=None
pdim2=None
sfac=None
degrad=pi/180.0
ncmp=0

def compare(vec,empty):	 # I add an empty argument just because the Simplex.py wants two input... -Muyuan
	"""Given an (az,alt,phi,x,y,z) vector, calculate the similarity
	of the probe to the map"""
	global cmp_probe,cmp_target,ncmp
	
#	print vec,pdim
#	print "\n%6.3f %6.3f %6.3f    %5.1f %5.1f %5.1f"%(vec[0],vec[1],vec[2],vec[3],vec[4],vec[5])
	t = Transform()
	t.set_pre_trans((vec[3]+tdim[0]/2,vec[4]+tdim[1]/2,vec[5]+tdim[2]/2))
	#t.set_trans((vec[3]+tdim[0]/2,vec[4]+tdim[1]/2,vec[5]+tdim[2]/2))
	t.set_rotation({'type':'eman', 'az':vec[0], 'alt':vec[1], 'phi':vec[2]})
	#t.set_trans((0,0,0))
	a=cmp_target.get_rotated_clip(t,pdim,1.0)
#	a=cmp_target.get_rotated_clip(Transform3D((vec[3]+tdim[0]/2,vec[4]+tdim[1]/2,vec[5]+tdim[2]/2),vec[0],vec[1],vec[2],(0,0,0)),pdim,1.0)
	
#	os.system("v2 clip.mrc")
	ncmp+=1

	#a.write_image("clip.mrc")
	#print cmp_probe.cmp("dot",a,{}),vec,a["mean"]

	return cmp_probe.cmp("dot",a,{"normalize":1})
	
def compares(vec,empty):
	"""Given an (az,alt,phi,x,y,z) vector, calculate the similarity
	of the probe to the map"""
	global cmp_probe,cmp_target,sfac
	
	t = Transform()
	t.set_pre_trans((vec[3]/float(sfac)+tdim2[0]/2,vec[4]/float(sfac)+tdim2[1]/2,vec[5]/float(sfac)+tdim2[2]/2))
	#t.set_trans((vec[3]/float(sfac)+tdim2[0]/2,vec[4]/float(sfac)+tdim2[1]/2,vec[5]/float(sfac)+tdim2[2]/2))
	t.set_rotation({'type':'eman', 'az':vec[0], 'alt':vec[1], 'phi':vec[2]})
	#t.set_trans((0,0,0))
	a=cmp_target.get_rotated_clip(t,pdim2,1.0)
	#print vec,pdim2
	#print a["mean"]
#	a=cmp_target.get_rotated_clip(Transform3D((vec[3]/float(sfac)+tdim2[0]/2,vec[4]/float(sfac)+tdim2[1]/2,vec[5]/float(sfac)+tdim2[2]/2),vec[0],vec[1],vec[2],(0,0,0)),pdim2,1.0)
	print 100*cmp_probe.cmp("dot",a,{}),vec,a["mean"]
	#cmp_probe.write_image("p.mrc")
	#a.write_image("a.mrc")
	#exit()
	return 100*cmp_probe.cmp("dot",a,{})
	#return -a["mean"]

def pdb_transform(infile,outfile,trans,center):
	pdb=open(infile,"r")
	out=open(outfile,"w")
	for line in pdb:
		if (line[:4]=='ATOM' or line[:6]=='HETATM') :
			x=float(line[30:38])-center[0]
			y=float(line[38:46])-center[1]
			z=float(line[46:54])-center[2]
			v2=trans.transform(x,y,z)+center
			line=list(line)
			line[30:38]= '{:>8}'.format('{:.3f}'.format(v2[0]))
			line[38:46]= '{:>8}'.format('{:.3f}'.format(v2[1]))
			line[46:54]= '{:>8}'.format('{:.3f}'.format(v2[2]))
			line=''.join(line)
		out.write(line)
	pdb.close()
	out.close()
			

def main():
	global tdim,pdim,tdim2,pdim2,sfac
	global cmp_probe,cmp_target
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] target.mrc probe.mrc
	
Locates the best 'docking' locations for a small probe in a large target map. Note that the probe
should be in a box barely large enough for it. The target may be arbitrarily padded. For best speed
both box sizes should be multiples of 8."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--shrink", "-S", type=int, help="shrink factor for initial search, default=auto", default=0)
	parser.add_argument("--num", "-N", type=int, help="Number of initial alternative positions, default=5", default=5)
	parser.add_argument("--epsilon","-E", type=float,help="final target accuracy, default=.01",default=.01)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	print "WARNING: This program is currently considered experimental. Contact sludtke@bcm.edu before using it for any serious project"
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	try: chains=options.chains
	except: chains=None
	logid=E2init(sys.argv,options.ppid)
	
	try : infile=open(args[0],"r")
	except : parser.error("Cannot open input file")
	

	
	# read the target and probe
	target=EMData()
	target.read_image(args[0])
	
	
	apix=target["apix_x"]
	probe=EMData()
	probefilename=args[1]
	# support pdb format
	if args[1].endswith(".pdb"):
		print "e2pdb2mrc.py {s} probe.mrc -R 10 --het --apix={a}>tmp.txt".format(s=args[1],a=apix)
		os.system("e2pdb2mrc.py {s} probe.mrc -R 10 --het --apix={a}>tmp.txt".format(s=args[1],a=apix))
		tmp=open('tmp.txt')
		lines=tmp.readlines()
		tmp.close()
		cent=[]
		for l in range(len(lines)):
			if lines[l].startswith("Bounding box"):
				break
		for q in range(3):
			#print lines[q+l][17:-1].split('-')
			t=[float(i) for i in (lines[q+l][17:-1].split(' - '))]
			#print t
			cent.append((t[0]+t[1])/2)
		
		probefilename="probe.mrc"
	else:
		probefilename=args[1]
	probe.read_image(probefilename)
	
		
	tdim=(target.get_xsize(),target.get_ysize(),target.get_zsize())
	pdim=(probe.get_xsize(),probe.get_ysize(),probe.get_zsize())
	
	if (pdim[0]>tdim[0] or pdim[1]>tdim[1] or pdim[2]>tdim[2]):
		print "Probe must fit within target"
		exit(1)
	target.process_inplace("normalize.unitsum")
	target.mult(10000)
	
	# shrink both by some factor which keeps the smallest axis of the probe at least 10 pixels
	# we'll have to reread the files if we want to recover the unscaled images
#	sfac=int(floor(min(pdim)/10.0))
	if options.shrink>0 : sfac=options.shrink
	else : sfac=int(floor(min(pdim)/12.0))
	print "Shrink by %d"%sfac
	target.process_inplace("math.meanshrink",{"n":sfac})
	probe.process_inplace("math.meanshrink",{"n":sfac})
	tdim2=(target.get_xsize(),target.get_ysize(),target.get_zsize())
	pdim2=(probe.get_xsize(),probe.get_ysize(),probe.get_zsize())
#	print (pdim2[0]-tdim2[0])/2,(pdim2[1]-tdim2[1])/2,(pdim2[2]-tdim2[2])/2,tdim2[0],tdim2[1],tdim2[2]
	probe.process_inplace("normalize.edgemean")

	probeclip=probe.get_clip(Region((pdim2[0]-tdim2[0])/2,(pdim2[1]-tdim2[1])/2,(pdim2[2]-tdim2[2])/2,tdim2[0],tdim2[1],tdim2[2]))
	#roughang=[(0,0)]

	#roughang=[(0,0),(45,0),(45,90),(45,180),(45,270),(90,0),(90,60),(90,120),(90,180),(90,240),(90,300),(135,0),(135,90),(135,180),(135,270),(180,0)]
	roughang=[(0,0),(30,0),(30,90),(30,180),(30,270),(60,0),(60,45),(60,90),(60,135),(60,180),(60,225),(60,270),(60,315),
	(90,0),(90,30),(90,60),(90,90),(90,120),(90,150),(90,180),(90,210),(90,240),(90,270),(90,300),(90,330),
	(180,0),(150,0),(150,90),(150,180),(150,270),(120,0),(120,45),(120,90),(120,135),(120,180),(120,225),(120,270),(120,315)]

#	Log.logger().set_level(Log.LogLevel.DEBUG_LOG)
	
	print "Searching for candidate locations in reduced map"
	edge=max(pdim2)/2		# technically this should be max(pdim), but generally there is some padding in the probe model, and this is relatively harmless
	print "edge ",edge
	best=[]
	sum=probeclip.copy_head()
	sum.to_zero()
	for a1,a2 in roughang:
		for a3 in range(0,360,30):
			prr=probeclip.copy()
			prr.rotate(a1,a2,float(a3))
			#prr.write_image('prr.%0d%0d%0d.mrc'%(a1,a2,a3))
			
			ccf=target.calc_ccf(prr,fp_flag.CIRCULANT,1)
			mean=float(ccf.get_attr("mean"))
			sig=float(ccf.get_attr("sigma"))
			ccf.process_inplace("mask.zeroedge3d",{"x0":edge,"x1":edge,"y0":edge,"y1":edge,"z0":edge,"z1":edge})
			sum+=ccf
			ccf.process_inplace("mask.onlypeaks",{"npeaks":0})		# only look at peak values in the CCF map
			#ccf.write_image('ccf.%0d%0d%0d.mrc'%(a1,a2,a3))
			vec=ccf.calc_highest_locations(mean+sig)
			
			for v in vec: best.append([v.value,a1,a2,a3,v.x-tdim2[0]/2,v.y-tdim2[1]/2,v.z-tdim2[2]/2,0])

#			print a1,a2,a3,mean+sig,float(ccf.get_attr("max")),len(vec)
	
	best.sort()		# this is a list of all reasonable candidate locations

	best.reverse()
	
	if len(best)<1:
		cm=target.calc_center_of_mass(0)
		best.append([0,0,0,0,cm[0]-tdim2[0]/2,cm[1]-tdim2[1]/2,cm[2]-tdim2[2]/2,0])
	print len(best)," possible candidates"

	# this is designed to eliminate angular redundancies in peak location
	print best[0]
	print best[-1]
	if len(best)>10000:
		best=best[0:10000]
	#print best
	for ii in range(len(best)):
		for jj in range(ii+1,len(best)):
			i=best[ii]
			j=best[jj]
			if (i[4]-j[4])**2+(i[5]-j[5])**2+(i[6]-j[6])**2>8.8 : continue
			if j[0]==i[0] : i[7]=1
	for i in best:
		for j in best:
			if (i[4]-j[4])**2+(i[5]-j[5])**2+(i[6]-j[6])**2>8.8 : continue
			if j[0]>i[0] : i[7]=1
			
	
	best2=[]
	for i in best:
		if not i[7]: best2.append([i[0],i[1],i[2],i[3],i[4]*sfac,i[5]*sfac,i[6]*sfac,i[7]])

	# now we find peaks in the sum of all CCF calculations, and find the best angle associated with each peak
	#sum.process_inplace("mask.onlypeaks",{"npeaks":0})
	#sum.write_image("sum.mrc")
	#vec=sum.calc_highest_locations(mean+sig+.0000001)
	#best2=[]
	#for v in vec:
		#print "%5.1f  %5.1f  %5.1f"%(v.x*sfac-tdim[0]/2,v.y*sfac-tdim[1]/2,v.z*sfac-tdim[2]/2)
		#for i in best:
			#if i[4]+tdim2[0]/2==v.x and i[5]+tdim2[1]/2==v.y and i[6]+tdim2[2]/2==v.z :
				#best2.append([i[0],i[1],i[2],i[3],i[4]*sfac,i[5]*sfac,i[6]*sfac,i[7]])
				#break

	best2.sort()
	best2.reverse()
	best2=best2[0:options.num]
	print len(best2), " final candidates"
	print "Qual     \talt\taz\tphi\tdx\tdy\tdz\t"
	for i in best2: 
		print "%1.5f  \t%1.3f\t%1.3f\t%1.3f\t%1.1f\t%1.1f\t%1.1f"%(-i[0],i[1],i[2],i[3],i[4],i[5],i[6])
	#exit()
	# try to improve the angles for each position
	print "\nOptimize each candidate in the reduced map with multiple angle trials"
	print "Qual     \talt\taz\tphi\tdx\tdy\tdz\t"
	cmp_target=target
	cmp_probe=probe
	for j in range(len(best2)):
		print j," --------"
		tries=[[0,0],[0,0],[0,0],[0,0]]
		testang=((0,0),(180.0,0),(0,180.0),(180.0,180.0))	# modify the 'best' angle a few times to try to find a better minimum
		for k in range(4):
			guess=best2[j][1:7]
			guess[0]+=testang[k][0]
			guess[1]+=testang[k][1]
			sm=Simplex(compares,guess,[15,15,15,5,5,5])
			m=sm.minimize(monitor=0,epsilon=.01)
			tries[k][0]=m[1]
			tries[k][1]=m[0]
			print "%1.3f  \t%1.2f\t%1.2f\t%1.2f\t%1.1f\t%1.1f\t%1.1f"%(-tries[k][0],tries[k][1][0],tries[k][1][1],tries[k][1][2],
				tries[k][1][3],tries[k][1][4],tries[k][1][5])
		best2[j][1:7]=min(tries)[1]		# best of the 4 angles we started with
	
	# reread the original images
	target.read_image(args[0])
	probe.read_image(probefilename)
	probe.process_inplace("normalize.unitsum")
	probe.mult(10000)
	
	cmp_target=target
	cmp_probe=probe
	
#	for i in best2:
#		c=probe.get_clip(Region((pdim[0]-tdim[0])/2,(pdim[1]-tdim[1])/2,(pdim[2]-tdim[2])/2,tdim[0],tdim[1],tdim[2]))
#		c.rotate_translate(*i[1:7])
#		c.write_image("z.%02d.mrc"%best2.index(i))
	
	print "Final optimization of each candidate"
	final=[]
	for j in range(len(best2)):
		sm=Simplex(compare,best2[j][1:7],[.5,.5,.5,2.,2.,2.])
		bt=sm.minimize(epsilon=options.epsilon)
		b=bt[0]
		print "\n%1.2f\t(%5.2f  %5.2f  %5.2f    %5.1f  %5.1f  %5.1f)"%(-bt[1],b[0],b[1],b[2],b[3],b[4],b[5])
		final.append((bt[1],b))
	
	print "\n\nFinal Results"
	print "Qual     \talt\taz\tphi\tdx\tdy\tdz\t"
	out=open("foldfitter.out","w")
	final.sort()
	for i,j in enumerate(final):
		b=j[1]
		print "%d. %1.3f  \t%1.2f\t%1.2f\t%1.2f\t%1.1f\t%1.1f\t%1.1f"%(i,-j[0],b[0],b[1],b[2],b[3],b[4],b[5])
		out.write("%d. %1.3f  \t%1.2f\t%1.2f\t%1.2f\t%1.1f\t%1.1f\t%1.1f\n"%(i,-j[0],b[0],b[1],b[2],b[3],b[4],b[5]))
		
		t=Transform()
		#t.set_pre_trans((b[3]+tdim[0]/2,b[4]+tdim[1]/2,b[5]+tdim[2]/2),b[0],b[1],b[2],(0,0,0))
		t.set_pre_trans((b[3]+tdim[0]/2,b[4]+tdim[1]/2,b[5]+tdim[2]/2))
		t.set_rotation({'type':'eman', 'az':b[0], 'alt':b[1], 'phi':b[2]})
		#t.set_trans((0,0,0))
		#t.set_trans((b[3]+tdim[0]/2,b[4]+tdim[1]/2,b[5]+tdim[2]/2))
		
		s=Transform()
		t=Transform()
		s.set_rotation({'type':'eman', 'az':b[0], 'alt':b[1], 'phi':b[2]})
		s.set_trans((b[3],b[4],b[5]))
		print (pdim[0]-tdim[0])/2,(pdim[1]-tdim[1])/2,(pdim[2]-tdim[2])/2
		pc=probe.get_clip(Region((pdim[0]-tdim[0])/2,(pdim[1]-tdim[1])/2,(pdim[2]-tdim[2])/2,tdim[0],tdim[1],tdim[2]))
		pc.transform(s)
		if target["MRC.nxstart"]==0 and target["MRC.nystart"]==0:
			shx=target['origin_x']
			shy=target['origin_y']
			shz=target['origin_z']
		else:
			shx=target["MRC.nxstart"]*target["apix_x"]
			shy=target["MRC.nystart"]*target["apix_x"]
			shz=target["MRC.nzstart"]*target["apix_x"]
		b[3]=b[3]*apix-pc['origin_x']+shx
		b[4]=b[4]*apix-pc['origin_y']+shy
		b[5]=b[5]*apix-pc["origin_z"]+shz
		#pc['origin_x']=target["MRC.nxstart"]*target["apix_x"]
		#pc['origin_y']=target["MRC.nystart"]*target["apix_x"]
		#pc["origin_z"]=target["MRC.nzstart"]*target["apix_x"]
		#pc.write_image("tst.mrc")
		t.set_rotation({'type':'eman', 'az':b[0], 'alt':b[1], 'phi':b[2]})
		t.set_trans((b[3],b[4],b[5]))
		pdb_transform(args[1],'final.%02d.pdb'%i,t,cent)


	print ncmp," total comparisons"
	out.close()
	
#	print compare(best2[0][1:7])
	
#	print best2[0]
#	print best2[-1]
		
if __name__ == "__main__":
    main()
