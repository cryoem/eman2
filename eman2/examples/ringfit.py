#!/usr/bin/env python
# 01/13/2014		Steven Ludtke
# This program tries to fit a closed loop of blobs to DNA-minicircle density pattern using
# a simple distance, angle, dihedral potential with a closed linear chain of balls

from EMAN2 import *
import sys
import numpy.linalg as LA
import random
import math
from numpy import *
import os
from os import listdir
from os.path import isfile, join

def main():
	#process_image(sys.argv[1])

	mypath=sys.argv[1]
	
	# Process a single mrc file
	if mypath.endswith(".mrc"):
		process_image(mypath,"")
	# Process all mrc files in a folder
	else:
		files = [ f for f in listdir(mypath) if f.endswith(".mrc")]
		#print files
		for fname in files:
			print join(mypath,fname)
			shape=process_image(join(mypath,fname),fname)
			
			# The statistics output part should be replaced by ringstat.py
			outfile=open("output","a")
			outfile.write("%s\t%1.3g\t%1.3g\t%1.3g\n"%(fname,shape[0],shape[1],shape[2]))
			outfile.close()
	
def process_image(imgname,imgprefix):
	p=EMData(imgname,0)
	print "Preprocessing"
	# This will (hopefully) isolate the mini-circle
	##p.process_inplace("normalize.edgemean")
	##p.write_image("000.mrc",0)
	##p.process_inplace("mask.dust3d",{"voxels":6000,"threshold":3})
	##p.process_inplace("threshold.belowtozero",{"minval":3})
	##p.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
	##p.process_inplace("normalize.edgemean")
	##p.process_inplace("mask.dust3d",{"voxels":30000,"threshold":3})
	##p.process_inplace("threshold.belowtozero",{"minval":3})
	##p.process_inplace("threshold.binary",{"value":p["mean"]+p["sigma"]*5.0})
	##p.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
	##p.process_inplace("normalize.unitsum")
	##p.mult(10000.0)
	##p.write_image("aaa.mrc",0)
	##p=EMData("aaa.mrc",0)
	##p.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.05})
	###p.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
	##p.process_inplace("normalize.unitsum")
	##p.mult(10000.0)
	##p.mult(30.0)
	##p.write_image("bbb.mrc",0)
	
	
	
	
	p.process_inplace("normalize.edgemean")
	p.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
	p.process_inplace("normalize.edgemean")
	p.process_inplace("mask.dust3d",{"voxels":4000,"threshold":2})
	p.process_inplace("mask.dust3d",{"voxels":2000,"threshold":3})
	p.write_image("aaa.mrc",0)
	p.process_inplace("mask.auto3d",{"nshells":2,"nshellsgauss":3,"radius":35,"threshold":1.6})
	p.process_inplace("threshold.belowtozero",{"minval":0})
	p.process_inplace("threshold.binary",{"value":p["mean"]+p["sigma"]*5.0})
	p.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
	p.process_inplace("normalize.unitsum")
	p.mult(10000.0)
	# Align the circle for setting initial points
	print "Inertia matrix and alignment"
	# compute the resulting inertia matrix
	an=Analyzers.get("inertiamatrix",{"verbose":0})
	an.insert_image(p)
	mxi=an.analyze()
	mx=EMNumPy.em2numpy(mxi[0])

	# Compute the eigenvalues/vectors
	eigvv=LA.eig(mx)		# a 3-vector with eigenvalues and a 3x3 with the vectors
	if min(eigvv[0])==0 :
		print "error on ",pf
		drn-=1
		sys.exit(1)

	eig=[(1.0/eigvv[0][i],eigvv[1][:,i]) for i in xrange(3)]  # extract for sorting
	eig=sorted(eig)		# now eig is sorted in order from major to minor axes
	#T=array([eig[0][1],eig[1][1],eig[2][1]])            # reassemble sorted matrix

	T=Transform((float(i) for i in (eig[0][1][0],eig[0][1][1],eig[0][1][2],0,eig[1][1][0],eig[1][1][1],eig[1][1][2],0,eig[2][1][0],eig[2][1][1],eig[2][1][2],0)))


	m=EMData(p)

	
	print "Initializing start points..."
	# Calculate the length of three axis of the circle
	SX=m.get_xsize()
	SY=m.get_ysize()
	SZ=m.get_zsize()
	xsta=SX;xend=0;ysta=SY;yend=0;zsta=SZ;zend=0
	for i in range(SX):
		for j in range(SY):
			for k in range(SZ):
				if(m.get_value_at(i,j,k)>0.1):
					v=T.transform(i,j,k)
					x=v[0]
					y=v[1]
					z=v[2]
					if(x<xsta): xsta=x
					if(x>xend): xend=x
					if(y<ysta): ysta=y
					if(y>yend): yend=y
					if(z<zsta): zsta=z
					if(z>zend): zend=z
	print xsta,xend,ysta,yend,zsta,zend
	#print SX,SY,SZ
	#print (xsta+xend)/2-SX/2,(ysta+yend)/2-SY/2,(zsta+zend)/2-SZ/2
	iT=T.inverse()
	bestres=10000	# Best result
	for ttt in range(10):	# Run 10 times with random start poingt for best result (and ambiguous score? # todo #)
		numofbp=21	# Maximum number of points (+1 for showing circle in chimera)
		nowbp=5		# Start number
		stepsz=500	# Interval of print status & check stablization
		plen=(nowbp*math.sin(math.pi/nowbp))/(math.pi)	# Penalty factor of length for representing a circle using polygon
		
		#totlen=2*336*3.3*plen # total length
		totlen=336*3.3*plen # total length

		#m=EMData(sys.argv[1],0)
		pa=PointArray()
		pa.set_number_points(nowbp)
		# set initial parameters
		pa.sim_set_pot_parms(totlen/nowbp, .01, .0, 35.9*pi/180.0, 0.0, 8000.0,m,totlen/nowbp*.6,800)
		startrange=40	# Range of random start positions
		startphase=random.uniform(0,2.0*pi/nowbp)	# Random starting phase
		#xsft=random.uniform(-startrange,startrange)
		#ysft=random.uniform(-startrange,startrange)
		#zsft=random.uniform(-startrange,startrange)
		
		# Initializing points on one plane, and then transform to the plane of the circle
		for i in xrange(nowbp):
			ang=2.0*pi*i/nowbp+startphase
		#	pa.set_vector_at(i,Vec3f(0,sin(ang)*100,cos(ang)*100),1.0)
			vx=0 + random.uniform(-startrange,startrange)
			vy=sin(ang)*(yend-ysta)*2.3 + random.uniform(-startrange,startrange) #- ((ysta+yend)*2-SY*2)  
			vz=cos(ang)*(zend-zsta)*2.3 + random.uniform(-startrange,startrange) #- ((zsta+zend)*2-SZ*2)  
			pa.set_vector_at(i,Vec3f(vx,vy,vz),1.0)
		#pa.save_to_pdb("uuu.pdb")
		pa.transform(T)
		pa.save_to_pdb("ttt.pdb")

		now_potential=pa.sim_potential()
		isstable=0	# When the result is stable for 3*stepsz iterations, add points or stop
		skipping=0	# When stop, skip the next iterations
		bestpintrail=100000	# Best potential score in this trail
		
		
		for i in xrange(200000):
			if skipping==0:	# run one simulation step
				pa.sim_minstep_seq(.1)
				if random.uniform(0,1.0)>.9996 and nowbp>10:
					print "swapping................"
					sn=random.randint(0,9)
					s=[sn,sn+9]
					for ii in range(s[0],(s[0]+s[1])/2):
						jj=s[1]+s[0]-ii
						tmp=pa.get_vector_at(ii)
						pa.set_vector_at(ii,pa.get_vector_at(jj),1.0)
						pa.set_vector_at(jj,tmp,1.0)
					pa.save_to_pdb("swap.pdb")

			if isstable>5:
				if nowbp*2<numofbp: # Add points when result is stable and number of current points is lower than the total number
					print "adding points...."
					pa.sim_add_point_double()	# Put one additional point on each edge
					nowbp=nowbp*2
					print nowbp
					plen=(nowbp*math.sin(math.pi/nowbp))/(math.pi)	# Recalculate the length penalty
					totlen=336*3.3*plen
					#totlen=2*336*3.3*plen
					pa.sim_set_pot_parms(totlen/nowbp, .5, 100, 35.9*pi/180.0, 0.0, 800.0,m,totlen/nowbp*.6,10000)
					#pa.sim_set_pot_parms(totlen/nowbp, 1, 150, 35.9*pi/180.0, 0.0, 8000.0,m,totlen/nowbp*.9,10000)
					isstable=0
				else:
					skipping=1

			#if i==500: print "aaa"
			if i%stepsz==0:
				
				print i/stepsz
				pa.sim_printstat()
				
				old_potential=now_potential
				now_potential=pa.sim_potential()
				if(abs((now_potential-old_potential)/old_potential)<.005):
					isstable+=1
					#print "aaa"
				else:
					isstable=0
				
				#print now_potential,bestpintrail
				if (now_potential<bestpintrail and nowbp*2>=numofbp):
					bestpintrail=now_potential
					bestpa=pa.copy()
				
				## output pdb file for md movies
				#fname="%d.pdb"%(i/stepsz)
				#pa.save_to_pdb(fname)
				
				#pdbfile = open(fname, "r")
				#lines = pdbfile.readlines()
				#pdbfile.close()
				
				#panum=pa.get_number_points()
				#if panum<numofbp:
					#outfile=open(fname,"a")
					#for k in range(panum,numofbp):
						#ln=lines[0]
						#s=str(k)+' '
						#if (len(s)<3): s=' '+s
						#ln=ln.replace(' 0 ',s)
						#outfile.write(ln)
					#outfile.close()

		if bestpintrail<bestres:
			bestres=bestpintrail
			fname=imgprefix+"_result.pdb"
			bestpa.save_to_pdb(fname)
			
			pdbfile = open(fname, "r")
			lines = pdbfile.readlines()
			pdbfile.close()
			
			panum=bestpa.get_number_points()
			if panum<numofbp:
				outfile=open(fname,"a")
				for k in range(panum,numofbp):
					ln=lines[0]
					s=str(k)+' '
					if (len(s)<3): s=' '+s
					ln=ln.replace(' 0 ',s)
					outfile.write(ln)
				outfile.close()
			
			bestpa.sim_add_point_double()
			bestpa.sim_add_point_double()
			bestpa.sim_add_point_double()
			img=bestpa.pdb2mrc_by_summation(SX,4.52,10.,-1)
			img.write_image("img.mrc",0)
			
			
			
		#

	print bestres

	#compute the resulting inertia matrix
	img=EMData("img.mrc")
	
	img.process_inplace("mask.addshells.gauss",{"val1":10,"val2":30})
	img.write_image("maskimg.mrc",0)
	finalimg=EMData(imgname,0)
	#finalimg=EMData("0143r_ori.mrc",0)
	#finalimg.process_inplace("xform.flip",{"axis":"y"})
	finalimg.process_inplace("mask.fromfile",{"filename":"maskimg.mrc","ismaskset":0})
	finalimg.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})

	finalimg.process_inplace("normalize.edgemean")
	finalimg.process_inplace("mask.dust3d",{"voxels":4000,"threshold":4})
	finalimg.process_inplace("threshold.belowtozero",{"minval":0})
	finalimg.write_image(imgprefix+"_finalimg.mrc",0)
	an=Analyzers.get("inertiamatrix",{"verbose":0})
	an.insert_image(finalimg)
	mxi=an.analyze()
	mx=EMNumPy.em2numpy(mxi[0])


	# Compute the eigenvalues/vectors
	eigvv=LA.eig(mx)		# a 3-vector with eigenvalues and a 3x3 with the vectors
	eig=[(1.0/eigvv[0][i],eigvv[1][:,i]) for i in xrange(3)]  # extract for sorting
	eig=sorted(eig)		# now eig is sorted in order from major to minor axes
	#T=array([eig[0][1],eig[1][1],eig[2][1]])            # reassemble sorted matrix

	T=Transform((float(i) for i in (eig[0][1][0],eig[0][1][1],eig[0][1][2],0,eig[1][1][0],eig[1][1][1],eig[1][1][2],0,eig[2][1][0],eig[2][1][1],eig[2][1][2],0)))


	finalimg.transform(T)

	#finalimg.write_image("trans_finalimg.mrc",0)
	# now the shape is aligned to Z/Y/X so the greatest axial extent should be along Z
	an=Analyzers.get("shape",{"verbose":0})
	an.insert_image(finalimg)
	shp=an.analyze()[0]

	# Z/Y - should always be >1, Y/X, Z/X
	#out.write("%1.3g\t%1.3g\t%1.3g\t# %s\n"%(shp[2]/shp[1],shp[1]/shp[0],shp[2]/shp[0],pf.split("/")[-1]))
	shape=sorted([abs(shp[0]),abs(shp[1]),abs(shp[2])])
	print shape

	print "%1.3g\t%1.3g\t%1.3g\t#"%(shape[2]/shape[1],shape[1]/shape[0],shape[2]/shape[0])
	return [shape[2]/shape[1],shape[1]/shape[0],shape[2]/shape[0]]


if __name__ == '__main__':
    main()