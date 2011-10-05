#!/usr/bin/env python

#
# Author: Matthew Baker, 06/30/2009 (mbaker@bcm.edu)
# Copyright (c) 2000-2009 Baylor College of Medicine
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
import numpy

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <volume file> <Angstroms per pixel> <resolution> <threshold>
	
	WARNING: This program is still under development

	Identifies alpha helices and beta sheets in maps at subnanometer resolutions"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--atoms", "-A", type=str, help="pseudoatoms file", default="none")
	parser.add_argument("--atomswt", "-P", type=float, help="pseudoatom weight", default=1.0)
	parser.add_argument("--coeff","-C", type=str,help="helix correlation file",default="none")
	parser.add_argument("--coeffwt","-H", type=float,help="helix correlation weight",default=1.0)
	parser.add_argument("--skeleton","-S", type=str,help="skeleton file",default="none")
#	parser.add_argument("--sketype","-T", type=str,help="skeleton type",default="none")
	parser.add_argument("--skeletonwt","-W", type=float,help="skeleton weight",default=1.0)
	parser.add_argument("--helixlength","-L", type=float,help="helix length om angstroms",default=16.2)
	parser.add_argument("--da","-D", type=float,help="helix angular search step",default=5.0)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()
		
	if len(args) < 4 : 
		parser.error("ERROR: supply all 4 required arguments!")

	logid=E2init(sys.argv,options.ppid)

	apix = float(args[1]) #Angstroms per pixel
	res = float(args[2]) #Resolution
	thr = float(args[3]) #Threshold
	helixsize = 3
	sheetsize = 4
	
	if options.coeffwt > 1 or options.coeffwt < 0:
		print "Helix correlation weight must be between 0 and 1"
		sys.exit()
	if options.skeletonwt > 1 or options.skeletonwt < 0:
		print "Skeleton weight must be between 0 and 1"
		sys.exit()
	if options.atomswt > 1 or options.atomswt < 0:
		print "Pseudoatom weight must be between 0 and 1"
		sys.exit() 
	
	try : 
		infile=open(args[0],"r")
		infile.close()
	except : 
		parser.error("Cannot open input file")
	
	target_volume_filepath = args[0]
	target_volume_name = os.path.basename(target_volume_filepath)
	target_volume_name = os.path.splitext(target_volume_name)[0]
	
	targetMRC=EMData()
	targetMRC.read_image(target_volume_filepath)
	(nx,ny,nz) = ( targetMRC.get_xsize(), targetMRC.get_ysize(), targetMRC.get_zsize() )
	if not (nx == ny and ny == nz and nx % 2 == 0):
		print "The density map must be even and cubic. Terminating SSEHunter."
		sys.exit()	
	
	# getting pseudoatoms
	if options.atoms=="none":
		patoms = generate_pseudoatoms(targetMRC.copy(), apix, res, thr)
		atomNumber = range(1, 1+len(patoms))
	else:
		(patoms, atomNumber) = read_pseudoatoms(options.atoms)
	
	# getting individual scores
	skeletonArray = skeleton_scores(patoms, targetMRC, thr, helixsize, sheetsize, options.skeleton)
	(coeffArray, avghhvalue, maxValue) = helix_correlation_scores(targetMRC, patoms, apix, res, options.coeff, options.helixlength, options.da)
	pseudoatomArray = geometry_scores(patoms, atomNumber, targetMRC, apix, thr, coeffArray, avghhvalue)

	# output composite scores
	write_composite_scores(target_volume_name, apix, patoms, atomNumber, coeffArray, options.coeffwt, avghhvalue, maxValue, skeletonArray, options.skeletonwt, pseudoatomArray, options.atomswt)

	E2end(logid)

def cross_product(a,b):
	"""Cross product of two 3-d vectors. from http://www-hep.colorado.edu/~fperez/python/python-c/weave_examples.html
	"""
	cross = [0]*3
	cross[0] = a[1]*b[2]-a[2]*b[1]
	cross[1] = a[2]*b[0]-a[0]*b[2]
	cross[2] = a[0]*b[1]-a[1]*b[0]
	return cross

def update_map(target, location,rangemin,rangemax):
	rmin=int(round(rangemin/2))
	rmax=int(round(rangemax/2))+1
	maxdistance=sqrt(3*rmin**2)
	for x in (range(rmin,rmax)):
		for y in (range(rmin,rmax)):
			for z in (range(rmin,rmax)):
				temp_value=target.get_value_at(location[0]+x,location[1]+y,location[2]+z)
				distance=sqrt(x**2+y**2+z**2)
				if x==0 and y==0 and z==0:
					pixel_value=0.0
				else:
					pixel_value=temp_value*(distance/maxdistance)
				target.set_value_at(location[0]+x,location[1]+y,location[2]+z, pixel_value)
				target.update()

def generate_pseudoatoms(target, apix, res, thr):
	print "No psuedoatom input; generating Pseudoatoms"
	patoms = []
	rangemin=-1*res/apix
	rangemax=res/apix
	out = open("pseudoatoms.pdb","w")
	Max_value = target["maximum"]
	Max_location = target.calc_max_location()
	i=1
	chain="A"
	point_array=[]
	while Max_value >= thr:
		patoms.append(Max_location)
		out.write("ATOM  %5d  C   GLY %s%4d    %8.3f%8.3f%8.3f  1.00     0      S_00  0 \n"%(i,chain,i,Max_location[0]*apix,Max_location[1]*apix,Max_location[2]*apix))
		update_map(target, Max_location, rangemin,rangemax)	
		Max_location = target.calc_max_location() #(x,y,z) for maximum value
		Max_value = target[Max_location]
		i=i+1
	out.close()
	return patoms

def read_pseudoatoms(pdb_filepath):
	print "Reading pseudatoms"
	patoms=[]
	atomNumber=[]
	
	file=open(atoms, "r")
	for line in file:
		if str(line[0:6].strip())=="ATOM":
			TempX=(float(line[30:38].strip()))/apix
			TempY=(float(line[38:46].strip()))/apix
			TempZ=(float(line[46:54].strip()))/apix
			patom=[TempX,TempY,TempZ]
			patoms.append(patom)
			atomNumber.append(int(line[6:11].strip()))
	
	file.close()
	return (patoms, atomNumber)
	

def get_angle(patoms, origin, p1, p2):
	v1=(patoms[p1][0]-patoms[origin][0],patoms[p1][1]-patoms[origin][1],patoms[p1][2]-patoms[origin][2])
	lengthV1=sqrt(v1[0]**2+v1[1]**2+v1[2]**2)
	v1n=(v1[0]/lengthV1,v1[1]/lengthV1,v1[2]/lengthV1)
	v2=(patoms[p2][0]-patoms[origin][0],patoms[p2][1]-patoms[origin][1],patoms[p2][2]-patoms[origin][2])
	lengthV2=sqrt(v2[0]**2+v2[1]**2+v2[2]**2)
	v2n=(v2[0]/lengthV2,v2[1]/lengthV2,v2[2]/lengthV2)

	dp = v1n[0]*v2n[0] + v1n[1]*v2n[1] + v1n[2]*v2n[2] #dot product
	if dp > 1:
		dp=1
	if dp<-1:
		dp=-1
	angle=acos(dp)*(180/pi)
	if angle>90:
		angle=180-angle
	return angle
def get_distance_matrix(patoms, apix):
	### all to all distance calculations
	print "Calculating all atom distances"
	distance=[]
	
	for atom1 in patoms:
		x1=atom1[0]
		y1=atom1[1]
		z1=atom1[2]
		d1=[]
	
		for atom2 in patoms:
			x2=atom2[0]
			y2=atom2[1]
			z2=atom2[2]
			dist=(sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))*apix
			d1.append(dist)
		distance.append(d1)
	return distance

def find_aspect_ratio(targetMRC, patoms, apix, thr):
	### set search area
	print "Grouping neighboring pixels"
	kernelwidth=int(round(5.0/apix))
	pixels=[]
	for atom in patoms:
		tempPixel=[]
		pixelCoords=[]
		intX=int(atom[0])
		intY=int(atom[1])
		intZ=int(atom[2])
		for tempX in range(intX-kernelwidth, intX+kernelwidth, 1):
			for tempY in range(intY-kernelwidth, intY+kernelwidth, 1):
				for tempZ in range(intZ-1*kernelwidth, intZ+kernelwidth, 1):
					if targetMRC.get_value_at(tempX,tempY,tempZ) > thr:
						tempPixel = [tempX, tempY, tempZ]
						if tempPixel not in pixelCoords:
							pixelCoords.append(tempPixel)
		pixels.append(pixelCoords)
	
	### Check 2nd moments to see if globular
	print "Assessing neighborhood shape"
	axisMatrix=[]
	for index1 in range(len(patoms)):
		NumPoints1=numpy.array(pixels[index1],'d')
		NumPoints1_mean=numpy.sum(NumPoints1,axis=0)/len(NumPoints1)
		NumPoints2=NumPoints1-NumPoints1_mean
		h = numpy.sum(map(numpy.outer,NumPoints2,NumPoints2),axis=0)
		[u1,x1,v1]=numpy.linalg.svd(h)
		if x1.all()==0:
			print index1,
			print "is bad"
			xmod=x1
		else:
			xmod=x1/max(x1)
			
		if xmod[2]==0:
			aspectratio=0
		else:
			aspectratio=xmod[1]/xmod[2]
			print aspectratio
		axisMatrix.append(aspectratio)
		
	return axisMatrix

def model_pdb_helix(length):
	print "Generating prototypical helix %f Angstroms in length"%length
	helix=PointArray()
	j=0
	p=0
	q=0
	r=0
	outfile=open("helix.pdb","w")
	while j<= round(length/1.54):
		Nxcoord=cos((100*j*pi)/180)*1.6
		Nycoord=sin((100*j*pi)/180)*1.6
		Nzcoord=j*1.52
		CAxcoord=cos(((28+(100*j))*pi)/180)*2.3
		CAycoord=sin(((28+(100*j))*pi)/180)*2.3
		CAzcoord=(j*1.54)+.83
		Cxcoord=cos(((61+(100*j))*pi)/180)*2.0
		Cycoord=sin(((61+(100*j))*pi)/180)*2.0
		Czcoord=(j*1.54)+1.7
		Oxcoord=cos(((61+(100*j))*pi)/180)*2.0
		Oycoord=sin(((61+(100*j))*pi)/180)*2.0
		Ozcoord=(j*1.54)+3.09
		p=j*4+1
		q=j*4+2
		r=j*4+3
		s=j*4+4
		j=j+1
		
		outfile.write("ATOM  %5d   N  GLY A%4d	%8.3f%8.3f%8.3f  1.00	 0	  S_00  0 \n"%(p,j,Nxcoord,Nycoord,Nzcoord))
		outfile.write("ATOM  %5d  CA  GLY A%4d	%8.3f%8.3f%8.3f  1.00	 0	  S_00  0 \n"%(q,j,CAxcoord,CAycoord,CAzcoord))
		outfile.write("ATOM  %5d   C  GLY A%4d	%8.3f%8.3f%8.3f  1.00	 0	  S_00  0 \n"%(r,j,Cxcoord,Cycoord,Czcoord))
		outfile.write("ATOM  %5d   O  GLY A%4d	%8.3f%8.3f%8.3f  1.00	 0	  S_00  0 \n"%(s,j,Oxcoord,Oycoord,Ozcoord))
	outfile.close()
	helix.read_from_pdb("helix.pdb")
	return helix

	
def model_mrc_helix(box_size, apix, res, length = 10.8, points = None, helixtype = "helix_pdb"):
	mrcHelix=EMData(box_size,box_size,box_size)
	helixtypeDict = {"gauss":0, "gauss_falloff":1, "polynomial":2, "helix_pdb":3, 0:0, 1:1, 2:2, 3:3}
	
	if helixtypeDict[helixtype] == 3:
		mrcHelix=points.pdb2mrc_by_summation(box_size,apix,res)
		mrcHelix.process_inplace("normalize.edgemean")
		mrcHelix.process_inplace("xform.centerofmass")
		
		mrcHelix.process_inplace("math.poly_radial_profile", {"length": length})
		mrcHelix.write_image("pdbhelix_polynomial.mrc")

	else:
		mrcHelix["apix_x"] = apix
		mrcHelix["apix_y"] = apix
		mrcHelix["apix_z"] = apix
		mrcHelix.process_inplace("math.model_em_cylinder", {"type": helixtypeDict[helixtype],"length": length})
	
	
#	mrcHelix=points.pdb2mrc_by_summation(box,apix,res)
#	mrcHelix.process_inplace("normalize.edgemean")
#	mrcHelix.process_inplace("xform.centerofmass")
#	aout=mrcHelix.copy()
#	aout.to_zero()
#	mrcHelix.process_inplace("filter.highpass.gauss",{"apix":apix,"cutoff_freq":100})

#	for i in range(box):
#		r = Region(0,0,i,box,box,1)
#		slice1 = mrcHelix.get_clip(r)
#		thresh = slice1.process("threshold.binary", {"value":10})
#		neg = thresh.process("mask.addshells",{"nshells":2})
#		real_neg = neg-thresh
#		real_neg.mult(-10)
#		pos = neg.process("mask.addshells",{"nshells":1})
#		real_pos = pos-neg
#		real_pos.mult(10)
#		solution = slice1.copy()
#		solution.mult(thresh)
#		solution += real_neg
#		solution += real_pos
#		aout.insert_clip(solution,[0,0,i])
#	aout.write_image("helix.mrc")
#	return aout
#	mrcHelix.write_image("helix.mrc")
#	mrcHelix.read_image("cyl.mrc",-1)

	return mrcHelix
	
def helixhunter_ccf(target, probe, da):
	print "Running helix correlation routine"
	bestCCF= target.copy() #To get metadata
	bestCCF.to_zero() #Clearing the volume
	
	s = Symmetries.get("c1")
	orients = s.gen_orientations("eman",{"delta":da,"inc_mirror":True})
	counter=1
	N = len(orients)
	for i in range(N):
		if i % 25 == 0:
			print "%5.2f" % (i*100.0/N) + " %"
		t = orients[i]
		probeMRC= probe.process("xform",{"transform":t}) 
		currentCCF=target.calc_mutual_correlation(probeMRC)
		bestCCF.process_inplace("math.max",{"with":currentCCF})

	bestCCF.write_image("int-hh-coeff.mrc")
	return bestCCF

def helix_correlation_scores(targetMRC, patoms, apix, res, coeff, helixlength, da):
	####This needs to be replaced 
	
	if (coeff=="none") :
		pdbHelix=model_pdb_helix(helixlength)	
		mrcHelix=model_mrc_helix(targetMRC.get_xsize(), apix, res, length=helixlength, points=pdbHelix, helixtype="helix_pdb")
		if da==0.0:
			da=2*asin(res/targetMRC.get_xsize())*(180.0/pi)
			print "da not set; setting da to %f degrees"%da
		hhMrc=helixhunter_ccf(targetMRC, mrcHelix, da)
	
	else:	
		hhMrc.read_image(coeff,-1)
	
	#hhMrc.set_attr_dict(targetMRC.get_attr_dict())
	coeffArray=[]
	avghhvalue=0.0
	maxValue=0.0
	atomCount = len(patoms)
	for patom in patoms: 
		hhvalue=float(hhMrc.get_value_at(patom[0],patom[1],patom[2]))
		if hhvalue > maxValue:
			maxValue=hhvalue
		avghhvalue=hhvalue+avghhvalue
		coeffArray.append(hhvalue)
	avghhvalue=avghhvalue/atomCount
	print "Correlation threshold: %f"%(avghhvalue)
	print "Correlation Maximum:   %f"%(maxValue)
	hhMrc.write_image("hhMrc.mrc")
	return (coeffArray, avghhvalue, maxValue)
	
def geometry_scores(patoms, atomNumber, targetMRC, apix, thr, coeffArray, avghhvalue):
	'''psuedoatom geometry calculations'''
	
	atomCount=len(patoms)
	distance = get_distance_matrix(patoms, apix) #distances between pseudoatoms
	pseudoatomArray=[]
	neighborhood=[]
	pointcloud=[]
	betadistance=8.25
	
	axisMatrix = find_aspect_ratio(targetMRC, patoms, apix, thr)
	
	for index1 in range(atomCount):
		### checks for nearest atoms and nearest helical atoms
		### Identifying related pseudoatoms
		
		################### Note by Ross #####################
		# coeffArray stores the cross-correlation based scores for the pseudoatoms
		# It seems that helixneighbor1 is the closest pseudoatom which is likely to be a helix atom based on coeffArray
		# The distance to helixneighbor1 is mindistance1
		# helixneighbor2 is the second closest and its distance is mindistance2
		# neighbor3 is the closest neighbor regardless of the cross-correlation score, and mindistance3 is the distance to it
		# neighbor4 is the second closest neighbor, and mindistance4 is the distance to it
		# Note: distances to nighbors must be at least 0.1
		#######################################################
		
		mindistance1=99998
		mindistance2=99999
		helixneighbor1=0
		helixneighbor2=0
		mindistance3=99998
		mindistance4=99999
		neighbor3=0
		neighbor4=0
		cloud=[]
		for index2 in range(atomCount):
			if distance[index1][index2]<=betadistance and index1!=index2:
				cloud.append(atomNumber[index2])
			if distance[index1][index2] <= mindistance2 and distance[index1][index2]>0.1 and coeffArray[index2]>=avghhvalue:
				if distance[index1][index2] <= mindistance1:
					mindistance2=mindistance1
					mindistance1=distance[index1][index2]
					helixneighbor2=helixneighbor1
					helixneighbor1=atomNumber[index2]
				else:
					mindistance2=distance[index1][index2]
					helixneighbor2=atomNumber[index2]
	
			if distance[index1][index2] <= mindistance4 and distance[index1][index2]>0.1:
				if distance[index1][index2] <= mindistance3:
					mindistance4=mindistance3
					mindistance3=distance[index1][index2]
					neighbor4=neighbor3
					neighbor3=atomNumber[index2]
				else:
					mindistance4=distance[index1][index2]
					neighbor4=atomNumber[index2]
		neighbors=(helixneighbor1,helixneighbor2,mindistance1,mindistance2,neighbor3, neighbor4)
		pointcloud.append(cloud)
		neighborhood.append(neighbors)
	
	sheetlist=[]
	generallist=[]
	helixlist=[]
	
	### Checking local gemoetry
	for index3 in range(atomCount):
	### check generic angles
		origin=index3
		p1=neighborhood[index3][4]-1
		p2=neighborhood[index3][5]-1
		genericAngle=get_angle(patoms,origin,p1,p2)
	
	### checks helix angles
		p1=neighborhood[index3][0]-1
		p2=neighborhood[index3][1]-1
		helixAngle=get_angle(patoms, origin, p1, p2)

	###checks sheet angles
		cloud=pointcloud[index3]
		arrayofxp=[]
		for firstpoint in cloud:
			point1=firstpoint-1
			pv1=(patoms[point1][0]-patoms[origin][0],patoms[point1][1]-patoms[origin][1],patoms[point1][2]-patoms[origin][2])
			lengthPV1=sqrt(pv1[0]**2+pv1[1]**2+pv1[2]**2)
			pv1n=(pv1[0]/lengthPV1,pv1[1]/lengthPV1,pv1[2]/lengthPV1)
	
			for secondpoint in cloud:
				point2=secondpoint-1
				if point2 != point1 and point2 != origin:
					pv2=(patoms[point2][0]-patoms[origin][0],patoms[point2][1]-patoms[origin][1],patoms[point2][2]-patoms[origin][2])
					lengthPV2=sqrt(pv2[0]**2+pv2[1]**2+pv2[2]**2)
					pv2n=(pv2[0]/lengthPV2,pv2[1]/lengthPV2,pv2[2]/lengthPV2)
					xp=cross_product(pv1,pv2)
					lengthxp=sqrt(xp[0]**2+xp[1]**2+xp[2]**2)
					if lengthxp>0:
						xpn=(xp[0]/lengthxp, xp[1]/lengthxp, xp[2]/lengthxp)
					else:
						xpn=(xp[0], xp[1], xp[2])
					arrayofxp.append(xpn)
		dpxpcounter=0
		dpxpsum=0

		if len(arrayofxp) >=2:
			for dpxp1 in arrayofxp:
				for dpxp2 in arrayofxp:
					if dpxp1 != dpxp2:
						dpxp = dpxp1[0]*dpxp2[0] + dpxp1[1]*dpxp2[1] + dpxp1[2]*dpxp2[2] #dot product
						if dpxp > 1:
							dpxpAngle=0
						elif dpxp<-1:
							dpxpAngle=180
						else:
							dpxpAngle=acos(dpxp)*(180/pi)
						if dpxpAngle>90:
							dpxpAngle=180-dpxpAngle
	
						dpxpcounter=dpxpcounter+1
						dpxpsum=dpxpsum+dpxpAngle
			if dpxpsum==0 and dpxpcounter==0:
				betaAngle=0
			else: 
				betaAngle=dpxpsum/dpxpcounter
		else:
			betaAngle=90	
	
		#generallist.append(genericAngle)
		#sheetlist.append(betaAngle)
		#helixlist.append(helixAngle)
		aspectratio=axisMatrix[index3]
		pascore=0
		print "%d: axis: %f, neighbor angle: %f, helix angle: %f, sheet angle: %f,  number of neighbors: %d"%(atomNumber[index3], aspectratio, genericAngle, helixAngle, betaAngle, len(cloud))
		if aspectratio <=2:
			pascore=pascore+1
		if aspectratio >=3:
			pascore=pascore-1
		if genericAngle <=40:
			pascore=pascore+1
		if genericAngle >=50:
			pascore=pascore-1		
		if helixAngle <=45 and mindistance1<12:
			pascore=pascore+0.5
		if helixAngle <=45 and mindistance2<12:
			pascore=pascore+0.5
		if betaAngle >=30 and betaAngle !=90:
			pascore=pascore-1
		if len(cloud) <= 3:
			pascore=pascore+1
		if len(cloud) > 3:
			pascore=pascore-1
		pseudoatomArray.append(float(pascore/4.0))	


	return pseudoatomArray
##################################################################


#skeleton_scores() by Ross Coleman
def skeleton_scores(pseudoatoms, vol, threshold, helix_size, sheet_size, skeleton):
	"""computes skeleton scores as in EMAN1 skeleton.C"""
	if skeleton != "none":
		print "Add the code to read a skeleton file here!"
	skeleton = vol.process("gorgon.binary_skel", {"threshold":threshold, "min_curve_width":helix_size, "min_surface_width":sheet_size, "mark_surfaces":True})
	skeletonScores = [] #corresponds with pseudoatoms
	
	MAX_DISTANCE = 4
	
	for (x,y,z) in pseudoatoms:
		low_x = max(0, x - MAX_DISTANCE)
		low_y = max(0, y - MAX_DISTANCE)
		low_z = max(0, z - MAX_DISTANCE)
		high_x = min(vol["nx"], x + MAX_DISTANCE + 1)
		high_y = min(vol["ny"], y + MAX_DISTANCE + 1)
		high_z = min(vol["nz"], z + MAX_DISTANCE + 1)
		count = 0
		curve_score = 0
		surface_score = 0
		for k in range(low_z, high_z):
			for j in range(low_y, high_y):
				for i in range(low_x, high_x):
					val = int(vol[(i,j,k)])
					i2 = i - x
					j2 = j - y
					k2 = k - z
					if val == 1:
						curve_score += 1 - (i2**2+j2**2+k2**2)/(3*MAX_DISTANCE**2)
						count += 1
					elif val == 2:
						surface_score += 1 - (i2**2+j2**2+k2**2)/(3*MAX_DISTANCE**2)
						count += 1
		if count == 0:
			score = 0
		else:
			score = float(curve_score - surface_score) / count
		skeletonScores.append( score )
		
		#Normalization
		min_score = min(skeletonScores)
		max_score = max(skeletonScores)
		for i in range(len(skeletonScores)):
			if skeletonScores[i] > 0:
				skeletonScores[i] /= float(max_score)
			else:
				skeletonScores[i] /= -1.0*min_score
	
	return skeletonScores	

def write_composite_scores(target_volume_name, apix, patoms, atomNumber, coeffArray, coeffwt, avghhvalue, maxValue, skeletonArray, skeletonwt, pseudoatomArray, atomwt):
	index4=0
	scorefile="score-%s.pdb"%(target_volume_name)
	chain="A"
	outscore=open(scorefile,"w")
	
	for atom in atomNumber:
		score=0
		
		print "Atom Number: %s   "%(atomNumber[index4])
		
		print "Coordinates: (%.2f,%.2f,%.2f)"%(patoms[index4][0]*apix,patoms[index4][1]*apix,patoms[index4][2]*apix) #Debugging: by Ross
		
		print "	  Correlation:	   %s"%(coeffArray[index4]),
		if coeffArray[index4] >= (0.9*avghhvalue):
			tmpscore=(coeffArray[index4]/maxValue)
			score=score+tmpscore*coeffwt
			print " (+%f)"%(tmpscore)
		else:
			tmpscore=(coeffArray[index4]-avghhvalue)/avghhvalue
			score=score+tmpscore*coeffwt
			print " (%f)"%(tmpscore)	
		
		print "	  Skeleton:  %s"%(skeletonArray[index4])
		score=score+skeletonArray[index4]*skeletonwt
	
		
		print "	  Pseudoatom geometry:	  %s"%(pseudoatomArray[index4])
		score=score+pseudoatomArray[index4]*atomwt
	
		
		print "	Score: %f"%(score) 
		line="ATOM  %5d  C   GLY %s%4d	%8.3f%8.3f%8.3f  1.00%6.2f0	  S_00  0 \n"%(atom,chain,atom,patoms[index4][0]*apix,patoms[index4][1]*apix,patoms[index4][2]*apix,score)
		outscore.write(line)
		#outscore.write(line[index4][:60]+"%6.2f"%(score)+line[index4][66:])
		index4=index4+1
	outscore.close()


# If executed as a program
if __name__ == '__main__':
	main() 
