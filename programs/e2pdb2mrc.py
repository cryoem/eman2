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

# e2pdb2mrc.py  07/23/3004  Steven Ludtke
# This program will generate an electron density map from a PDB file. Unlike the earlier versions
# of this program, operations like applying symmetry or centering are now offloaded onto
# e2procpdb.py. Note that atomic form factors are not included in this program. It is designed 
# for intermediate resolutions (~4 A and higher). Each atom is represented by a Gaussian with
# the approrpiate number of electrons.

# PDB sample line
#           1         2         3         4         5         6         7
# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
# ATOM      1  N   ASP L   1      57.169  62.936   7.756  1.00 30.65      1ACY 129
# HEADER    COMPLEX(ANTIBODY/HIV-1 FRAGMENT)        10-FEB-94   1ACY      1ACY   2


from EMAN2 import *
from math import *
import os
import sys
import traceback

# HO is a hydrogen attached to an oxygen. 'W' is water (infrequently found)
atomdefs={'H':(1.0,1.00794),'HO':(1.0,1.00794),'C':(6.0,12.0107),'A':(7.0,14.00674),'N':(7.0,14.00674),'O':(8.0,15.9994),'P':(15.0,30.973761),'K':(19.0,39.0983),
	'S':(16.0,32.066),'W':(18.0,1.00794*2.0+15.9994),'AU':(79.0,196.96655) }

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] input.pdb output.mrc
	
	Converts a pdb file into an electron density map. 0,0,0 in PDB space will 
	map to the center of the volume. Use e2procpdb.py to adjust coordinates,
	apply symmetry, etc. Resolution is equivalent to standard cryoEM definition, 
	using 1/2 width of Gaussian in Fourier space."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--apix", "-A", type=float, help="A/voxel", default=1.0)
	parser.add_argument("--res", "-R", type=float, help="Resolution in A, equivalent to Gaussian lowpass with 1/e width at 1/res",default=2.8)
	parser.add_argument("--box", "-B", type=str, help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	parser.add_argument("--het", action="store_true", help="Include HET atoms in the map", default=False)
	parser.add_argument("--center", action="store_true", help="Move the atomic center to the center of the box", default=False)
	parser.add_argument("--chains",type=str,help="String list of chain identifiers to include, eg 'ABEFG'")
	parser.add_argument("--info", action="store_true", help="If this is specified, information on the PDB file is displayed, no conversion is performed.",default=False)
	#parser.add_argument("--full", action="store_true", help="Apply non-crystallographic symmetry (MTRIXs or BIOMTs) to get full structure.",default=False)
	parser.add_argument("--quiet",action="store_true",default=False,help="Verbose is the default")
	parser.add_argument("--model", type=int,default=None, help="Extract only a single numbered model from a multi-model PDB")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--quick", action="store_true", help="Use a slight approximation to the Gaussian during insertion. Does not support B-factors.",default=False)
	parser.add_argument("--addpdbbfactor", action="store_true", help="Use the bfactor/temperature factor as the atom blurring radius, equivalent to Gaussian lowpass with 1/e width at 1/bfactor",default=False)
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	try: chains=options.chains
	except: chains=None
	try: box=options.box
	except: box=None	

	logger = E2init(sys.argv,options.ppid)
		
	if not options.quick:
		
		try: boxs=options.box.split(',')
		except:
			pass 
		try: 
			boxs=[int(i) for i in boxs]
			boxsizeinput=max(boxs)
		except:
			pass				
		
		if options.addpdbbfactor==False:
			addpdbbfactor=-1
		else:
			addpdbbfactor=1
		
		try : infile=open(args[0],"r")
		except : raise IOError("%s is an invalid file name" %file_name)
	
		if options.apix > options.res:
			print("The specified --apix must be <= the specified --res. Generally, res should be >= 2*apix.")
			sys.exit(1)

		if options.solvate and (options.apix > 1.0 or options.res > 2.5):
			print("For solvation to work properly, --apix must be <=1.0 and --res must be <=2.5.")
			sys.exit(1)

		aavg=[0,0,0]	# calculate atomic center
		amin=[1.0e20,1.0e20,1.0e20]		# coords of lower left front corner of bounding box
		amax=[-1.0e20,-1.0e20,-1.0e20]	# coords
		natm=0
		atoms=[]		# we store a list of atoms to process to avoid multiple passes
		nelec=0
		mass=0

	# parse the pdb file and pull out relevant atoms
		stm=False
		
		ihelix = 0
		isheet = 0
		ires = 0
		lastres = ""

		#if options.full:
		#	tfs = []

		for line in infile:
			if options.model!=None:
				if line[:5]=="MODEL":
					if int(line.split()[1])==options.model: stm=True
				if stm and line[:6]=="ENDMDL" : break
				if not stm: continue

		#	if options.full:
		#		if line[:5] == "MTRIX" or line[:18] == "REMARK 350   BIOMT":
		#			_,_,_,rx,ry,rz,trans = line.strip().split()
		#			tf = Transform()
		#			tf.set_matrix([rx,0,0,0,0,ry,0,0,0,0,rz,0])
		#			tf.set_trans(trans)
		#			tfs.append(tf)

			if line[:5] =="HELIX":
				ihelix += int(line[72:77])
				
			elif line[:5]=="SHEET":
				isheet += abs(int(line[22:26])-int(line[33:37]))+1

			elif (line[:4]=='ATOM' or (line[:6]=='HETATM' and options.het)) :
				
				if lastres != line[17:21]:
					ires += 1
					lres = line[17:21]

				if chains and not (line[21] in chains) : continue

				try:
					aseq=int(line[6:11].strip())
				except:
					aseq=-1

				try:
					resol=int(line[22:26].strip())
				except:
					resol=10

				try:
					a=line[12:14].strip().translate(None,"0123456789")
	
					x=float(line[30:38])
					y=float(line[38:46])
					z=float(line[46:54])
				except:
					print "PDB Parse error:\n%s\n'%s','%s','%s'  '%s','%s','%s'\n"%(
						line,line[12:14],line[6:11],line[22:26],line[30:38],line[38:46],line[46:54])
					print a,aseq,resol,x,y,z

				atoms.append((a,x,y,z))
						
				aavg[0]+=x
				aavg[1]+=y
				aavg[2]+=z
				natm+=1
				
				amin[0]=min(x,amin[0])
				amin[1]=min(y,amin[1])
				amin[2]=min(z,amin[2])
				amax[0]=max(x,amax[0])
				amax[1]=max(y,amax[1])
				amax[2]=max(z,amax[2])
	
				try:
					nelec+=atomdefs[a[0].translate(None,"0123456789").upper()][0]
					mass+=atomdefs[a[0].translate(None,"0123456789").upper()][1]
					#nelec+=atomdefs[a.upper()][0]
					#mass+=atomdefs[a.upper()][1]
				except:
					print("Unknown atom %s ignored at %d"%(a,aseq))
		
		infile.close()
		
		print "%d atoms used with a total charge of %d e- and a mass of %d kDa"%(natm,nelec,mass/1000)
		print "Atomic center at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"%(aavg[0]/natm,aavg[1]/natm,aavg[2]/natm)
		print "Bounding box: x: %7.2f - %7.2f"%(amin[0],amax[0])
		print "              y: %7.2f - %7.2f"%(amin[1],amax[1])
		print "              z: %7.2f - %7.2f"%(amin[2],amax[2])
		print("Helix (%): {}".format(100*(ihelix/ires)))
		print("Sheet (%): {}".format(100*(isheet/ires)))
	
		if options.info: sys.exit(1)

		#print "boxsize=", boxsize
		try:
			boxsize=boxsizeinput
		except:
			boxsize=max([(amax[0]-amin[0]),(amax[1]-amin[1]),(amax[2]-amin[2])])
			boxsize=int(1.9*boxsize+(2.0*options.res/options.apix))

		pa=PointArray()
		pa.read_from_pdb(args[0])
		if options.center:
			pa.center_to_zero()
		#if options.full:
		#	for tf in tfs:
		#		new_pa = pa.rotate(tf) + PointArray(tf.get_trans())
		#		out = pa.pdb2mrc_by_summation(boxsize,options.apix,options.res,addpdbbfactor)
		#else:
		out=pa.pdb2mrc_by_summation(boxsize,options.apix,options.res,addpdbbfactor)
		out.write_image(args[1])
		
		print '\r   %d\nConversion complete'%len(atoms)

	else:
		if options.addpdbbfactor : print "WARNING: B-factors not supported in quick mode"
		#if options.full : print "WARNING: --full not supported in quick mode"
		outmap = pdb_2_mrc(args[0],options.apix,options.res,options.het,box,chains,options.model,options.center,options.quiet)
		outmap.write_image(args[1])

	E2end(logger)
						
# this function originally added so that it could be accessed independently (for Junjie Zhang by David Woolford)
def pdb_2_mrc(file_name,apix=1.0,res=2.8,het=False,box=None,chains=None,model=None,center=False,quiet=False):
	'''
	file_name is the name of a pdb file
	apix is the angstrom per pixel
	res is requested resolution, quivalent to Gaussian lowpass with 1/e width at 1/res
	het is a flag inidicating whether HET atoms should be included in the map
	box is the boxsize, can be a single int (e.g. 128), a tuple (e.g. [128,64,54]), or a string (e.g. "128" or "128,64,57")
	chains is a string list of chain identifiers, eg 'ABEFG'
	quiet can be used to turn of helpful print outs
	'''
	
	try : infile=open(file_name,"r")
	except : raise IOError("%s is an invalid file name" %file_name)
	
	
	if res<=apix : print "Warning: res<=apix. Generally res should be 2x apix or more"
	
	aavg=[0,0,0]	# calculate atomic center
	amin=[1.0e20,1.0e20,1.0e20]		# coords of lower left front corner of bounding box
	amax=[-1.0e20,-1.0e20,-1.0e20]	# coords
	natm=0
	atoms=[]		# we store a list of atoms to process to avoid multiple passes
	nelec=0
	mass=0

	# parse the pdb file and pull out relevant atoms
	stm=False
	for line in infile:
		if model!=None:
			if line[:5]=="MODEL":
				if int(line.split()[1])==model: stm=True
			if stm and line[:6]=="ENDMDL" : break
			if not stm: continue

		if (line[:4]=='ATOM' or (line[:6]=='HETATM' and het)) :
			if chains and not (line[21] in chains) : continue
			
			try:
				aseq=int(line[6:11].strip())
			except:
				aseq=-1

			try:
				resol=int(line[22:26].strip())
			except:
				resol=10

			try:
				a=line[12:14].strip().translate(None,"0123456789'")
	
				x=float(line[30:38])
				y=float(line[38:46])
				z=float(line[46:54])
			except:
				print "PDB Parse error:\n%s\n'%s','%s','%s'  '%s','%s','%s'\n"%(
					line,line[12:14],line[6:11],line[22:26],line[30:38],line[38:46],line[46:54])
				print a,aseq,resol,x,y,z

			atoms.append((a,x,y,z))
						
			aavg[0]+=x
			aavg[1]+=y
			aavg[2]+=z
			natm+=1
			
			amin[0]=min(x,amin[0])
			amin[1]=min(y,amin[1])
			amin[2]=min(z,amin[2])
			amax[0]=max(x,amax[0])
			amax[1]=max(y,amax[1])
			amax[2]=max(z,amax[2])

			try:
				nelec+=atomdefs[a.upper()][0]
				mass+=atomdefs[a.upper()][1]
			except:
				print("Unknown atom %s ignored at %d"%(a,aseq))
							
	infile.close()
	
	if not quiet:
		print "%d atoms used with a total charge of %d e- and a mass of %d kDa"%(natm,nelec,mass/1000)
		print "atomic center at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"%(aavg[0]/natm,aavg[1]/natm,aavg[2]/natm)
		print "Bounding box: x: %7.2f - %7.2f"%(amin[0],amax[0])
		print "              y: %7.2f - %7.2f"%(amin[1],amax[1])
		print "              z: %7.2f - %7.2f"%(amin[2],amax[2])
	
	# precalculate a prototypical Gaussian to resample
	# 64^3 box with a real-space 1/2 width of 12 pixels
	gaus=EMData()
	gaus.set_size(64,64,64)
	gaus.to_one()
	
	gaus.process_inplace("mask.gaussian",{"outer_radius":12.0})

	# find the output box size, either user specified or from bounding box
	outbox=[0,0,0]
	try:
		# made
		if isinstance(box,int):
			outbox[0]=outbox[1]=outbox[2]=box
		elif isinstance(box,list):
			outbox[0]=box[0]
			outbox[1]=box[1]
			outbox[2]=box[2]
		else:
			spl=box.split(',')
			if len(spl)==1 : outbox[0]=outbox[1]=outbox[2]=int(spl[0])
			else :
				outbox[0]=int(spl[0])
				outbox[1]=int(spl[1])
				outbox[2]=int(spl[2])
	except:
		pad=int(2.0*res/apix)
		outbox[0]=int((amax[0]-amin[0])/apix)+pad
		outbox[1]=int((amax[1]-amin[1])/apix)+pad
		outbox[2]=int((amax[2]-amin[2])/apix)+pad
		outbox[0]+=outbox[0]%2
		outbox[1]+=outbox[1]%2
		outbox[2]+=outbox[2]%2
		
	if not quiet: print "Box size: %d x %d x %d"%(outbox[0],outbox[1],outbox[2])
	
	# initialize the final output volume
	outmap=EMData()
	outmap.set_size(outbox[0],outbox[1],outbox[2])
	outmap.to_zero()	
	for i in range(len(aavg)): aavg[i] = aavg[i]/float(natm)	
	# fill in the atom gaussians
	xt = outbox[0]/2 - (amax[0]-amin[0])/(2*apix)
	yt = outbox[1]/2 - (amax[1]-amin[1])/(2*apix)
	zt = outbox[2]/2 - (amax[2]-amin[2])/(2*apix)
	for i,a in enumerate(atoms):
		if not quiet and i%1000==0 : 
			print '\r   %d'%i,
			sys.stdout.flush()
		try:
			# This insertion strategy ensures the output is centered.
			elec=atomdefs[a[0].translate(None,"0123456789").upper()][0]
# This was producing different results than the "quick" mode, and did not match the statement printed above!!!
#			outmap.insert_scaled_sum(gaus,(a[1]/apix+xt-amin[0]/apix,a[2]/apix+yt-amin[1]/apix,a[3]/apix+zt-amin[2]/apix),res/(pi*12.0*apix),elec)
			if center: outmap.insert_scaled_sum(gaus,((a[1]-aavg[0])/apix+outbox[0]/2,(a[2]-aavg[1])/apix+outbox[1]/2,(a[3]-aavg[2])/apix+outbox[2]/2),res/(pi*12.0*apix),elec)
			else: outmap.insert_scaled_sum(gaus,(a[1]/apix+outbox[0]/2,a[2]/apix+outbox[1]/2,a[3]/apix+outbox[2]/2),res/(pi*12.0*apix),elec)
		except: print "Skipping %d '%s'"%(i,a[0])		
	if not quiet: print '\r   %d\nConversion complete'%len(atoms)		
	outmap.set_attr("apix_x",apix)
	outmap.set_attr("apix_y",apix)
	outmap.set_attr("apix_z",apix)	
	outmap.set_attr("origin_x",-xt*apix+amin[0])	
	outmap.set_attr("origin_y",-yt*apix+amin[1])	
	outmap.set_attr("origin_z",-zt*apix+amin[2])
	return outmap
	
if __name__ == "__main__":
    main()
