#!/usr/bin/env python

#
# Author: Pawel Penczek, 4/4/2007 (Pawel.A.Penczek@uth.tmc.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

#  Implemented by PAP 03/02/09.  Based on Agarwal Acta Cryst. 1978, A34, 791.

# PDB sample line
#           1         2         3         4         5         6         7
# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
# ATOM      1  N   ASP L   1      57.169  62.936   7.756  1.00 30.65      1ACY 129
# HEADER    COMPLEX(ANTIBODY/HIV-1 FRAGMENT)        10-FEB-94   1ACY      1ACY   2



from EMAN2 import *
from sparx import *
from optparse import OptionParser
from math import *
from global_def import *
import global_def
import os
import sys

atomdefs={'H':(1.0,1.00794),'C':(6.0,12.0107),'A':(7.0,14.00674),'N':(7.0,14.00674),'O':(8.0,15.9994),'P':(15.0,30.973761),
	'S':(16.0,32.066),'W':(18.0,1.00794*2.0+15.9994),'AU':(79.0,196.96655) }

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] input.pdb output.hdf
	
Converts a pdb file into an electron density map. 0,0,0 in PDB space will 
map to the center of the volume."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--apix", "-A", type="float",        help="Angstrom/voxel", default=1.0)
	parser.add_option("--box",  "-B", type="string",       help="Box size in pixels, <xyz> or <x,y,z>")
	parser.add_option("--het",        action="store_true", help="Include HET atoms in the map", default=False)
	#parser.add_option("--chains",type="string",help="String list of chain identifiers to include, eg 'ABEFG'")
	parser.add_option("--center",     type="string",       default="n", help="center: c - coordinates; a - center of gravity; <x,y,z> - a vector (in Angstrom) to substract from all coordinates; default: n - no" )
	parser.add_option("--O",          action="store_true", default=False, help="use O system of coordinates")
	parser.add_option("--quiet",      action="store_true", default=False, help="Verbose is the default")
	parser.add_option("--tr0",        type="string",       default="none", help="Filename of initial 3x4 transformation matrix")

	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	#try: chains=options.chains
	#except: 
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	chains=None
	
	try : infile=open(args[0],"r")
	except : parser.error("Cannot open input file")
	
	aavg=[0,0,0]	# calculate atomic center
	natm=0
	atoms=[]		# we store a list of atoms to process to avoid multiple passes
	nelec=0
	mass=0

	# read in initial-transformation file:
	if(options.tr0 != "none"):
		cols = read_text_file(options.tr0,-1)
		txlist=[]
		for i in range(3):
			txlist.append(cols[0][i])
			txlist.append(cols[1][i])
			txlist.append(cols[2][i])
			txlist.append(cols[3][i])
		tr0 = Transform(txlist)


	# parse the pdb file and pull out relevant atoms
	for line in infile:
		if (line[:4]=='ATOM' or (line[:6]=='HETATM' and options.het)) :
			if chains and not (line[21] in chains) : continue
			
			try:
				a=line[12:14].strip()
				aseq=int(line[6:11].strip())
				res=int(line[22:26].strip())
				if(options.O):
					x=float(line[38:46])
					y=float(line[30:38])
					z=-float(line[46:54])
				else:
					x=float(line[30:38])
					y=float(line[38:46])
					z=float(line[46:54])
			except:
				print "PDB Parse error:\n%s\n'%s','%s','%s'  '%s','%s','%s'\n"%(
					line,line[12:14],line[6:11],line[22:26],line[30:38],line[38:46],line[46:54])
				print a,aseq,res,x,y,z
			try:
				nelec += atomdefs[a.upper()][0]
				mass  += atomdefs[a.upper()][1]
			except:
				print("Unknown atom %s ignored at %d"%(a,aseq))

			atoms.append([a,x,y,z])

			if(options.center == "a"):
				aavg[0] += x*atomdefs[a.upper()][1]
				aavg[1] += y*atomdefs[a.upper()][1]
				aavg[2] += z*atomdefs[a.upper()][1]
			else:
				aavg[0] += x
				aavg[1] += y
				aavg[2] += z
			natm += 1
			
	infile.close()
	
	if not options.quiet:
		print "%d atoms used with a total charge of %d e- and a mass of %d kDa"%(natm,nelec,mass/1000)

	# center PDB according to option:
	if(options.center == "a"):
		if not options.quiet:
			print "center of gravity at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"%(aavg[0]/mass,aavg[1]/mass,aavg[2]/mass)
		for i in xrange( len(atoms) ) :
			atoms[i][1] -= aavg[0]/mass
			atoms[i][2] -= aavg[1]/mass
			atoms[i][3] -= aavg[2]/mass
	if(options.center == "c"):
		if not options.quiet:
			print "atomic center at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"%(aavg[0]/natm,aavg[1]/natm,aavg[2]/natm)
		for i in xrange( len(atoms) ) :
			atoms[i][1] -= aavg[0]/natm
			atoms[i][2] -= aavg[1]/natm
			atoms[i][3] -= aavg[2]/natm
	spl = options.center.split(',')
	if len(spl)==3:   # substract the given vector from all coordinates
		if not options.quiet:
			print "vector to substract: %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"%(float(spl[0]),float(spl[1]),float(spl[2]))
		for i in xrange( len(atoms) ) :
			atoms[i][1] -= float(spl[0])
			atoms[i][2] -= float(spl[1])
			atoms[i][3] -= float(spl[2])

	# apply given initial transformation (this used to be done before centering,
	# thereby loosing the translation. This is the right place to apply tr0):
	if(options.tr0 != "none"):
		if not options.quiet:
			print "Applying initial transformation to PDB coordinates... "
		for i in xrange(len(atoms)):
			atom_coords = Vec3f(atoms[i][1],atoms[i][2],atoms[i][3])
			new_atom_coords = tr0*atom_coords
			atoms[i][1] = new_atom_coords[0]
			atoms[i][2] = new_atom_coords[1]
			atoms[i][3] = new_atom_coords[2]
		if not options.quiet:
			print "done.\n"

	# bounding box:
	amin=[atoms[0][1],atoms[0][2],atoms[0][3]]
	amax=[atoms[0][1],atoms[0][2],atoms[0][3]]
	for i in xrange(len(atoms)):
		amin[0]=min(atoms[i][1],amin[0])
		amin[1]=min(atoms[i][2],amin[1])
		amin[2]=min(atoms[i][3],amin[2])
		amax[0]=max(atoms[i][1],amax[0])
		amax[1]=max(atoms[i][2],amax[1])
		amax[2]=max(atoms[i][3],amax[2])

	if not options.quiet:
		print "Bounding box: x: %7.2f - %7.2f"%(amin[0],amax[0])
		print "              y: %7.2f - %7.2f"%(amin[1],amax[1])
		print "              z: %7.2f - %7.2f"%(amin[2],amax[2])
	
	# find the output box size, either user specified or from bounding box
	box=[0,0,0]
	try:
		spl=options.box.split(',')
		if len(spl)==1 : box[0]=box[1]=box[2]=int(spl[0])
		else :
			box[0]=int(spl[0])
			box[1]=int(spl[1])
			box[2]=int(spl[2])
	except:
		for i in xrange(3):
			box[i]=int(2*max(fabs(amax[i]), fabs(amin[i]))/options.apix)
			#  Increase the box size by 1/4.
			box[i]+=box[i]//4
	# figure oversampled box size
	bigb = max(box[0],box[1],box[2])
	fcbig = 1
	"""
	while(bigb*fcbig <= 512):   # maximum boxsize = 512
		fcbig +=1
	fcbig -= 1
	if(fcbig == 0):
		print  "Pixel size too small resulting in a box size >512"
		sys.exit()
	"""
	if not options.quiet: print "Box size: %d x %d x %d"%(box[0],box[1],box[2]),",  oversampling ",fcbig

	# Calculate working dimensions
	pixelbig = options.apix/fcbig
	bigbox = []
	for i in xrange(3): bigbox.append(box[i]*fcbig)

	# initialize the final output volume
	outmap=EMData(bigbox[0],bigbox[1],bigbox[2], True)
	nc = []
	for i in xrange(3):  nc.append( bigbox[i]//2 )
	# fill in the atoms
	for i in xrange(len(atoms)):
		#print "Adding %d '%s'"%(i,atoms[i][0])
		if not options.quiet and i%1000==0 :
			print '\r   %d'%i,
			sys.stdout.flush()
		try:
			elec = atomdefs[atoms[i][0].upper()][0]
			#outmap[int(atoms[i][1]/pixelbig+bigbox[0]//2),int(atoms[i][2]/pixelbig+bigbox[1]//2),int(atoms[i][3]/pixelbig+bigbox[2]//2)] += elec
			for k in xrange(2):
			        pz = atoms[i][3]/pixelbig+nc[2]
			        dz = pz - int(pz)
			        uz = ((1-k) + (2*k-1)*dz)*elec
			        for l in xrange(2):
			        	py = atoms[i][2]/pixelbig+nc[1]
			        	dy = py - int(py)
			        	uy = ((1-l) + (2*l-1)*dy)*uz
			        	for m in xrange(2):
			        		px = atoms[i][1]/pixelbig+nc[0]
			        		dx = px - int(px)
			        		outmap[int(px)+m,int(py)+l,int(pz)+k] += ((1-m) + (2*m-1)*dx)*uy
		except: print "Skipping %d '%s'"%(i,atoms[i][0])
		
	if not options.quiet: print '\r   %d\nConversion complete.'%len(atoms)  #,"    Now shape atoms."
	"""
	fftip(outmap)
	# Atom in Fourier space has sigma = 0.41 [1/A]
	outmap = filt_gaussl(outmap,pixelbig*0.41)
	# prepare for Fourier truncation
	if(fcbig > 1):
		outmap = filt_tanl(outmap,0.5/fcbig, 0.02)
		outmap = outmap.FourTruncate(box[0],box[1],box[2], False)
	outmap = fft(outmap)
	"""
	(filename_path, filextension) = os.path.splitext(args[1])
	if filextension == ".hdf" :
		outmap.set_attr("apix_x",options.apix)
		outmap.set_attr("apix_y",options.apix)
		outmap.set_attr("apix_z",options.apix)
		outmap.set_attr("pixel_size",options.apix)
		outmap.write_image(args[1],0, EMUtil.ImageType.IMAGE_HDF)
	elif filextension == ".spi": outmap.write_image(args[1],0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
	else:   ERROR("unknown image type","sxpdb2em",1)
				
if __name__ == "__main__":
    main()

