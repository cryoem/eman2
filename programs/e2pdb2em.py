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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

#  Implemented by PAP 02/27/09.  Based on Agarwal Acta Cryst. 1978, A34, 791.

# PDB sample line
#           1         2         3         4         5         6         7
# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
# ATOM      1  N   ASP L   1      57.169  62.936   7.756  1.00 30.65      1ACY 129
# HEADER    COMPLEX(ANTIBODY/HIV-1 FRAGMENT)        10-FEB-94   1ACY      1ACY   2



from EMAN2 import *
from sparx import *
from math import *
from global_def import *
import os
import sys

atomdefs={'H':(1.0,1.00794),'C':(6.0,12.0107),'A':(7.0,14.00674),'N':(7.0,14.00674),'O':(8.0,15.9994),'P':(15.0,30.973761),
	'S':(16.0,32.066),'W':(18.0,1.00794*2.0+15.9994),'AU':(79.0,196.96655) }

def create_Fgaus(nx,ny,nz,pixel_size,C1,B1,C2,B2):
	# precalculate a Gaussian function in Fourier space (data from Table 1)
	gaus=EMData(nx,ny,nz, False)  # Fourier object!
	nxp = nx//2
	nyp = ny//2
	nzp = nz//2
	dx2 = (1.0/float(nxp))**2
	dy2 = (1.0/nyp)**2
	dz2 = (1.0/nzp)**2
	for i in xrange((nx+2)//2):
		argx = float(i*i)*dx2
		for j in xrange(ny):
			jy=j
			if(jy>nyp): jy=jy-nyp
			argy = argx + float(jy*jy)*dy2
			for k in xrange(nz):
				jz = k
				if (jz>nzp): jz=jz-nzp
				argz = (float(jz*jz)*dz2 + argy)/pixel_size/pixel_size
				gaus[i,j,k] = C1*exp(-B1*argz/4) + C2*exp(-B2*argz/4)   #Eq.40
	return gaus

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] input.pdb output.hdf
	
	Converts a pdb file into an electron density map. 0,0,0 in PDB space will 
	map to the center of the volume. Each atom is represented by a Gaussian with a width 
	defined by the specified resolution, and an amplitude equal to the atomic number.

	While this is only an approximation, for typical cryo-EM work it is more than adequate.
	Atomic form factors will have little impact until close to atomic resolution at
	which point molecular orbitals really should also be considered."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--apix", "-A", type=float, help="A/voxel", default=1.0)
	parser.add_argument("--box", "-B", type=str, help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	parser.add_argument("--het", action="store_true", help="Include HET atoms in the map", default=False)
	#parser.add_argument("--chains",type=str,help="String list of chain identifiers to include, eg 'ABEFG'")
	parser.add_argument("--center",  type=str,  default="a", help="center: c - coordinates, a - center of gravity, n - no" )
	parser.add_argument("--O",  action="store_true", default=False, help="use O system of coordinates")
	parser.add_argument("--quiet",action="store_true",default=False,help="Verbose is the default")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	#try: chains=options.chains
	#except: 
	chains=None
	
	try : infile=open(args[0],"r")
	except : parser.error("Cannot open input file")
	
	aavg=[0,0,0]	# calculate atomic center
	amin=[1.0e20,1.0e20,1.0e20]		# coords of lower left front corner of bounding box
	amax=[-1.0e20,-1.0e20,-1.0e20]	# coords
	natm=0
	atoms=[]		# we store a list of atoms to process to avoid multiple passes
	nelec=0
	mass=0

	encountered_atoms = []
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
				if( not(a.upper() in encountered_atoms)):  encountered_atoms.append(a.upper())
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
			
			amin[0]=min(x,amin[0])
			amin[1]=min(y,amin[1])
			amin[2]=min(z,amin[2])
			amax[0]=max(x,amax[0])
			amax[1]=max(y,amax[1])
			amax[2]=max(z,amax[2])

							
	infile.close()
	
	if not options.quiet:
		print "%d atoms used with a total charge of %d e- and a mass of %d kDa"%(natm,nelec,mass/1000)
		if(options.center == "a"):
			print "center of gravity at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"%(aavg[0]/mass,aavg[1]/mass,aavg[2]/mass)
			for i in xrange( len(atoms) ) :
				atoms[i][1] -= aavg[0]/mass
				atoms[i][2] -= aavg[1]/mass
				atoms[i][3] -= aavg[2]/mass
			for i in xrange(3): amin[i] -= aavg[i]/mass
		else:
			print "atomic center at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"%(aavg[0]/natm,aavg[1]/natm,aavg[2]/natm)
			if(options.center == "c"):
				for i in xrange( len(atoms) ) :
					atoms[i][1] -= aavg[0]/natm
					atoms[i][2] -= aavg[1]/natm
					atoms[i][3] -= aavg[2]/natm
				for i in xrange(3): amin[i] -= aavg[i]/natm
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
		box[0]=int(2*max(fabs(amin[0]),fabs(amax[0]))/options.apix)
		box[1]=int(2*max(fabs(amin[1]),fabs(amax[1]))/options.apix)
		box[2]=int(2*max(fabs(amin[2]),fabs(amax[2]))/options.apix)
		box[0]+=box[0]%2
		box[1]+=box[1]%2
		box[2]+=box[2]%2

	box[0]=box[1]=box[2]=max(box[0], box[1], box[2])

	if not options.quiet: print "Box size: %d x %d x %d"%(box[0],box[1],box[2])

	# Generate atom prototypes (constants from Table 1).
	for atom in encountered_atoms:
		if(atom == "H"):   gaussH = create_Fgaus(box[0],box[1],box[2], options.apix, 0.486, 34.1284, 0.5098, 8.8996)
		elif(atom == "C"):  gaussN = create_Fgaus(box[0],box[1],box[2], options.apix, 3.0102,29.9132, 2.9705, 2.8724)
		elif(atom == "N"):  gaussO = create_Fgaus(box[0],box[1],box[2], options.apix, 3.0492, 25.0383, 3.9432, 3.4059)
		elif(atom == "O"):  gaussC = create_Fgaus(box[0],box[1],box[2], options.apix, 3.2942,20.0401,4.6968,3.1184)
		else:              gaussS = create_Fgaus(box[0],box[1],box[2], options.apix, 5.6604, 33.04, 10.314, 1.816)

	del encountered_atoms

	# initialize the final output volume
	outmap=EMData(box[0],box[1],box[2], False)

	# fill in the atoms
	for i in xrange(len(atoms)):
		#print "Adding %d '%s'"%(i,atoms[i][0])
		if not options.quiet and i%1000==0 :
			print '\r   %d'%i,
			sys.stdout.flush()
		try:
			atom = atoms[i][0].upper()
			if(atom == "H"):   gaus = gaussH
			elif(atom == "O"):  gaus = gaussO
			elif(atom == "C"):  gaus = gaussC
			elif(atom == "N"):  gaus = gaussN
			else:              gaus = gaussS
			Util.add_img(outmap, fshift(gaus, atoms[i][1]/options.apix, atoms[i][2]/options.apix, atoms[i][3]/options.apix) )
			#elec = atomdefs[atoms[i][0].upper()][0]
			#Util.add_img(outmap, fshift(gaus,atoms[i][1]/options.apix,atoms[i][2]/options.apix,atoms[i][3]/options.apix)*elec)
			#outmap.insert_scaled_sum(gaus,(atoms[i][1]/options.apix+box[0]/2,atoms[i][2]/options.apix+box[1]/2,atoms[i][3]/options.apix+box[2]/2),options.res/(pi*12.0*options.apix),elec)
		except: print "Skipping %d '%s'"%(i,atoms[i][0])
		
	if not options.quiet: print '\r   %d\nConversion complete'%len(atoms)
	info(outmap)
	outmap = fft(fshift(outmap,box[0]//2,box[1]//2,box[2]//2))
	info(outmap)
	(filename_path, filextension) = os.path.splitext(args[1])
	if filextension == ".hdf" :
		outmap.set_attr("apix_x",options.apix)
		outmap.set_attr("apix_y",options.apix)
		outmap.set_attr("apix_z",options.apix)
		outmap.set_attr("pixel_size",options.apix)
		outmap.write_image(args[1],0, EMUtil.ImageType.IMAGE_HDF)
	elif filextension == ".spi": outmap.write_image(args[1],0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER)
	else:   ERROR("unknown image type","e2pdb2em",1)
				
if __name__ == "__main__":
    main()

