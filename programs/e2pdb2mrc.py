#!/bin/env python
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
from optparse import OptionParser
from math import *
import os
import sys

atomdefs={'H':(1.0,1.00794),'C':(6.0,12.0107),'A':(7.0,14.00674),'N':(7.0,14.00674),'O':(8.0,15.9994),'P':(15.0,30.973761),
	'S':(16.0,32.066),'W':(18.0,1.00794*2.0+15.9994),'AU':(79.0,196.96655) }

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: %prog [options] input.pdb output.mrc
	
Converts a pdb file into an electron density map. 0,0,0 in PDB space will 
map to the center of the volume. Use e2procpdb.py to adjust coordinates,
apply symmetry, etc. Resolution is equivalent to standard cryoEM definition, 
using 1/2 width of Gaussian in Fourier space."""

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
	
	if options.res<=options.apix : print "Warning: res<=apix. Generally res should be 2x apix or more"
	
	aavg=[0,0,0]	# calculate atomic center
	amin=[1.0e20,1.0e20,1.0e20]		# coords of lower left front corner of bounding box
	amax=[-1.0e20,-1.0e20,-1.0e20]	# coords
	natm=0
	atoms=[]		# we store a list of atoms to process to avoid multiple passes
	nelec=0
	mass=0

	# parse the pdb file and pull out relevant atoms
	for line in infile:
		if (line[:4]=='ATOM' or (line[:6]=='HETATM' and options.het)) :
			if chains and not (line[21] in chains) : continue
			
			try:
				a=line[12:14].strip()
				aseq=int(line[6:11].strip())
				res=int(line[22:26].strip())
	
				x=float(line[30:38])
				y=float(line[38:46])
				z=float(line[46:54])
			except:
				print "PDB Parse error:\n%s\n'%s','%s','%s'  '%s','%s','%s'\n"%(
					line,line[12:14],line[6:11],line[22:26],line[30:38],line[38:46],line[46:54])
				print a,aseq,res,x,y,z

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
	
	if not options.quiet:
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
	gaus.filter("eman1.mask.gaussian",{"outer_radius":12.0})

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
		pad=int(2.0*options.res/options.apix)
		box[0]=int(2*max(fabs(amin[0]),fabs(amax[0]))/options.apix)+pad
		box[1]=int(2*max(fabs(amin[1]),fabs(amax[1]))/options.apix)+pad
		box[2]=int(2*max(fabs(amin[2]),fabs(amax[2]))/options.apix)+pad
		box[0]+=box[0]%2
		box[1]+=box[1]%2
		box[2]+=box[2]%2
		
	if not options.quiet: print "Box size: %d x %d x %d"%(box[0],box[1],box[2])
	
	# initialize the final output volume
	outmap=EMData()
	outmap.set_size(box[0],box[1],box[2])
	outmap.to_zero()
	
	# fill in the atom gaussians
	for i,a in enumerate(atoms):
		if not options.quiet and i%1000==0 : 
			print '\r   %d'%i,
			sys.stdout.flush()
		elec=atomdefs[a[0].upper()][0]
		outmap.insert_scaled_sum(gaus,(a[1]/options.apix+box[0]/2,a[2]/options.apix+box[1]/2,a[3]/options.apix+box[2]/2),options.res/(pi*12.0*options.apix),elec)
		
	if not options.quiet: print '\r   %d\nConversion complete'%len(atoms)
	outmap.set_attr("apix_x",options.apix)
	outmap.set_attr("apix_y",options.apix)
	outmap.set_attr("apix_z",options.apix)
	outmap.write_image(args[1])
				
if __name__ == "__main__":
    main()
