#!/bin/env python

# PDB sample line
#           1         2         3         4         5         6         7
# 01234567890123456789012345678901234567890123456789012345678901234567890123456789
# ATOM      1  N   ASP L   1      57.169  62.936   7.756  1.00 30.65      1ACY 129
# HEADER    COMPLEX(ANTIBODY/HIV-1 FRAGMENT)        10-FEB-94   1ACY      1ACY   2


from EMAN2 import *
from optparse import OptionParser

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

	parser.add_option("--apix", "-A", type="float", help="Angstroms/voxel")
	parser.add_option("--res", "-R", type="float", help="Resolution, equivalent to Gaussian lowpass with 1/e width at 1/res")
	parser.add_option("--box", "-B", type="string", help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	parser.add_option("--het", action="store_true", help="Include HET atoms in the map")
	parser.add_option("--chains",type="string",help="String list of chain identifiers to include, eg 'ABEFG'")
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	try: chains=options.chains
	except: chains=None
	
	try : infile=open(args[0],"r")
	except : parser.error("Cannot open input file")
	
	# find the output box size
	try:
		spl=options.box.split(',')
		if len(spl==1) : box[0]=box[1]=box[2]=int(spl[0])
		else :
			box[0]=int(spl[0])
			box[1]=int(spl[1])
			box[2]=int(spl[2])
	except:
		box[0]=box[1]=box[2]=0
				
	aavg=[0,0,0]	# calculate atomic center
	amin=[1.0e20,1.0e20,1.0e20]		# coords of lower left front corner of bounding box
	amax=[-1.0e20,-1.0e20,-1.0e20]	# coords
	natm=0
	atoms=[]		# we store a list of atoms to process to avoid multiple passes
	nelec=0
	mass=0

	# parse the pdb file and pull out relevant atoms
	for line in infile:
		if (line[:4]=='ATOM' or line[:6]=='HETATM') :
			if chains and not (line[21] in chains) : continue
			
			a=line[12:14].strip()
			aseq=int(line[6:10].strip())
			res=int(line[22:26].strip())
	
			x=float(line[30:38])
			y=float(line[38:46])
			z=float(line[46:54])
			
			atoms.append((a,x,y,z))
						
			aavg[0]+=x
			aavg[1]+=y
			aavg[2]+=z
			natm+=1
			
			amin[0]=min(x,amin[0])
			amin[1]=min(y,amin[1])
			amin[2]=min(z,amin[2])
			amax[0]=max(x,amax[0])
			amax[1]=max[y,amax[1])
			amax[2]=max(z,amax[2])

			try:
				nelec+=atomdefs[a.upper()][0]
				mass+=atomdefs[a.upper()][1]
			except:
				print("Unknown atom %s ignored at %d"%(a,aseq))
							
	infile.close()
	
	print "%d atoms used with a total charge of %d e- and a mass of %d kDa"%(natm,nelec,mass/1000)
	print "atomic center at %1.1f,%1.1f,%1.1f (center of volume at 0,0,0)"%(aavg[0]/natm,aavg[1]/natm,aavg[2]/natm)
	
	gaus=EMData()
	gaus.set_size(64,64,64)
	gaus.to_one()
	gaus.filter("GaussMask",{"outer_radius":20.0})
	
if __name__ == "__main__":
    main()
