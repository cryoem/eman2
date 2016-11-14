# simtomo.py		Steven Ludtke	2/2011
# This program will generate a set of 'tomographic single particles' in random orientation
# with tomographic sampling. Generates a set of projections for each orientation, adds flatband noise, then reconstructs
# from the projections.

from EMAN2 import *

from sys import argv
from os import system
import random

optcl=EMData(argv[1],0)
log=file("tomo.ort","w")

# 8 particles
for m in range(32):
	
	xf=Transform({"type":"eman","az":random.uniform(0,360.0),"phi":random.uniform(0,360.0),"alt":random.uniform(0,360.0)})
	print "Model ",m,xf
	log.write("%d\t%s\n"%(m,str(xf)))
	ptcl=optcl.process("xform",{"transform":xf})

	# altitude range and step
	i=0
	for alt in range(-54,54,2):
		proj=ptcl.project("standard",Transform({"type":"eman","alt":alt}))
		proj.process_inplace("normalize.edgemean")
		
		noise=test_image(1,size=(ptcl["nx"],ptcl["ny"]))
		noise2=noise.process("filter.lowpass.gauss",{"cutoff_abs":.25})
		noise.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.75})
		noise=noise*3+noise2*3

		proj.add(noise)

		proj.write_image("proj.hdf",i)
		i+=1

	system("e2make3d.py --input=proj.hdf --output=recon.hdf --pad=%d --keep=1 --no_wt --iter=2"%proj["nx"])

	recon=EMData("recon.hdf",0)

	recon.write_image("tomo.hdf",m)
	
