#!/usr/bin/env python

#
# Author: Steven Ludtke, 08/26/14 (sludtke@bcm.edu)
# Copyright (c) 2014- Baylor College of Medicine
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

import os
import sys
import random
import time
import string
import math
from os import system
from os import unlink
from sys import argv
from EMAN2 import *
from EMAN2db import db_open_dict

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <refine_xx folder> <iteration number> <3D excluded region mask>

This program modifies raw particle data to subtract-away undesired portions of the density on a per-particle basis.
For example, one could take GroEL particles with D7 symmetry and extract the top and bottom heptamers as separate
particles for reconstruction. In this example there would be 2x as many output particles as input particles. The 
portion of the map identified by the excluded region mask will be subtracted from each particle as optimally as
possible, taking CTF, normalization and other factors into account. 

To use this program you should first run a complete refinement to the best possible resolution. Refinement parameters
will be examined automatically to extract the corresponding particles and projections.
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--sym", dest = "sym", default="c1",help = "Specify symmetry - choices are: c<n>, d<n>, tet, oct, icos.", guitype='strbox', row=10, col=1, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--debug",action="store_true",default=False,help="Enable debugging mode with verbose output and image display. Not suitable for real runs.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	#parser.add_argument("--ncls","-N",type=int,help="Number of classes to generate",default=-1)
	#parser.add_argument("--average","-A",action="store_true",help="Average the particles within each class",default=False)

	(options, args) = parser.parse_args()
#	if len(args)<4 : parser.error("Please specify <raw particle file> <class mx> <projections> <mask> <output file> ")

	# each of these classmx variables becomes a 2 element list with even and odd particles respectively
	try:
		pathmx="{}/classmx_{:02d}_even.hdf".format(args[0],args[1])
		classmx=[EMData(pathmx,0)]
		nptcl=[classmx["ny"]]
		cmxtx=[EMData(pathmx,2)]
		cmxty=[EMData(pathmx,3)]
		cmxalpha=[EMData(pathmx,4)]
		cmxmirror=[EMData(pathmx,5)]

		pathmx="{}/classmx_{:02d}_odd.hdf".format(args[0],args[1])
		classmx.append(EMData(pathmx,0))
		nptcl.append(classmx["ny"])
		cmxtx.append(EMData(pathmx,2))
		cmxty.append(EMData(pathmx,3))
		cmxalpha.append(EMData(pathmx,4))
		cmxmirror.append(EMData(pathmx,5))
	except:
		print "Error reading classification matrix. Must be full classification matrix with alignments"
		sys.exit(1)

	# this reads all of the EMData headers from the projections, should be same for even and odd
	pathprj="{}/projections_{:02d}_even.hdf".format(args[0],args[1])
	nref=EMUtil.get_image_count(pathprj)
	eulers=[EMData(pathprj,i,1)["xform.projection"] for i in range(nref)]

	# The 3D reference volume we are using for subtraction
	threed=EMData("{}/threed_{:02d}.hdf".format(args[0],args[1]),0)

	# The mask in "neutral" orientation
	mask=EMData(args[2],0)

	logid=E2init(sys.argv, options.ppid)

	# We will need masked versions of our volume for all symmetric orientations
	# generate the symmetry transforms here for convenience
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(options.sym)
	symxfs=[xf.get_sym(options.sym,i) for i in xrange(nsym)]

	# now we loop over each class, and assess the masked region for each particle in terms of
	# sigma of the image. Note that we don't have a list of which particle is in each class,
	# but rather a list of which class each particle is in, so we do this a bit inefficiently for now
	for i in xrange(nref):
		if options.verbose>1 : print "--- Class %d"%i

		# projection with the correct overall Euler angles for particles in one class
		proj=threed.project("standard",{"transform":eulers[i]})
		
		# We regenerate the masked volumes for each class to avoid using too much RAM
		projs=[]
		for xf in symxfs:
			maskx=mask.process("xform",{"transform":xf})
			masked=threed.copy()
			masked.mult(maskx)
			projs.append(masked.project("standard",{"transform":eulers[i]}))
			
#		proj.process_inplace("normalize.circlemean")
#		proj.mult(softmask)

		for eo in range(2):
			for j in range(nptcl[eo]):
				if classmx[eo][0,j]!=i : continue		# only proceed if the particle is in this class

				ptcl=EMData(args[0],j)
				
				# Find the transform for this particle (2d) and apply it to the unmasked/masked projections
				ptclxf=Transform({"type":"2d","alpha":cmxalpha[0,j],"mirror":int(cmxmirror[0,j]),"tx":cmxtx[0,j],"ty":cmxty[0,j]}).inverse()
				projc=[i.process("xform",{"transform":ptclxf}) for i in projs]		# we transform the projections, not the particle (as in the original classification)

				# make a filtered particle with 
				ctf=ptcl["ctf"]
				ds=1.0/(ctf.apix*ptcl["ny"])
				ssnr=ctf.compute_1d(ptcl["ny"],ds,Ctf.CtfType.CTF_SNR_SMOOTH,None)		# The smoothed curve
				plot(ssnr)

				
				
				projc.process_inplace("normalize.toimage",{"to":ptcl2})
				projc2.process_inplace("normalize.toimage",{"to":ptcl2})
				cmp1=ptcl2.cmp(simcmp[0],projc, simcmp[1])
				cmp2=ptcl2.cmp(simcmp[0],projc2,simcmp[1])
				result=cmp1-cmp2
				if options.debug: display((ptcl2,projc,projc2,proj))

				statr.append(result)
				statr2.append((cmp1+cmp2,cmp1-cmp2,0))

				if cmp1+cmp2<options.alistacks :
					ptcl2=ptcl.process("xform",{"transform":ptclxf.inverse()})
					ptcl2.mult(projmask)
					ptcl2.write_image("aligned_{}.hdf".format(i),nalis+2)

					if nalis==0 :
						projc=proj.process("normalize.toimage",{"to":ptcl2})
						projc2=proj2.process("normalize.toimage",{"to":ptcl2})
						projc.write_image("aligned_{}.hdf".format(i),0)
						projc2.write_image("aligned_{}.hdf".format(i),1)

					nalis+=1

				if options.tstcls==i :
					ptcl2=ptcl.process("xform",{"transform":ptclxf.inverse()})
					ptcl2.mult(projmask)
					if cmp1>cmp2:
						try: avgim1.add(ptcl2)
						except: avgim1=ptcl2
					else:
						try: avgim2.add(ptcl2)
						except: avgim2=ptcl2





if __name__ == "__main__":
	main()
