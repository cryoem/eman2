#!/usr/bin/env python

#
# Author: Steven Ludtke, 4/2/2010
# Copyright (c) 2010 Baylor College of Medicine
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
from EMAN2db import db_open_dict
from math import *
import os
import sys
import datetime
from numpy import array
import traceback
import json

try:
	import numpy as np
	import matplotlib
	matplotlib.use("AGG")
#	matplotlib.use("PDF")
	import matplotlib.pyplot as plt
	pltcolors=["k","b","g","r","m","c","darkblue","darkgreen","darkred","darkmagenta","darkcyan","0.5"]
except:
	print "Matplotlib not available, plotting options will not be available"


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] [multi_xx]
	This program is still in its early stages. Eventually will provide a variety of tools for
	evaluating a single particle reconstruction refinement run."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--iter", type=int, default=None, help="If a refine_XX folder is being used, this selects a particular refinement iteration. Otherwise the last complete iteration is used.")
	parser.add_argument("--sym",type=str,help="Symmetry to be used in searching adjacent unit cells", default="c1")
	parser.add_argument("--maps",type=str,help="Comma separated list of maps against which particles will be compared.", default=None)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.maps:
		

	if args[0][:6]=="multi_":
		jsparm=js_open_dict(args[0]+"/0_refine_parms.json")
		
		if options.iter==None:
			try:
				options.iter=int(jsparm["last_map"][0].split("_")[-2][:2])
				options.sym=jsparm["sym"]
			except:
				print "Could not find a completed iteration in ",args[0]
				sys.exit(1)
	
	try:
		pathmx="{}/classmx_{:02d}.hdf".format(args[0],options.iter)
		classmx=EMData(pathmx,0)
		nptcl=classmx["ny"]
		cmxtx=EMData(pathmx,2)
		cmxty=EMData(pathmx,3)
		cmxalpha=EMData(pathmx,4)
		cmxmirror=EMData(pathmx,5)
	except:
		traceback.print_exc()
		print "====\nError reading classification matrix. Must be full classification matrix with alignments"
		sys.exit(1)
	
	if options.verbose: print "{} particles in classmx".format(nptcl)

	# path to the even/odd particles used for the refinement
	cptcl=str(jsparm["input"])
	
	# this reads all of the EMData headers from the projections, should be same for even and odd
	pathprj="{}/projections_{:02d}.hdf".format(args[0],options.iter)
	nref=EMUtil.get_image_count(pathprj)
	eulers=[EMData(pathprj,i,1)["xform.projection"] for i in range(nref)]

	# The 3D reference volumes we are using for subtraction
	maps=[]
	for lm in jsparm["last_map"]:
		d = EMData("{}/{}".format(args[0],lm)) #options.iter,len(jsparm["last_map"]),i)
		maps.append(d)
	
	nx = maps[0]["nx"]
	apix = maps[0]["apix_x"]
	restarget = 2*apix+0.25 # maybe incorrect? could get approx target from fsc_mutual_avg_02.txt
	automaskexpand = int(nx/20)
		
	# The masks applied to the reference volumes, used for 2-D masking of particles for better power spectrum matching
	masks=[]
	for threed in maps:
		# New version of automasking based on a more intelligent interrogation of the volume
		vol=threed.copy()
		
		vol.process_inplace("filter.lowpass.gauss",{"cutoff_freq":min(0.1,1.0/restarget)})		# Mask at no higher than 10 A resolution
		md=vol.calc_radial_dist(nx/2,0,1,3)	# radial max value per shell in real space
		rmax=int(nx/2.2)		# we demand at least 10% padding
		vmax=max(md[:rmax])			# max value within permitted radius
		# this finds the first radius where the max value @ r falls below overall max/4
		# this becomes the new maximum mask radius
		act=0
		mv=0,0
		for i in xrange(rmax):
			if md[i]>mv[0] : mv=md[i],i		# find the radius of the  max val in range
			if not act and md[i]<0.9*vmax : continue
			act=True
			if md[i]<0.2*vmax :
				rmax=i
				break
		rmaxval=mv[1]
		vmax=mv[0]
		
		# excludes any spurious high values at large radius
		vol.process_inplace("mask.sharp",{"outer_radius":rmax})
		
		# automask
		mask=vol.process("mask.auto3d",{"threshold":vmax*.2,"radius":0,"nshells":int(nx*0.05+.5+automaskexpand),"nshellsgauss":int(restarget*1.5/apix),"nmaxseed":24,"return_mask":1})
		
		# We expand the mask a bit, since we want to consider problems with "touching" particles
		mask.process_inplace("threshold.binary",{"value":0.2})
		mask.process_inplace("mask.addshells",{"nshells":nx//15})
		mask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})
		
		masks.append(mask)
	
	rings=[int(2*nx*apix/res) for res in (100,30,15,8,4)]
	print rings
	nbands = len(rings)-1

	fout=open("ptclfsc_multi_{}.txt".format(args[0][-2:]),"w")
	# generate a projection for each particle so we can compare

	pj = 0
	for i in xrange(nref):
		if options.verbose>1 : print "--- Class %d/%d"%(i,nref-1)
		
		phi=eulers[i].get_rotation("eman")["phi"]
		alt=eulers[i].get_rotation("eman")["alt"]
		az=eulers[i].get_rotation("eman")["az"]

		for j in xrange(nptcl):
			if classmx[0,j]!=i : continue	# only proceed if the particle is in this class
			if options.verbose > 6: print "{}\t{}".format(i,j)

			#from IPython import embed
			#embed()

			# the particle itself
			try: ptcl=EMData(cptcl,j)
			except:
				print "Unable to read particle: {} ({})".format(cptcl,j)
				sys.exit(1)
			try: defocus=ptcl["ctf"].defocus
			except: defocus=-1.0
			
			data = []
			for threed,ptclmask in zip(maps,masks):
				
				# The first projection is unmasked, used for scaling
				proj=threed.project("standard",{"transform":eulers[i]})
				projmask=ptclmask.project("standard",eulers[i])	# projection of the 3-D mask for the reference volume to apply to particles
				
				# Find the transform for this particle (2d) and apply it to the unmasked/masked projections
				ptclxf=Transform({"type":"2d","alpha":cmxalpha[0,j],"mirror":int(cmxmirror[0,j]),"tx":cmxtx[0,j],"ty":cmxty[0,j]}).inverse()
				projc=proj.process("xform",{"transform":ptclxf})	# we transform the projection, not the particle (as in the original classification)
			
				projmaskc=projmask.process("xform",{"transform":ptclxf})
				ptcl.mult(projmaskc)
	
				# Particle vs projection FSC
				fsc = ptcl.calc_fourier_shell_correlation(projc)
	
				third = len(fsc)/3
				fsc=array(fsc[third:third*2])
				for k in xrange(nbands): # sum the fsc into 5 range values
					s = sum(fsc[rings[k]:rings[k+1]])/(rings[k+1]-rings[k])
					data.append(str(s))
			
			# to which model does this particle belong according to similarity at each resolution range
			best = [np.argmax(data[k::nbands]) for k in xrange(nbands)]
			
			# which model "wins" majority of times?
			counts = [0,0,0,0]
			for b in best:
				if b == 0: counts[0]+=1
				elif b == 1: counts[1]+=1
				elif b == 2: counts[2]+=1
				elif b == 3: counts[3]+=1
			winner = np.argmax(counts)
			
			for d in [best[0],best[1],best[2],best[3],winner,phi,alt,az,i,defocus]:
				data.append(str(d))
			
			dat = "\t".join(data)
			cmt = "\t# {};{}\n".format(j,cptcl)
			fout.write(dat+cmt)
			
			pj+=1

	print "Results in ptclfsc_multi_{}.txt".format(args[0][-2:])
	sys.exit(0)

if __name__ == "__main__":
    main()
