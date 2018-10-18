#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

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

from past.utils import old_div
from builtins import range
import os
import sys
import random
import time
import string
import math
import queue
from os import system
from os import unlink
from sys import argv
import traceback
from EMAN2 import *
from numpy import arange

def monoextract(jsd,cls,threed,nptcl,mask,subunitmask,outermask,symxfs,eulers,cptcl,classmx,cmxtx,cmxty,cmxalpha,cmxmirror,domask,newbox,verbose):
	# for a single class-average, subtract away a projection of the reference with the exclusion mask in
	# each symmetry-related orientation, after careful scaling. Note that we don't have a list of which particle is in each class,
	# but rather a list of which class each particle is in, so we do this a bit inefficiently for now
	# subunitmask - when applied, includes only the density for the subunit to be extracted
	# mask - when applied, includes everything inside the normal 3-D whole model mask EXCEPT the subunit to be extracted
	# ptclmask - 
	if verbose : print("--- Class %d"%cls)
	ret=[]

	# We regenerate the masked volumes for each class to avoid using too much RAM

	# Projection of the full, unmasked 3-D reference
	# followed by projections excluding the subunit to be retained in each symmetry-based orientation
	projs=[threed.project("standard",{"transform":eulers[cls]})]
	
	# 2-D masks to be applied to the subparticles after subtraction
	projmask=[]
	for xf in symxfs:
		maskx=mask.process("xform",{"transform":xf})
		masked=threed.copy()
		masked.mult(maskx)
		projs.append(masked.project("standard",{"transform":eulers[cls]}))

		# these are the projection masks. There is one less of these since we don't have one for the unmasked ref volume
		# These include any region in 2-D which might conceivably contain subunit information
		# we need this for centering even if we don't apply the final mask
		if outermask!=None:
			omaskx=outermask.process("xform",{"transform":xf})
#				projmask.append(omaskx.project("standard",{"transform":eulers[cls]}))
			projmask.append(omaskx.project("standard",{"transform":eulers[cls]}).process("threshold.binary",{"value":0.01}))
		else:
			omaskx=subunitmask.process("xform",{"transform":xf})
			projmask.append(omaskx.project("standard",{"transform":eulers[cls]}).process("threshold.binary",{"value":0.01}))
		
	if verbose : print("--- Class %d masks prepared"%cls)
	
	for eo in range(2):
		for j in range(nptcl[eo]):
			if classmx[eo][0,j]!=cls : 
				#if options.debug: print("XXX {}\t{}\t{}\t{}".format(cls,("even","odd")[eo],j,classmx[eo][0,j]))
				continue		# only proceed if the particle is in this class
			if verbose>1: print("{}\t{}\t{}".format(cls,("even","odd")[eo],j))

			ptcl=EMData(cptcl[eo],j)
#				ptcl.write_image(options.output,-1)
			#lowth=ptcl["mean"]-ptcl["sigma"]*3.0
			#highth=ptcl["mean"]+ptcl["sigma"]*3.0

			# Find the transform for this particle (2d) and apply it to the unmasked/masked projections
			ptclxf=Transform({"type":"2d","alpha":cmxalpha[eo][0,j],"mirror":int(cmxmirror[eo][0,j]),"tx":cmxtx[eo][0,j],"ty":cmxty[eo][0,j]}).inverse()
			projc=[im.process("xform",{"transform":ptclxf}) for im in projs]		# we transform the projections, not the particle (as in the original classification)
			projm=[im.process("xform",{"transform":ptclxf}) for im in projmask]
			maskctr=[j.calc_center_of_mass(0.5) for j in projm]

			# now subtract the masked versions processed using the scaling/normalization
			# we got from the whole/unmasked particle
			for k,pr in enumerate(projc[1:]):
				ptcl3=ptcl.process("math.sub.optimal",{"ref":projc[0],"actual":pr})
				if domask :
					ptcl3.mult(projm[k])
#					projm[k].write_image(options.output,-1)
#					ptcl3.write_image(options.output,-1)
				ptcl3.translate(-maskctr[k][0]+ptcl3["nx"]//2,-maskctr[k][1]+ptcl3["ny"]//2,0)
				if newbox>3 : 
					ptcl4=ptcl3.get_clip(Region((ptcl3["nx"]-newbox)//2,(ptcl3["ny"]-newbox)//2,newbox,newbox))
					ret.append(ptcl4)
				else: ret.append(ptcl3)
#					print ptcl.cmp("optsub",projc[0])
	jsd.put((cls,ret))

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <refine_xx folder> <iteration number>

WARNING - particularly with a lot of threads and larger box-sizes, this program can eat a LOT of RAM. Try running on a small number of threads 
briefly, and monitor memory usage before ramping up to the full number of available threads.

This program modifies raw particle data to subtract-away undesired portions of the density on a per-particle basis.
For example, one could take GroEL particles with D7 symmetry and extract the top and bottom heptamers as separate
particles for reconstruction. In this example there would be 2x as many output particles as input particles. The
portion of the map identified by the excluded region mask will be subtracted from each particle as optimally as
possible, taking CTF, normalization and other factors into account.

To use this program you should first run a complete refinement to the best possible resolution. Refinement parameters
will be examined automatically to extract the corresponding particles and projections.

Since each input particle may possibly become multiple output particles (with --sym), no relationship to the original micrographs
is maintained, and a single output stack is generated by this program.
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--output",type=str,default=None,help="Specify output filename, which will contain nptcl*sym particles. Default=sets/subptcl_<refine>_<iter>.hdf")
	parser.add_argument("--sym", dest = "sym", default="c1",help = "Specify symmetry - choices are: c<n>, d<n>, tet, oct, icos. If specified each input 'particle' will become N extracted subparticles.")
	parser.add_argument("--subunitmask",type=str,default=None,help = "Required. When applied to the 3-D volume, this mask contains the subunit being extracted. 'soft' edges are permitted")
	parser.add_argument("--outermask",type=str,default=None,help = "Optional. If specified, this mask is projected into 2-D and used to mask out noise outside the subunit. If not specified a thresholded subunitmask is used. Only useful with --masked. ")
	parser.add_argument("--masked",action="store_true",default=False,help="If specified, each output subparticle will be masked based on the projection mask. Recommended.")
	parser.add_argument("--newbox", type=int, help="New box size for extracted regions",default=-1)
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--debug",action="store_true",default=False,help="Enable debugging mode with verbose output and image display. Not suitable for real runs.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	#parser.add_argument("--ncls","-N",type=int,help="Number of classes to generate",default=-1)
	#parser.add_argument("--average","-A",action="store_true",help="Average the particles within each class",default=False)

	(options, args) = parser.parse_args()
	args[1]=int(args[1])
	
	if len(args)!=2 : 
		print("The usage of this program has changed from previous versions. You now need to specify a 3-D mask indication the portion of the map you wish to KEEP (previously the mask idenfied a region to exclude). You may also optionally specify a 3-D mask indicating a slightly larger mask indicating the data to include in the final subtracted projection.")
		sys.exit(1)

	if options.subunitmask==None : 
		print("The --subunitmask option is required. Please see --help.")
		sys.exit(1)

	if options.output==None:
		options.output="sets/subptcl_{}_{}.hdf".format(args[0].split("_")[-1],args[1])
	outlst="sets/subptcl_{}_{}.lst".format(args[0].split("_")[-1],args[1])

	# each of these classmx variables becomes a 2 element list with even and odd particles respectively
	try:
		pathmx="{}/classmx_{:02d}_even.hdf".format(args[0],args[1])
		classmx=[EMData(pathmx,0)]
		nptcl=[classmx[0]["ny"]]
		cmxtx=[EMData(pathmx,2)]
		cmxty=[EMData(pathmx,3)]
		cmxalpha=[EMData(pathmx,4)]
		cmxmirror=[EMData(pathmx,5)]

		pathmx="{}/classmx_{:02d}_odd.hdf".format(args[0],args[1])
		classmx.append(EMData(pathmx,0))
		nptcl.append(classmx[1]["ny"])
		cmxtx.append(EMData(pathmx,2))
		cmxty.append(EMData(pathmx,3))
		cmxalpha.append(EMData(pathmx,4))
		cmxmirror.append(EMData(pathmx,5))
	except:
		traceback.print_exc()
		print("====\nError reading classification matrix. Must be full classification matrix with alignments")
		sys.exit(1)

	if options.verbose: print("{} even and {} odd particles in classmx".format(nptcl[0],nptcl[1]))

	# path to the even/odd particles used for the refinement
	cptcl=js_open_dict("{}/0_refine_parms.json".format(args[0]))["input"]
	cptcl=[str(i) for i in cptcl]

	# this reads all of the EMData headers from the projections, should be same for even and odd
	pathprj="{}/projections_{:02d}_even.hdf".format(args[0],args[1])
	nref=EMUtil.get_image_count(pathprj)
	eulers=[EMData(pathprj,i,1)["xform.projection"] for i in range(nref)]

	# The 3D reference volume we are using for subtraction
	threed=EMData("{}/threed_{:02d}.hdf".format(args[0],args[1]),0)

	# When applied this mask will include only the subunit to be extracted
	# note that it is not binarized, so "soft" edges are permitted
	subunitmask=EMData(options.subunitmask)
	
	# this is the "inverted" subunit mask, including everything except the subunit, again, "soft" edges are preserved
	mask=subunitmask.process("math.linear",{"scale":-1.0,"shift":1.0})
	ptclmask=EMData(args[0]+"/mask.hdf",0) # The full model mask from the reference volume
	mask.mult(ptclmask)
	# mask now excludes the subunit and volume outside the masked 3-D particle 

	# Projections of the 3-D outermask will be used in 2-D to eliminate noise
	# if not specified, the subunitmask will be projected (and thresholded) instead
	if options.outermask!=None:
		outermask=EMData(options.outermask)
	else: outermask=None

	logid=E2init(sys.argv, options.ppid)

	# We will need masked versions of our volume for all symmetric orientations
	# generate the symmetry transforms here for convenience
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(options.sym)
	symxfs=[xf.get_sym(options.sym,i) for i in range(nsym)]

	# clean up the output file
	try: os.unlink(options.output)
	except: pass

		
	NTHREADS=max(options.threads+1,2)		# we only run in threaded mode, and need one for collecting the results
	jsd=queue.Queue(0)

	# we start out with just a list of parameters to avoid initializing so many Thread objects at once
	thrds=[(jsd,cls,threed,nptcl,mask,subunitmask,outermask,symxfs,eulers,cptcl,classmx,cmxtx,cmxty,cmxalpha,cmxmirror,options.masked,options.newbox,options.verbose-1) for cls in range(nref)]

	# here we run the threads and save the results, no actual alignment done here
	if options.verbose: print(len(thrds)," threads")
	thrtolaunch=0
	ncmpl=0
	starttime=time.time()
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			if options.verbose>1 : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
			thrds[thrtolaunch]=threading.Thread(target=monoextract,args=thrds[thrtolaunch])
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)		# the final threads are running, we can wait a while
		
		# We do a running estimate in min/sec for the job to finish
		if options.verbose>1 and ncmpl>0: 
			frac=ncmpl/float(len(thrds))
			ETR=int(((tlf-starttime)/frac)-(time.time()-starttime))
			print("{:0.2f}% complete.  ETC- {:d}:{:02d}".format(100.0*frac,ETR//60,ETR%60))
	
		while not jsd.empty():
			ncmpl+=1			# number of completed tasks
			tlf=time.time()			# used for remaining time estimate
			cls,ptcls=jsd.get()
			for p in ptcls: 
				p.write_image(options.output,-1)
			thrds[cls].join()
			thrds[cls]=None

	os.system("e2proclst.py --create {} {}".format(outlst,options.output))
	print("Extraction complete")

	E2end(logid)






if __name__ == "__main__":
	main()
