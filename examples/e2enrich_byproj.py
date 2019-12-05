#!/usr/bin/env python
# Author: Steven Ludtke, 07/17/19 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
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

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import time
import os
import threading
import queue
from sys import argv,exit
import numpy as np

def enrich(thr,jsd,img,proj,mask,verbose):
	"cross correlates proj under all in-plane rotations with img"

	avgr=Averagers.get("minmax",{"max":1})
	for ang in range(0,359,5):
		# we set the mean to 0 and the standard deviation for unit vector length
		projr=proj.process("xform",{"alpha":ang})
		projr.add(-projr["mean"])
		projr.process_inplace("normalize.unitlen")
		# The mask, however, should be left with its normal 0-1 range
		maskr=mask.process("xform",{"alpha":ang}).process("normalize.unitlen")
		ccf=img.calc_ccf(projr)
		ccf.process_inplace("filter.highpass.gauss",{"cutoff_pixels":5})
		ccf.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/12.0})
		avgr.add_image(ccf.process("xform.phaseorigin.tocenter"))
#		ccf.process_inplace("normalize")	# now our values are more like Z-scores
		ccf.process_inplace("math.exp",{"low":ccf["sigma"]/2.0})
		ccf.process_inplace("mask.onlypeaks")
		wt=ccf.convolute(maskr)
		ccf=ccf.convolute(projr)
		try:
			ccfsum.add(ccf)
			wtsum.add(wt)
		except:
			ccfsum=ccf
			wtsum=wt
#		avgr.add_image(ccf)
#		print("{},{}\t{}\t{}\t{}".format(thr,ang,ccf["minimum"],ccf["maximum"],ccf["mean"]))

	if verbose: print("Thread complete for thread {}",thr)

	ret=[ccfsum,wtsum,avgr.finish()]
	# signal that we're done
	jsd.put((thr,ret))

def get_ptcl(arg):
	"Iterator generates EMData objects for input particle stack. Supports ',' specification of a single image"
	if "," in arg : yield EMData(arg.split(",")[0],int(arg.split(",")[1]))
	else:
		n=EMUtil.get_image_count(arg)
		for i in range(n): yield EMData(arg,i)

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2enrich_byproj.py <particles> <projections> <mask_projections> [options]
This program will use projections of a (usually single subunit) volume to 'reinterpret' each particle in the context of the subunit.
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#parser.add_argument("--nenrich", default=54,type=int,help="Number of additional particles to average with each particle. Default=5")
	#parser.add_argument("--redoinvar",choices=["bispec","harmonic"],help="Recomputes invariants",default=None)
	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer.", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
#	parser.add_argument("--gui",action="store_true",help="Permits interactive adjustment of mask parameters",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode="tuning[True]")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	NTHREADS=max(options.threads+1,2)		# we have one controlling thread
	
	proj=EMData.read_images(args[1])
	mask=EMData.read_images(args[2])
	if len(proj)!=len(mask) :
		print("Error: number of projections doesn't match number of mask projections")
		sys.exit(1)
	
	# centering
	mask=[m.process("xform.centerofmass") for m in mask]		# center of mass on the projection
	for i,p in enumerate(proj): p.transform(mask[i]["xform.align2d"])	# impose the same translation on the projections
	
	if options.verbose: print(len(proj)," projections")
	
	logid=E2init(sys.argv, options.ppid)

	for ptcl in get_ptcl(args[0]):
		if ptcl["nx"]!=proj[0]["nx"]:
			proj=[p.get_clip(Region((p["nx"]-ptcl["nx"])//2,(p["ny"]-ptcl["ny"])//2,ptcl["nx"],ptcl["ny"])) for p in proj]
			mask=[m.get_clip(Region((m["nx"]-ptcl["nx"])//2,(m["ny"]-ptcl["ny"])//2,ptcl["nx"],ptcl["ny"])) for m in mask]
		 
		jsd=queue.Queue(0)
		thrds=[(i,jsd,ptcl,proj[i],mask[i],max(0,options.verbose-1)) for i in range(len(proj))]
		
		
		if options.verbose: print("Beginning threads")
		# standard thread execution loop
		thrtolaunch=0
		while thrtolaunch<len(thrds) or threading.active_count()>1:
			if thrtolaunch<len(thrds):
				while (threading.active_count()>=options.threads) : time.sleep(0.1)
				if options.verbose>0 : 
					print("\r Starting thread {}/{}      ".format(thrtolaunch,len(thrds)), end=' ')
					sys.stdout.flush()
				thrds[thrtolaunch]=threading.Thread(target=enrich,args=thrds[thrtolaunch])		# replace args
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else: time.sleep(0.1)
			
			# return is [N,dict] a dict of image# keyed processed images
			while not jsd.empty():
				rd=jsd.get()
				rd[1][1].add(0.1)	# limits any problems with unweighted regions of the image
				
				try:
					sumim.add(rd[1][0])
					sumnorm.add(rd[1][1])
				except:
					sumim=rd[1][0]
					sumnorm=rd[1][1]
					
				rd[1][2].write_image("outz.hdf",rd[0])
#				rd[1][1].write_image("out2.hdf",rd[0])
				im=rd[1][0]/rd[1][1]
				im.write_image("out.hdf",rd[0])
									
				thrds[rd[0]].join()
				thrds[rd[0]]=None
				
		if options.verbose: print("All jobs complete")

		im=sumim/sumnorm
		im.write_image("final.hdf",-1)
		ptcl.write_image("input.hdf",-1)

	E2end(logid)


if __name__ == "__main__":
	main()

