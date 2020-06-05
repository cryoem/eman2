#!/usr/bin/env python
# Author: Steven Ludtke, 06/28/2018 (sludtke@bcm.edu)
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

def enrich(thr,jsd,imfile,lsx,proj,nenrich,redoinvar,i0,i1,verbose):
	
	ret=[]
	# we loop over the specified range of input images
	if verbose: print("Thread running on particles {} - {}".format(i0,i1))
	for i in range(i0,i1):
		an,af,ac=lsx.read(i)
		qary=np.matmul(proj,proj[i])		# matrix multiplication, highest values should be the most similar particles
		best=list(qary.argsort()[-nenrich-1:-1])
		best.reverse()			# first image should be best, working progressively may improve alignment?
		if verbose>1: print("{}: {}".format(i,best))
		avg=EMData(imfile,i)
		aliref=avg.copy()
		sim=[]
		for j,k in enumerate(best):
#			img=EMData(imfile,k).align("rotate_translate_tree",aliref,{"flip":1})
			img=EMData(imfile,k)
			ali=img.align("rotate_translate_tree",aliref,{"flip":1})
			ali=img.align("refine",aliref,{"xform.align2d":ali["xform.align2d"]},"frc",{"minres":80,"maxres":20})				
			sim.append(ali.cmp("frc",aliref,{"minres":80,"maxres":20}))
			avg.add(ali)
		avg.mult(old_div(1.0,(nenrich+1)))
		avg["class_ptcl_src"]=imfile
		avg["enrich_quals"]=sim
		avg["class_ptcl_idxs"]=[i]+best
		if avg["sigma"]==0 : print("Error: {} : {}".format(i,best))

		if redoinvar!=None:
			if redoinvar=="bispec":
				bspec=avg.process("math.bispectrum.slice",{"fp":bispec_invar_parm[1],"size":bispec_invar_parm[0]})
			elif redoinvar=="harmonic":
				bspec=avg.process("math.harmonic",{"fp":4})
			ret.append((avg,af.split("__")[0]+"__ctf_flip_enrich{}.hdf".format(nenrich),an,bspec,af.split("__")[0]+"__ctf_flip_invar.hdf"))
		else:
			ret.append((avg,af.split("__")[0]+"__ctf_flip_enrich{}.hdf".format(nenrich),an))


	if verbose: print("Thread complete for particles {} - {}".format(i0,i1))

	# signal that we're done
	jsd.put((thr,ret))

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2enrich.py <sets/particles.lst> [options]
This program will enrich each particle in a SPA project by averaging it with N similar particles. 
The additional particles are identified using an MSA subspace of the bispectra, which makes the process
fast enough to be practical. This program is best used with a set containing reasonably reliable 
particles, to avoid "bad particle" contamination. Particles must have run through standard ctf
autoprocessing prior to using this program, but no other processing is required. Suggest running this on
'fullres' particles.
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--nenrich", default=54,type=int,help="Number of additional particles to average with each particle. Default=5")
	parser.add_argument("--redoinvar",choices=["bispec","harmonic"],help="Recomputes invariants",default=None)
	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer.", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
#	parser.add_argument("--gui",action="store_true",help="Permits interactive adjustment of mask parameters",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode="tuning[True]")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	NTHREADS=max(options.threads+1,2)		# we have one controlling thread
	
	try: bispec=args[0].split("__ctf_flip")[0]+"__ctf_flip_invar.lst"
	except: bispec=args[0].rsplit(".",1)[0]+"_invar.hdf"
	n=EMUtil.get_image_count(args[0])
	if options.verbose: print(n," particles")
	step=max(1,n//10000)
	if step>1 : step="--step 0,{}".format(step)
	else: step=""
	
	logid=E2init(sys.argv, options.ppid)

	if options.verbose: print("Computing MSA of particle invariants ({})".format(bispec))
	# we start by running MSA on the full set of bispectra
	ret=launch_childprocess("e2msa.py {bispec} enrich_basis.hdf enrich_proj.hdf  --mode factan --nomean --nbasis=10 {step}".format(bispec=bispec,step=step))

	tmp=EMData("enrich_proj.hdf",0)
	ptcl_proj=to_numpy(tmp).copy()
	lsx=LSXFile(args[0],True)
	
	jsd=queue.Queue(0)
	nstep=n//50+1		# 50 particles at a time. Arbitrary
	thrds=[(i,jsd,args[0],lsx,ptcl_proj,options.nenrich,options.redoinvar,i*50,min(len(ptcl_proj),(i+1)*50),max(0,options.verbose-1)) for i in range(nstep)]
	
	if options.verbose: print("Beginning reprojections")
	# standard thread execution loop
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1 or not jsd.empty():
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
			
			# file writing only in the controlling thread
			if options.redoinvar:
				for im,fsp,i,bsp,bspfsp in rd[1]:
					im.write_image(fsp,i)
					bsp.write_image(bspfsp,i)
			else:
				for im,fsp,i in rd[1]:
					im.write_image(fsp,i)
				
			thrds[rd[0]].join()
			thrds[rd[0]]=None
			

	if options.verbose: print("All jobs complete")

	E2end(logid)


if __name__ == "__main__":
	main()

