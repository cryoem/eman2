#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
# align all particles to reference and store alignment results
# Author: Steven Ludtke (sludtke@bcm.edu)
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

def alifn(jsd,fsp,i,a,options):
	t=time.time()
	b=EMData(fsp,i).do_fft()
	b.process_inplace("xform.phaseorigin.tocorner")

	# we align backwards due to symmetry
	if options.verbose>2 : print("Aligning: ",fsp,i)
	c=a.xform_align_nbest("rotate_translate_3d_tree",b,{"verbose":0,"sym":options.sym,"sigmathis":0.1,"sigmato":1.0, "maxres":options.maxres,"wt_ori":options.wtori},options.nsoln)
	for cc in c : cc["xform.align3d"]=cc["xform.align3d"].inverse()

	jsd.put((fsp,i,c[0]))
	if options.verbose>1 : print("{}\t{}\t{}\t{}".format(fsp,i,time.time()-t,c[0]["score"]))

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_align.py [options] <subvolume_stack> <reference>
Note that this program is not part of the original e2spt hierarchy, but is part of an experimental refactoring.

This program will take an input stack of subtomograms and a reference volume, and perform a missing-wedge aware alignment of each particle to the reference. If --goldstandard is specified, then even and odd particles will be aligned to different perturbed versions of the reference volume, phase-randomized past the specified resolution."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--goldstandard",type=float,help="If specified, will phase randomize the even and odd references past the specified resolution (in A, not 1/A)",default=0)
	parser.add_argument("--goldcontinue",action="store_true",help="Will use even/odd refs corresponding to specified reference to continue refining without phase randomizing again",default=False)
	parser.add_argument("--saveali",action="store_true",help="Save a stack file (aliptcls.hdf) containing the aligned subtomograms.",default=False)
	parser.add_argument("--savealibin",type=int,help="shrink aligned particles before saving",default=1)
	parser.add_argument("--path",type=str,default=None,help="Path to a folder where results should be stored, following standard naming conventions (default = spt_XX)")
	parser.add_argument("--sym",type=str,default="c1",help="Symmetry of the input. Must be aligned in standard orientation to work properly.")
	parser.add_argument("--maxres",type=float,help="Maximum resolution to consider in alignment (in A, not 1/A)",default=0)
	parser.add_argument("--wtori",type=float,help="Weight for using the prior orientation in the particle header. default is -1, i.e. not used.",default=-1)
	parser.add_argument("--nsoln",type=int,help="number of solutions to keep at low resolution for the aligner",default=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="spt_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : fls=[0]
		options.path = "spt_{:02d}".format(max(fls)+1)
		try: os.mkdir(options.path)
		except: pass

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : options.iter=1
		else: options.iter=max(fls)+1

	reffile=args[1]
	NTHREADS=max(options.threads+1,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)

	if options.goldcontinue:
		ref=[]
		try:
			ref.append(EMData(reffile[:-4]+"_even.hdf",0))
			ref.append(EMData(reffile[:-4]+"_odd.hdf",0))
		except:
			print("Error: cannot find one of reference files, eg: ",EMData(reffile[:-4]+"_even.hdf",0))
	else:
		ref=[]
		ref.append(EMData(reffile,0))
		ref.append(EMData(reffile,0))

		if options.goldstandard>0 :
			ref[0].process_inplace("filter.lowpass.randomphase",{"cutoff_freq":old_div(1.0,options.goldstandard)})
			ref[1].process_inplace("filter.lowpass.randomphase",{"cutoff_freq":old_div(1.0,options.goldstandard)})
			ref[0].write_image("{}/align_ref.hdf".format(options.path),0)
			ref[1].write_image("{}/align_ref.hdf".format(options.path),1)

	ref[0]=ref[0].do_fft()
	ref[0].process_inplace("xform.phaseorigin.tocorner")
	ref[1]=ref[1].do_fft()
	ref[1].process_inplace("xform.phaseorigin.tocorner")

	angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))
	jsd=queue.Queue(0)

	n=-1
	N=EMUtil.get_image_count(args[0])
	thrds=[threading.Thread(target=alifn,args=(jsd,args[0],i,ref[i%2],options)) for i in range(N)]

	# here we run the threads and save the results, no actual alignment done here
	print(len(thrds)," threads")
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			if options.verbose : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)

		while not jsd.empty():
			fsp,n,d=jsd.get()
			angs[(fsp,n)]=d
			if options.saveali:
				v=EMData(fsp,n)
				v.transform(d["xform.align3d"])
				if options.savealibin>1:
					v.process_inplace("math.meanshrink",{"n":options.savealibin})
				v.write_image("{}/aliptcls_{:02d}.hdf".format(options.path, options.iter),n)


	for t in thrds:
		t.join()

	E2end(logid)


if __name__ == "__main__":
	main()

