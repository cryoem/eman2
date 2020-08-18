#!/usr/bin/env python
# align all particles to reference and store alignment results
# Author: Steven Ludtke (sludtke@bcm.edu)
# Copyright (c) 2019- Baylor College of Medicine
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
# 09/29/2019  Steve Ludtke


from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import time
import os
import threading
import queue
from sys import argv,exit
from EMAN2jsondb import JSTask
import numpy as np
import sklearn.decomposition as skdc

def ptclextract(jsd,db,ks,shrink,layers,sym,verbose):
	#  we have to get the 3d particles in the right orientation
	lasttime=time.time()
	for i,k in ks:
		# this only happens in the first thread
		if verbose and time.time()-lasttime>3:
			print("\r  {}/{}       ".format(i,ks[-1][0]),end="")
			sys.stdout.flush()
			lasttime=time.time()
			
		parm=db[k[2]]		# alignment for one particle
		ptcl=EMData(k[0],k[1])
		xf=parm["xform.align3d"]
		ptcl.process_inplace("xform",{"transform":xf})
		if shrink>1 : ptcl.process_inplace("math.meanshrink",{"n":shrink})
		if sym!="" and sym!="c1" : ptcl.process_inplace("xform.applysym",{"sym":sym})
		ptcl["score"]=parm["score"]
		
		# these are the range limited orthogonal projections
		x=ptcl.process("misc.directional_sum",{"axis":"x","first":ptcl["nx"]/2-layers,"last":ptcl["nx"]/2+layers+1})
		y=ptcl.process("misc.directional_sum",{"axis":"y","first":ptcl["nx"]/2-layers,"last":ptcl["nx"]/2+layers+1})
		z=ptcl.process("misc.directional_sum",{"axis":"z","first":ptcl["nx"]/2-layers,"last":ptcl["nx"]/2+layers+1})

		# different directions sometimes have vastly different standard deviations, independent normalization may help balance
		x.process_inplace("normalize")
		y.process_inplace("normalize")
		z.process_inplace("normalize")
		
		# we pack the 3 projections into a single 2D image
		all=EMData(x["nx"]*3,x["ny"],1)
		all.insert_clip(x,(0,0))
		all.insert_clip(y,(x["nx"],0))
		all.insert_clip(z,(x["nx"]*2,0))
		all["score"]=parm["score"]
		
		all["orig_file"]=ptcl["data_source"]
		all["orig_n"]=ptcl["data_n"]
		all["orig_key"]=k[2]
		jsd.put((i,k,all))


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_classify_byproj.py --ncls <n> --path <spt_XX> --iter <N> 
This program is part of the 'new' hierarchy of e2spt_ programs.

This program will generate a set of 3 orthogonal projections of the central slices of each aligned particle from an existing spt_XX result, 
then use this triplet of projections to perform a k-means classification on the particles. In addition to the class-averages, will also 
produce new sets/ for each class, which could be further-refined.
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--path",type=str,default=None,help="Path to an existing spt_XX folder with the alignment results to use, defualt = highest spt_XX",guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=0, col=0, rowspan=1, colspan=2,mode="gui")
	parser.add_argument("--iter",type=int,help="Iteration number within path, default = last iteration",default=-1,guitype="intbox", row=1, col=0, rowspan=1, colspan=1,mode="gui")
	parser.add_argument("--ncls",type=int,help="Number of classes to generate",default=3,guitype="intbox", row=1, col=1, rowspan=1, colspan=1,mode="gui")
	parser.add_argument("--nbasis",type=int,help="Number of basis vectors for the MSA phase, default=4",default=4)
	parser.add_argument("--layers",type=int,help="number of slices about the center to use for the projection in each direction, ie 0->1, 1->3, 2->5. Default=2",default=2,guitype="intbox", row=2, col=1, rowspan=1, colspan=1,mode="gui")	
	parser.add_argument("--sym",type=str,default="c1",help="Symmetry of the input. Must be aligned in standard orientation to work properly.")
	parser.add_argument("--shrink", default=1,type=int,help="shrink the particles before processing",guitype="intbox", row=2, col=0, rowspan=1, colspan=1,mode="gui")
	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.", guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--saveali",action="store_true",help="In addition to the unaligned sets/ for each class, generate aligned particle stacks per class",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
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
		if options.verbose : print("Using iteration ",options.iter)

	db=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))

	logid=E2init(sys.argv, options.ppid)

	# the keys are name,number pairs for the individual particles
	ks=[eval(k)+(k,) for k in db.keys()]
	ks.sort()
	ks=list(enumerate(ks))
	
	prjs=[0]*len(ks)
	if options.verbose: print("Generating sections")
	lasttime=time.time()
	jsd=queue.Queue(0)

	NTHREADS=max(options.threads,2)		# we have one thread just writing results
	thrds=[threading.Thread(target=ptclextract,args=(jsd,db,ks[i::NTHREADS-1],options.shrink,options.layers,options.sym,options.verbose>1 and i==0)) for i in range(NTHREADS-1)]

	try: os.unlink("{}/alisecs_{:02d}.hdf".format(options.path,options.iter))
	except: pass
	# here we run the threads and save the results, no actual alignment done here
	if options.verbose: print(len(thrds)," threads")
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1 or not jsd.empty():
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			if options.verbose>1 : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
		#if options.verbose>1 and thrtolaunch>0:
			#frac=thrtolaunch/float(len(thrds))
			#print("{}% complete".format(100.0*frac))
	
		while not jsd.empty():
			i,k,all=jsd.get()
			prjs[i]=all
			all.write_image("{}/alisecs_{:02d}.hdf".format(options.path,options.iter),i)

	for t in thrds:
		t.join()

	if options.verbose: print("Performing PCA")

	mask=sum(prjs)
	mask.process_inplace("threshold.notzero")
	prjsary=np.zeros((len(prjs),int(mask["square_sum"]+0.5)))		# input to PCA
	for i,p in enumerate(prjs):
		pp=p.process("misc.mask.pack",{"mask":mask})	# pack 3D unmasked values into 1D
		prjsary[i]=to_numpy(pp)

	# run PCA and reproject in one step, then we will classidy the projections
	P=skdc.PCA(options.nbasis)
	prjsprjs=[from_numpy(i) for i in P.fit_transform(prjsary)]
		
	if options.verbose: print("Classifying")
	# classification
	an=Analyzers.get("kmeans")
	an.set_params({"ncls":options.ncls,"maxiter":100,"minchange":len(ks)//(options.ncls*25),"verbose":options.verbose-1,"slowseed":1,"outlierclass":0,"mininclass":2})
	
	an.insert_images_list(prjs)
	centers=an.analyze()

	if options.verbose: print("\nGenerating 3-D class-averages")
	# Generate new averages from the original particles
	# Write new sets/ for each class
	for i in range(options.ncls):
		try: os.unlink("sets/{}_{:02d}_{:02d}.lst".format(options.path,options.iter,i))
		except: pass

	hdr=EMData(ks[0][1][0],ks[0][1][1],True)
	nx=hdr["nx"]
	ny=hdr["ny"]
	nz=hdr["nz"]
	sets=[LSXFile("sets/{}_{:02d}_{:02d}.lst".format(options.path,options.iter,i)) for i in range(options.ncls)]
	
	# Initialize averages
	avgs=[EMData(nx,ny,nz) for i in range(options.ncls)]
	for a in avgs:
		a["apix_x"]=hdr["apix_x"]
		a["apix_y"]=hdr["apix_y"]
		a["apix_z"]=hdr["apix_z"]
	
	# do the actual averaging
	for n,im in enumerate(prjs):
		if options.verbose and time.time()-lasttime>3:
			print("\r  {}/{}       ".format(n+1,len(prjs)),end="")
			sys.stdout.flush()
			lasttime=time.time()
			
		cls=im["class_id"]
		# read the original volume again, and fix its orientation (again)
		ptcl=EMData(im["orig_file"],im["orig_n"])
		xf=db[im["orig_key"]]["xform.align3d"]
		ptcl.transform(xf)
		avgs[cls].add(ptcl)		# add to the correct class average
		
		# If requested, we save the aligned (shrunken) volume to the stack for this class
		if options.saveali: 
			if options.shrink>1 : ptcl.process_inplace("math.meanshrink",{"n":options.shrink})
			ptcl.write_image("{}/aliptcl_{:02d}_{:02d}.hdf".format(options.path,options.iter,cls),-1)
		
		sets[cls].write(-1,im["orig_n"],im["orig_file"])
	
	if options.verbose: print("\nSaving classes")
	try: os.unlink("{}/classes_sec_{:02d}.hdf".format(options.path,options.iter))
	except: pass
	try: os.unlink("{}/classes_{:02d}.hdf".format(options.path,options.iter))
	except: pass
	# write class averages
	for i in range(options.ncls):
		n=centers[i]["ptcl_repr"]
		print("Class {}: {}".format(i,n))
		centers[i].write_image("{}/classes_sec_{:02d}.hdf".format(options.path,options.iter),i)
		if n>0 : avgs[i].mult(1.0/n)
		avgs[i].write_image("{}/classes_{:02d}.hdf".format(options.path,options.iter),i)
		
	if options.verbose: print("Done")

	E2end(logid)

if __name__ == "__main__":
	main()









