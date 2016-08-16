#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/11/15 
# Copyright (c) 2015- Baylor College of Medicine
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
import math
from copy import deepcopy
import os
import sys
import random
from random import choice
import traceback

from EMAN2jsondb import JSTask,jsonclasses

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <output> [options]

	The goal of this program is to reduce the heterogeneity of a reconstruction by splitting a single map
	into two maps, each more homogeneous. You must run e2refine_easy to completion before using this program.
	It will take the class-averaging results from the final iteration, and split the particles from each 
	class-average into 2 groups, producing 2 class-averages for each. The program then attempts to construct
	a maximally self-consistent grouping of these pairs of class averages into 2 3-D maps. 
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--path", default=None, type=str,help="The name of an existing refine_xx folder, where e2refine_easy ran to completion",guitype='filebox', filecheck=False,browser="EMBrowserWidget(withmodal=True,multiselect=False)", row=3, col=0, rowspan=1, colspan=3)
	parser.add_argument("--usebasis", default=0,type=int,help="Select which Eigenimage to use for separation. With novarimax, n=0 is highest energy.", guitype='intbox', row=5, col=0, rowspan=1, colspan=1)
	parser.add_argument("--nbasis", default=-1,type=int,help="Number of basis vectors to compute. Must be at least usebasis+1. Default 6 or usebasis+1.", guitype='intbox', row=4, col=0, rowspan=1, colspan=1)
	parser.add_argument("--novarimax", action="store_true",default=False, help="Disable varimax rotation among computed basis vectors.",guitype='boolbox', row=7, col=0, rowspan=1, colspan=1)
	parser.add_argument("--mask", default=None, help="Optional 3D mask to focus the classification", guitype='filebox', browser='EMSetsTable(withmodal=True,multiselect=False)', filecheck=False, row=6, col=0, rowspan=1, colspan=3, mode="refinement")
	parser.add_argument("--parallel", default="thread:2", help="Standard parallelism option. Default=thread:2", guitype='strbox', row=8, col=0, rowspan=1, colspan=2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.nbasis<=1 :
		options.nbasis=6
		if options.nbasis<=options.usebasis+1 :
			options.nbasis=options.usebasis+1
			print "--nbasis adjusted to ",options.nbasis

	if options.path==None:
		paths=[i for i in os.listdir(".") if "refine_" in i and len(i)==9]
		paths.sort()
		options.path=paths[-1]

	pathnum=options.path[-2:]

	# check the specified path for the files we need
	try:
		olddb = js_open_dict(options.path+"/0_refine_parms.json")
		last_map=olddb["last_map"]
		targetres=olddb["targetres"]
		last_iter=int(last_map.split("_")[-1][:2])
		try: 
			ptcls=olddb["inputavg"]
			if ptcls==None : raise Exception
		except: ptcls=olddb["input"]
		
		sym=olddb["sym"]
		if olddb["breaksym"]:
			sym="c1"
		if options.verbose : print "Found iteration {} in {}, using {}".format(last_iter,options.path," & ".join(ptcls))
	except:
		traceback.print_exc()
		print "Error: Cannot find necessary files in ",options.path
		sys.exit(1)
		
	logger=E2init(sys.argv,options.ppid)

	# classmx is a list with 2 elements. Each element is a list of EMData from the corresponding cls_result file
	classmx=[]
	classmx.append(EMData.read_images("{}/cls_result_{:02d}_even.hdf".format(options.path,last_iter)))
	classmx.append(EMData.read_images("{}/cls_result_{:02d}_odd.hdf".format(options.path,last_iter)))
	ncls=max(int(classmx[0][0]["maximum"])+1,int(classmx[1][0]["maximum"])+1)

	# Rearrange the info in classmx
	classlists=[[] for i in xrange(ncls)]	# empty list for each class
	
	# This will produce a list of particles with Transforms for each class
	for eo in (0,1):
		for y in xrange(classmx[eo][0]["ny"]):
			ptcl=[eo,y,Transform({"type":"2d","tx":classmx[eo][2][0,y],"ty":classmx[eo][3][0,y],"alpha":classmx[eo][4][0,y],"mirror":int(classmx[eo][5][0,y])})]
			#print ptcl, 
			#print int(classmx[eo][0][0,y])
			classlists[int(classmx[eo][0][0,y])].append(ptcl)
	
	#if len(classlists[0])>100 :
		#print "Warning: this program is normally intended for use with downsampled data and fairly coarse angular sampling. If you try to use it with a large number of class-averages you may have a variety of problems, and should insure that your machine has sufficient RAM."
		

	# Initialize parallelism
	from EMAN2PAR import EMTaskCustomer
	etc=EMTaskCustomer(options.parallel)

	# Empty image to pad classes file
	zero=EMData(str(ptcls[0]),0)
	zero.to_zero()
	zero["ptcl_repr"]=0
	
	# Euler angles for averages
	projin="{}/projections_{:02d}_even.hdf".format(options.path,last_iter)
	eulers=[EMData(projin,i,True)["xform.projection"] for i in xrange(ncls)]
	
	# Prepare mask if specified
	if options.mask!=None:
		mask=EMData(options.mask)
		
	else : mask=None
	
	# prepare tasks
	tasks=[]
	gc=0
	ns=[classmx[eo][0]["ny"] for eo in (0,1)]
	for c,cl in enumerate(classlists):
		if len(cl)<20 : 							# we require at least 20 particles in a class to make the attempt
#			zero.write_image(classout[0],c)
#			zero.write_image(classout[1],c)
			continue
		if mask!=None :
			maskp=mask.project("standard",eulers[c])
		else: maskp=None
		tasks.append(ClassSplitTask(ptcls,ns,cl,c,eulers[c],maskp,options.usebasis,options.nbasis,options.novarimax,options.verbose-1))
		gc+=1
	
#	for t in tasks: t.execute()

	# execute task list
	taskids=etc.send_tasks(tasks)
	alltaskids=taskids[:]

	classes=[]
	while len(taskids)>0 :
		curstat=etc.check_task(taskids)
		for i,j in enumerate(curstat):
			if j==100 :
				rslt=etc.get_results(taskids[i])
				rsltd=rslt[1]
				cls=rslt[0].options["classnum"]
				if rsltd.has_key("failed") :
					print "Bad average in ",cls
				else:
					#rsltd["avg1"].write_image(classout[0],cls)
					#rsltd["avg2"].write_image(classout[1],cls)
					ncls=rsltd["avg1"]["ptcl_repr"]+rsltd["avg2"]["ptcl_repr"]
					# note that the 2 results we get back are in arbitrary order!
					# the next section of code with 3D reconstruction is designed to sort out
					# which average should be paired with which
					classes.append([ncls,rsltd["avg1"]["xform.projection"],rsltd["avg1"],rsltd["avg2"],rsltd["basis"],cls])	# list of (ptcl_repr,xform,avg1,avg2)
				
		taskids=[j for i,j in enumerate(taskids) if curstat[i]!=100]

		if options.verbose and 100 in curstat :
			print "%d/%d tasks remain"%(len(taskids),len(alltaskids))
		if 100 in curstat :
			E2progress(logger,1.0-(float(len(taskids))/len(alltaskids)))

	if options.verbose : print "Completed all tasks\nGrouping consistent averages"

	classes.sort(reverse=True)		# we want to start with the largest number of particles
	apix=classes[0][2]["apix_x"]

	boxsize=classes[0][2]["ny"]
	pad=good_size(boxsize*1.5)
	if options.verbose: print "Boxsize -> {}, padding to {}".format(boxsize,pad)
		
	# a pair of reconstructors. we will then simultaneously reconstruct in the pair, and use each to decide on the best target for each particle
	recon=[Reconstructors.get("fourier",{"size":[pad,pad,pad],"sym":sym,"mode":"gauss_5"}) for i in (0,1)]
	for r in recon: r.setup()
	
	# We insert the first class-average (with the most particles) randomly into reconstructor 1 or 2
	p2=classes[0][2].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
	p3=recon[0].preprocess_slice(p2,classes[0][1])
	recon[0].insert_slice(p3,classes[0][1],classes[0][2].get_attr_default("ptcl_repr",1.0))

	p2=classes[0][3].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
	p3=recon[1].preprocess_slice(p2,classes[0][1])
	recon[1].insert_slice(p3,classes[0][1],classes[0][3].get_attr_default("ptcl_repr",1.0))
	
	classes[0].append(0)

	if options.verbose : print "Reconstruction: pass 1"
	for i,c in enumerate(classes[1:]):
		proj=EMData(projin,c[5])		# the projection corresponding to this average
		# while this does cost us a final interpolation, high resolution isn't the primary aim anyway, and getting the alignment consistent is important
		# also gives us a chance to normalize
		c[2]["xform.align2d"]=Transform()
		ali2=c[2].align("refine",proj)
		ali2.process_inplace("normalize.toimage",{"to":proj,"ignore_zero":1})
		c[3]["xform.align2d"]=Transform()
		ali3=c[3].align("refine",proj)
		ali3.process_inplace("normalize.toimage",{"to":proj,"ignore_zero":1})
		
#		print "ROT:\t",ali2["xform.align2d"].get_params("2d"),"\t",ali3["xform.align2d"].get_params("2d")
		
		# note that ali2 and c[2] are the same except for a final alignment
		a2=ali2.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))		# first class-average
		a3=recon[0].preprocess_slice(a2,classes[0][1])
		a3n=c[2].get_attr_default("ptcl_repr",1.0)
		
		# similarly ali3 and c[3] are the same
		b2=ali3.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
		b3=recon[1].preprocess_slice(b2,classes[0][1])						# I don't believe it matters if we use recon[0] or 1 here, but haven't checked
		b3n=c[3].get_attr_default("ptcl_repr",1.0)
		
		recon[0].determine_slice_agreement(a3,c[1],a3n,False)
#		print a3.get_attr_dict()
		q0a=a3["reconstruct_absqual_lowres"]		# quality for average a in reconstruction0
#		n0a=a3["reconstruct_norm"]			# normalization for same
		
		recon[1].determine_slice_agreement(a3,c[1],a3n,False)
		q1a=a3["reconstruct_absqual_lowres"]		# quality for average a in reconstruction0
#		n1a=a3["reconstruct_norm"]			# normalization for same
		
		recon[0].determine_slice_agreement(b3,c[1],b3n,False)
		q0b=b3["reconstruct_absqual_lowres"]		# quality for average a in reconstruction0
#		n0b=b3["reconstruct_norm"]			# normalization for same
		
		recon[1].determine_slice_agreement(b3,c[1],b3n,False)
		q1b=b3["reconstruct_absqual_lowres"]		# quality for average a in reconstruction0
#		n1b=b3["reconstruct_norm"]			# normalization for same
		
		if options.verbose>1 : print i,q0a,q1a,q0b,q1b,q0a+q1b,q1a+q0b
		if options.verbose>2 : print "\t\t",n0a,n1a,n0b,n1b
			
		if q0a+q1b>q1a+q0b :		# if true, a -> recon0 and b -> recon1 
			c.append(0)				# we put a 0 at the end of the classes element if we use a->0,b->1 ordering, 1 if swapped
#			a3.mult(n0a)
			recon[0].insert_slice(a3,c[1],a3n)
#			b3.mult(n1b)
			recon[1].insert_slice(b3,c[1],b3n)
		else:
			c.append(1)
#			a3.mult(n1a)
			recon[1].insert_slice(a3,c[1],a3n)
#			b3.mult(n0b)
			recon[0].insert_slice(b3,c[1],b3n)

	if options.verbose : print "Reconstruction: pass 2"
	
	# another pass with the filled reconstruction to make sure our initial assignments were ok
#	for i,c in enumerate(classes[1:]):
#		a2=c[2].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))		# first class-average
#		a3=recon[0].preprocess_slice(a2,classes[0][1])
#		a3n=c[2].get_attr_default("ptcl_repr",1.0)
#		
#		b2=c[3].get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
#		b3=recon[1].preprocess_slice(b2,classes[0][1])						# I don't believe it matters if we use recon[0] or 1 here, but haven't checked
#		b3n=c[3].get_attr_default("ptcl_repr",1.0)
#		
#		recon[0].determine_slice_agreement(a3,c[1],a3n,0) # c[-1]==0
#		q0a=a3["reconstruct_absqual"]			# quality for average a in reconstruction0
#		n0a=a3["reconstruct_norm"]			# normalization for same
#		
#		recon[1].determine_slice_agreement(a3,c[1],a3n,0) # c[-1]==1
#		q1a=a3["reconstruct_absqual"]			# quality for average a in reconstruction0
#		n1a=a3["reconstruct_norm"]			# normalization for same
#		
#		recon[0].determine_slice_agreement(b3,c[1],b3n,0) # c[-1]==1
#		q0b=b3["reconstruct_absqual"]			# quality for average a in reconstruction0
#		n0b=b3["reconstruct_norm"]			# normalization for same
#		
#		recon[1].determine_slice_agreement(b3,c[1],b3n,0) # c[-1]==0
#		q1b=b3["reconstruct_absqual"]			# quality for average a in reconstruction0
#		n1b=b3["reconstruct_norm"]			# normalization for same
#		
#		if options.verbose>1 : print i,q0a,q1a,q0b,q1b,q0a+q1b,q1a+q0b
#			
#		if q0a+q1b>q1a+q0b :		# if true, a -> recon0 and b -> recon1 
#			if c[-1]==1 :
#				c[-1]=0
#				print i," 1->0"
#			
#			c.append(0)				# we put a 0 at the end of the classes element if we use a->0,b->1 ordering, 1 if swapped
#			a3.mult(n0a)
#			recon[0].insert_slice(a3,c[1],a3n)
#			b3.mult(n1b)
#			recon[1].insert_slice(b3,c[1],b3n)
#		else:
#			if c[-1]==0 :
#				c[-1]=1
#				print i," 0->1"
#
#			c.append(1)
#			a3.mult(n1a)
#			recon[1].insert_slice(a3,c[1],a3n)
#			b3.mult(n0b)
#	
#		
	if options.verbose : print "All done, writing output"

	if mask!=None: msk="_msk"
	else: msk=""
	classout=["{}/classes_{:02d}_bas{}{}_split0.hdf".format(options.path,last_iter,options.usebasis,msk),"{}/classes_{:02d}_bas{}{}_split1.hdf".format(options.path,last_iter,options.usebasis,msk)]
	basisout="{}/classes_{:02d}{}_basis".format(options.path,last_iter,msk)
	threedout="{}/threed_{:02d}{}_split.hdf".format(options.path,last_iter,msk)
	threedout2="{}/threed_{:02d}{}_split_filt_bas{}.hdf".format(options.path,last_iter,msk,options.usebasis)
	setout=["sets/split_{}{}_bas{}_0.lst".format(pathnum,msk,options.usebasis),"sets/split_{}{}_bas{}_1.lst".format(pathnum,msk,options.usebasis)]
	split=[r.finish(True).get_clip(Region((pad-boxsize)/2,(pad-boxsize)/2,(pad-boxsize)/2,boxsize,boxsize,boxsize)) for r in recon]
	split[0]["apix_x"]=apix
	split[0]["apix_y"]=apix
	split[0]["apix_z"]=apix
	split[1]["apix_x"]=apix
	split[1]["apix_y"]=apix
	split[1]["apix_z"]=apix
	split[0].process_inplace("mask.soft",{"outer_radius":-8,"width":4})
	split[1].process_inplace("mask.soft",{"outer_radius":-8,"width":4})
	split[0].write_image(threedout,0)
	split[1].write_image(threedout,1)

	# now we write the class-averages and the new (split) particle files
	lstin =[LSXFile(ptcls[0],True),LSXFile(ptcls[1],True)]
	try:
		os.unlink("sets/split0.lst")
		os.unlink("sets/split1.lst")
	except: pass
	lstout=[LSXFile("sets/split0.lst"),LSXFile("sets/split1.lst")]
	for i,c in enumerate(classes):
		c[2].write_image(classout[c[-1]],i)	# class-average
		ptcln=c[2]["class_eoidxs"]		# eofile/ptcl# pairs
		for p in xrange(0,len(ptcln),2):
			lstout[0][-1]=lstin[ptcln[p]][ptcln[p+1]]		# wierd syntax, but the -1 here appends
			
		c[3].write_image(classout[c[-1]^1],i)
		ptcln=c[3]["class_eoidxs"]		# eofile/ptcl# pairs
		for p in xrange(0,len(ptcln),2):
			lstout[1][-1]=lstin[ptcln[p]][ptcln[p+1]]		# wierd syntax, but the -1 here appends

		if options.verbose>2:
			c[4][0].write_image(basisout+"1.hdf",i)
			c[4][1].write_image(basisout+"2.hdf",i)
			c[4][2].write_image(basisout+"3.hdf",i)

	launch_childprocess("e2proclst.py sets/split0.lst --mergesort {}".format(setout[0]))
	launch_childprocess("e2proclst.py sets/split1.lst --mergesort {}".format(setout[1]))

	try:
		os.unlink("sets/split0.lst")
		os.unlink("sets/split1.lst")
	except:
		pass

	if os.path.exists("strucfac.txt"):
		launch_childprocess("e2proc3d.py {} {} --setsf strucfac.txt --process filter.wiener.byfsc:fscfile={}/fsc_masked_{:02d}.txt:snrmult=2:sscale=1.1:maxfreq={} --process mask.soft:outer_radius=-9:width=4".format(threedout,threedout2,options.path,last_iter,1.0/targetres))
	else:
		print "Missing structure factor, cannot filter properly"
		launch_childprocess("e2proc3d.py {} {} --process filter.wiener.byfsc:fscfile={}/fsc_masked_{:02d}.txt:snrmult=2:sscale=1.1:maxfreq={} --process mask.soft:outer_radius=-9:width=4".format(threedout,threedout2,options.path,last_iter,1.0/targetres))

	E2end(logger)

class ClassSplitTask(JSTask):
	"""This task will create a single class-average"""

	def __init__(self,ptclfiles,ns,ptcls,nc,euler,mask,usebasis,nbasis,novarimax,verbose):
		"""ptclfiles is a list of 2 (even/odd) particle stacks. ns is the number of particles in each of ptcfiles. ptcls is a list of lists containing [eo,ptcl#,Transform]"""
#		sys.stderr=file("task.err","a")
		data={"particles1":["cache",ptclfiles[0],(0,ns[0])],"particles2":["cache",ptclfiles[1],(0,ns[1])]}
		JSTask.__init__(self,"ClassSplit",data,{},"")

		self.options={"particles":ptcls,"classnum":nc,"euler":euler,"usebasis":usebasis,"novarimax":novarimax,"nbasis":nbasis,"mask":mask,"verbose":verbose}

	def execute(self,callback=None):
		"""This does the actual class-averaging, and returns the result"""
		options=self.options

		if options["verbose"]>0 : print "Start averaging class {} with {} particles ".format(options["classnum"],len(options["particles"]))
		
		files=self.data["particles1"][1],self.data["particles2"][1]
		
		ptcls=options["particles"]			# just a shortcut

		#if len(options["particles"])<5 :
			#z=EMData(str(files[ptcls[0][0]]),ptcls[0][1]).to_zero()
			#return {"avg":z,"basis":[z,z,z,z,z]}
		 		
#		print files,ptcls[0]
		# read in all particles and append each to element to ptcls
		avgr=Averagers.get("mean")
		for p in ptcls: 
			p.append(EMData(str(files[p[0]]),p[1]).process("xform",{"transform":p[2]}))
			p.append(p[-1].process("filter.highpass.gauss",{"cutoff_freq":0.01})
				.process("filter.lowpass.gauss",{"cutoff_freq":0.05})
				.process("normalize.circlemean",{"radius":-6})
				.process("mask.soft",{"outer_radius":-8,"width":4}))
			avgr.add_image(p[4])
		
		# Copy each particle minus it's mean value
		avg=avgr.finish()
		
		# PCA on the mean-subtracted particles
		# At this point p[3] will be the particle in the correct orientation
		# p[4] will be the filtered/masked particle
		# p[5] will be the filtered/masked/bg subtr particle
		for p in ptcls: 
			p.append(p[4].copy())
			p[5].sub(avg)

		if options["mask"]==None:
			mask=ptcls[0][-1].copy()
			mask.to_one()
#			mask.process("mask.soft",{"outer_radius":-8,"width":4})
			mask.process("mask.sharp",{"outer_radius":-10})
		else:
			mask=options["mask"]

#		print "basis start"
		pca=Analyzers.get("pca_large",{"nvec":options["nbasis"],"mask":mask,"tmpfile":"tmp{}".format(options["classnum"])})
		for p in ptcls: 
			pca.insert_image(p[5])		# filter to focus on lower resolution differences
		basis=pca.analyze()

		# Varimax rotation... good idea?
		if not options["novarimax"]:
			pca2=Analyzers.get("varimax",{"mask":mask})
			for im in basis:
				pca2.insert_image(im)
		
			basis=pca2.analyze()

		
		# if you turn this on multithreaded it will crash sometimes
		#avg.mult(0.05)	#arbitrary for debugging
		#avg.write_image("pca.hdf",-1)
		#basis[0].write_image("pca.hdf",-1)
		#basis[1].write_image("pca.hdf",-1)
		#basis[2].write_image("pca.hdf",-1)
		#basis[3].write_image("pca.hdf",-1)
#		print "basis"
		
		# at the moment we are just splitting into 2 classes, so we'll use the first eigenvector. A bit worried about defocus coming through, but hopefully ok...
		dots=[p[5].cmp("ccc",basis[self.options["usebasis"]]) for p in ptcls]	# NOTE: basis number is passed in as an option, may not be #1 or #3 (default)
		if len(dots)==0:
			return {"failed":True}
		dota=sum(dots)/len(dots)
		
#		print "average"
		# we will just use the sign of the dot product to split
		avgr=[Averagers.get("mean"),Averagers.get("mean")]
		incl=[[],[]]
		for i,d in enumerate(dots):
			if d<dota: 
				avgr[0].add_image(ptcls[i][3])
				incl[0].append(ptcls[i][0])
				incl[0].append(ptcls[i][1])
			else: 
				avgr[1].add_image(ptcls[i][3])
				incl[1].append(ptcls[i][0])
				incl[1].append(ptcls[i][1])
		
		#for p in ptcls: 
			#avgr.add_image(p[3].process("xform",{"transform":p[2]}))
		
		if options["verbose"]>0: print "Finish averaging class {}".format(options["classnum"])
#		if callback!=None : callback(100)
#		return {"avg":avg,"basis":basis}
		try:
			avg1=avgr[0].finish()
			avg1["xform.projection"]=options["euler"]
			avg1["class_eoidxs"]=incl[0]			# contains alternating file # and particle #
#			print avg1
			avg2=avgr[1].finish()
			avg2["xform.projection"]=options["euler"]
			avg2["class_eoidxs"]=incl[1]
#			print avg2
		except:
			return {"failed":True}

#		print basis
		return {"avg1":avg1,"avg2":avg2,"basis":basis[:3]}

jsonclasses["ClassSplitTask"]=ClassSplitTask.from_jsondict


if __name__ == "__main__":
    main()
