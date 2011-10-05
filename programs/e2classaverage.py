#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/27/2010 - rewritten almost from scratch
# Author: David Woolford, 9/7/2007 (woolford@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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

READ_HEADER_ONLY = True

from EMAN2db import EMTask, db_check_dict

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <output> [options]

	This program produces iterative class-averages, one of the secrets to EMAN's rapid convergence.
	Normal usage is to provide a stack of particle images and a classification matrix file defining
	class membership. Members of each class are then iteratively aligned to each other and averaged
	together with (optional) CTF correction.  It is also possible to use this program on all of the
	images in a single stack.

	"""
		
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input", type=str, help="The name of the input particle stack", default=None)
	parser.add_argument("--output", type=str, help="The name of the output class-average stack", default=None)
	parser.add_argument("--oneclass", type=int, help="Create only a single class-average. Specify the number.",default=None)
	parser.add_argument("--classmx", type=str, help="The name of the classification matrix specifying how particles in 'input' should be grouped. If omitted, all particles will be averaged.", default=None)
	parser.add_argument("--ref", type=str, help="Reference image(s). Used as an initial alignment reference and for final orientation adjustment if present. Also used to assign euler angles to the generated classes. This is typically the projections that were used for classification.", default=None)
	parser.add_argument("--storebad", action="store_true", help="Even if a class-average fails, write to the output. Forces 1->1 numbering in output",default=False)
	parser.add_argument("--resultmx",type=str,help="Specify an output image to store the result matrix. This contains 5 images where row is particle number. Rows in the first image contain the class numbers and in the second image consist of 1s or 0s indicating whether or not the particle was included in the class. The corresponding rows in the third, fourth and fifth images are the refined x, y and angle (respectively) used in the final alignment, these are updated and accurate, even if the particle was excluded from the class.", default=None)
	parser.add_argument("--iter", type=int, help="The number of iterations to perform. Default is 1.", default=1)
	parser.add_argument("--prefilt",action="store_true",help="Filter each reference (c) to match the power spectrum of each particle (r) before alignment and comparison",default=False)
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average. Default is None.", default=None)
	parser.add_argument("--aligncmp",type=str,help="The comparitor used for the --align aligner. Default is ccc.",default="ccc")
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	parser.add_argument("--raligncmp",type=str,help="The comparitor used by the second stage aligner.",default="ccc")
	parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average.",default="mean")
	parser.add_argument("--setsfref",action="store_true",help="This will impose the 1-D structure factor of the reference on the class-average (recommended when a reference is available)",default=False)
	parser.add_argument("--cmp",type=str,help="The comparitor used to generate quality scores for the purpose of particle exclusion in classes, strongly linked to the keep argument.", default="ccc")
	parser.add_argument("--keep",type=float,help="The fraction of particles to keep in each class.",default=1.0)
	parser.add_argument("--keepsig", action="store_true", help="Causes the keep argument to be interpreted in standard deviations.",default=False)
	parser.add_argument("--automask",action="store_true",help="Applies a 2-D automask before centering. Can help with negative stain data, and other cases where centering is poor.")
	parser.add_argument("--bootstrap",action="store_true",help="Ignored. Present for historical reasons only.")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is normalize.edgemean. If you want to turn this option off specify \'None\'", default="normalize.edgemean")
	parser.add_argument("--usefilt", dest="usefilt", default=None, help="Specify a particle data file that has been low pass or Wiener filtered. Has a one to one correspondence with your particle data. If specified will be used to align particles to the running class average, however the original particle will be used to generate the actual final class average")
	parser.add_argument("--idxcache", default=False, action="store_true", help="Ignored. Present for historical reasons.")
	parser.add_argument("--dbpath", help="Ignored. Present for historical reasons.", default=".")
	parser.add_argument("--resample",action="store_true",help="If set, will perform bootstrap resampling on the particle data for use in making variance maps.",default=False)
	parser.add_argument("--odd", default=False, help="Used by EMAN2 when running eotests. Includes only odd numbered particles in class averages.", action="store_true")
	parser.add_argument("--even", default=False, help="Used by EMAN2 when running eotests. Includes only even numbered particles in class averages.", action="store_true")
	parser.add_argument("--parallel", default=None, help="parallelism argument")
	parser.add_argument("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--debug","-d",action="store_true",help="Print debugging infromation while the program is running. Default is off.",default=False)
	parser.add_argument("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_argument("--check","-c",action="store_true",help="Performs a command line argument check only.",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if (options.check): options.verbose = 9 # turn verbose on if the user is only checking...

	error = check(options,True)

	if options.align : options.align=parsemodopt(options.align)
	if options.ralign : options.ralign=parsemodopt(options.ralign)
	if options.aligncmp : options.aligncmp=parsemodopt(options.aligncmp)
	if options.raligncmp : options.raligncmp=parsemodopt(options.raligncmp)
	if options.averager : options.averager=parsemodopt(options.averager)
	if options.cmp : options.cmp=parsemodopt(options.cmp)
	if options.normproc : options.normproc=parsemodopt(options.normproc)

	if options.resultmx!=None : options.storebad=True

	if (options.verbose>0):
		if (error):
			print "e2classaverage.py command line arguments test.... FAILED"
		else:
			print "e2classaverage.py command line arguments test.... PASSED"

	# returning a different error code is currently important to e2refine.py - returning 0 tells e2refine.py that it has enough
	# information to execute this script
	if error : exit(1)
	if options.check: exit(0)
	
	logger=E2init(sys.argv,options.ppid)
	
	try: 
		classmx=EMData.read_images(options.classmx)		# we keep the entire classification matrix in memory, since we need to update it in most cases
		ncls=int(classmx[0]["maximum"])
	except:
		ncls=1
		if options.resultmx!=None :
			print "resultmx can only be specified in conjunction with a valid classmx input."
			sys.exit(1)

	nptcl=EMUtil.get_image_count(options.input)

	# Initialize parallelism
	if options.parallel :
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
		pclist=[options.input]
		if options.ref: pclist.append(options.ref)
		if options.usefilt: pclist.append(options.usefilt)
		etc.precache(pclist)

	# prepare tasks
	tasks=[]
	if ncls>1:
		if options.oneclass==None : clslst=range(ncls)
		else : clslst=[options.oneclass]
		
		for cl in clslst:
			ptcls=classmx_ptcls(classmx[0],cl)			
			if options.resample : ptcls=[random.choice(ptcls) for i in ptcls]	# this implements bootstrap resampling of the class-average
			if options.odd : ptcls=[i for i in ptcls if i%2==1]
			if options.even: ptcls=[i for i in ptcls if i%2==0]
			tasks.append(ClassAvTask(options.input,ptcls,options.usefilt,options.ref,options.iter,options.normproc,options.prefilt,
			  options.align,options.aligncmp,options.ralign,options.raligncmp,options.averager,options.cmp,options.keep,options.keepsig,
			  options.automask,options.verbose,cl))
	
	else:
		ptcls=range(nptcl)
		if options.resample : ptcls=[random.choice(ptcls) for i in ptcls]
		if options.odd : ptcls=[i for i in ptcls if i%2==1]
		if options.even: ptcls=[i for i in ptcls if i%2==0]
		tasks.append(ClassAvTask(options.input,range(nptcl),options.usefilt,options.ref,options.iter,options.normproc,options.prefilt,
			  options.align,options.aligncmp,options.ralign,options.raligncmp,options.averager,options.cmp,options.keep,options.keepsig,
			  options.automask,options.verbose,0))

	# execute task list
	if options.parallel:				# run in parallel
		taskids=etc.send_tasks(tasks)
		alltaskids=taskids[:]

		while len(taskids)>0 :
			curstat=etc.check_task(taskids)
			for i,j in enumerate(curstat):
				if j==100 :
					rslt=etc.get_results(taskids[i])
					if rslt[1]["average"]!=None:
						rslt[1]["average"]["class_ptcl_src"]=options.input
						if options.ref!=None : rslt[1]["average"]["projection_image"]=options.ref
						if options.storebad : rslt[1]["average"].write_image(options.output,rslt[1]["n"])
						else: rslt[1]["average"].write_image(options.output,-1)
						
						
						# Update the resultsmx if requested
						if options.resultmx!=None:
							allinfo=rslt[1]["info"]				# the info result array list of (qual,xform,used) tuples
							pnums=rslt[0].data["images"][2]		# list of image numbers corresponding to information

							for n,info in enumerate(allinfo):
								y=pnums[n]		# actual particle number
								
								# find the matching class in the existing classification matrix
								for x in range(classmx[0]["nx"]):
									if classmx[0][x,y]==rslt[1]["n"] :		# if the class number in the classmx matches the current class-average number
										break
								else :
									print "Resultmx error: no match found ! (%d %d %d)"%(x,y,rslt[1]["n"])
									continue
								xform=info[1].get_params("2d")
								classmx[1][x,y]=info[2]					# used
								classmx[2][x,y]=xform["tx"]				# dx
								classmx[3][x,y]=xform["ty"]				# dy
								classmx[4][x,y]=xform["alpha"]			# da
								classmx[5][x,y]=xform["mirror"]			# flip
								try: classmx[6][x,y]=xform["scale"]
								except: pass
					# failed average
					elif options.storebad :
						blk=EMData(options.ref,0)
						blk=EMData(blk["nx"],blk["ny"],1)
						blk.to_zero()
						blk.set_attr("ptcl_repr", 0)
						blk.write_image(options.output,rslt[1]["n"])
					
			taskids=[j for i,j in enumerate(taskids) if curstat[i]!=100]
			
			if options.verbose and 100 in curstat : print "%d/%d tasks remain"%(len(taskids),len(alltaskids))
			time.sleep(3)
		
		
		if options.verbose : print "Completed all tasks"

	# single thread
	else:
		for t in tasks:
			rslt=t.execute()
			if rslt==None : sys.exit(1)
			
			if rslt["average"]!=None :
				rslt["average"]["class_ptcl_src"]=options.input
				if options.ref!=None : rslt["average"]["projection_image"]=options.ref
				if options.storebad : rslt["average"].write_image(options.output,t.options["n"])
				else: rslt["average"].write_image(options.output,-1)
				
				# Update the resultsmx if requested
				if options.resultmx!=None:
					allinfo=rslt["info"]				# the info result array list of (qual,xform,used) tuples
					pnums=t.data["images"][2]		# list of image numbers corresponding to information
					for n,info in enumerate(allinfo):
						y=pnums[n]		# actual particle number
						
						# find the matching class in the existing classification matrix
						for x in range(classmx[0]["nx"]):
							if classmx[0][x,y]==rslt["n"] :		# if the class number in the classmx matches the current class-average number
								break
						else :
							print "Resultmx error: no match found ! (%d %d %d)"%(x,y,rslt[1]["n"])
							continue
						xform=info[1].get_params("2d")
						classmx[1][x,y]=info[2]					# used
						classmx[2][x,y]=xform["tx"]				# dx
						classmx[3][x,y]=xform["ty"]				# dy
						classmx[4][x,y]=xform["alpha"]			# da
						classmx[5][x,y]=xform["mirror"]			# flip
						try: classmx[6][x,y]=xform["scale"]
						except: pass
						
			# Failed average
			elif options.storebad :
				blk=EMData(options.ref,0)
				blk=EMData(blk["nx"],blk["ny"],1)
				blk.to_zero()
				blk.set_attr("ptcl_repr", 0)
				blk.write_image(options.output,t.options["n"])
				
	if options.resultmx!=None:
		if options.verbose : print "Writing results matrix"
		for i,j in enumerate(classmx) : j.write_image(options.resultmx,i)

	E2end(logger)

class ClassAvTask(EMTask):
	"""This task will create a single task-average"""

	def __init__(self,imagefile,imagenums,usefilt=None,ref=None,niter=1,normproc=("normalize.edgemean",{}),prefilt=0,align=("rotate_translate_flip",{}),
		  aligncmp=("ccc",{}),ralign=None,raligncmp=None,averager=("mean",{}),scmp=("ccc",{}),keep=1.5,keepsig=1,automask=0,verbose=0,n=0):
		if usefilt==None : usefilt=imagefile
		data={"images":["cache",imagefile,imagenums],"usefilt":["cache",usefilt,imagenums]}
		if ref!=None : data["ref"]=["cache",ref,n]
		EMTask.__init__(self,"ClassAv",data,{},"")

		self.options={"niter":niter, "normproc":normproc, "prefilt":prefilt, "align":align, "aligncmp":aligncmp,
			"ralign":ralign,"raligncmp":raligncmp,"averager":averager,"scmp":scmp,"keep":keep,"keepsig":keepsig,
			"automask":automask,"verbose":verbose,"n":n}
	
	def execute(self,callback=None):
		"""This does the actual class-averaging, and returns the result"""
		options=self.options

		if options["verbose"]>0 : print "Start averaging class ",options["n"]

		try: ref=EMData(self.data["ref"][1],self.data["ref"][2])
		except: ref=None

#		print [self.data["images"][1]]+self.data["images"][2]

		# make the class-average
		try:
			avg,ptcl_info=class_average([self.data["usefilt"][1]]+self.data["usefilt"][2],ref,options["niter"],options["normproc"],options["prefilt"],options["align"],
				options["aligncmp"],options["ralign"],options["raligncmp"],options["averager"],options["scmp"],options["keep"],options["keepsig"],
				options["automask"],options["verbose"],callback)
		except KeyboardInterrupt: return None
		except SystemExit: return None
		except:
			return {"average":None,"info":None,"n":self.options["n"]}

		try: ref_orient=avg["xform.projection"]
		except: ref_orient=None

		try: ref_model=avg["model_id"]
		except: ref_model=0

		# Final alignment to the reference (if there is one)
		if ref!=None :
			#ref.process_inplace("normalize.edgemean")
			ali=align_one(avg,ref,True,self.options["align"],self.options["aligncmp"],self.options["ralign"],self.options["raligncmp"])
			fxf=ali["xform.align2d"]
			avg1=avg
			if options["verbose"]>0 : print "Final realign:",fxf
#			avg=class_average_withali([self.data["images"][1]]+self.data["images"][2],ptcl_info,Transform(),options["averager"],options["normproc"],options["verbose"])
#			avg.write_image("bdb:xf",-1)
			avg=class_average_withali([self.data["images"][1]]+self.data["images"][2],ptcl_info,fxf,options["averager"],options["normproc"],options["verbose"])
#			avg.write_image("bdb:xf",-1)

			#self.data["ref"].write_image("tst.hdf",-1)
			#avg1.write_image("tst.hdf",-1)
			#avg.write_image("tst.hdf",-1)
		else :
			# Nothing to align to, so we just try to center the final average
			ali=avg.process("xform.centerofmass",{"threshold":avg["mean"]+avg["sigma"]})
			fxf=ali["xform.align2d"]
			if options["verbose"]>0 : print "Final center:",fxf
			avg=class_average_withali([self.data["images"][1]]+self.data["images"][2],ptcl_info,Transform(),options["averager"],options["normproc"],options["verbose"])
			
		if ref_orient!=None: 
			avg["xform.projection"]=ref_orient
			avg["model_id"]=ref_model
			try: avg["projection_image_idx"]=self.data["ref"][2]
			except: pass
		
		return {"average":avg,"info":ptcl_info,"n":options["n"]}

def get_image(images,n,normproc=("normalize.edgemean",{})):
	"""used to get an image from a descriptor as provided to class_average function. Always a copy of the actual image."""
	if isinstance(images[0],EMData) : ret=images[n].copy()
	elif n>len(images)-2 : raise Exception, "get_image() outside range"
	else: ret=EMData(images[0],images[n+1])

	if normproc!=None : ret.process_inplace(normproc[0],normproc[1])
	
	return ret
	
def align_one(ptcl,ref,prefilt,align,aligncmp,ralign,raligncmp):
	"""Performs the multiple steps of a single particle-alignment"""

	if prefilt : ref.process_inplace("filter.matchto",{"to":ptcl})
	
	# initial alignment
	if align!=None :
		ali=ptcl.align(align[0],ref,align[1],aligncmp[0],aligncmp[1])

	# refine alignment if requested
	if ralign!=None:
		ralign[1]["xform.align2d"] = ali.get_attr("xform.align2d")
		ali=ptcl.align(ralign[0],ref,ralign[1],raligncmp[0],raligncmp[1])

	return ali

def class_average_withali(images,ptcl_info,xform,averager=("mean",{}),normproc=("normalize.edgemean",{}),verbose=0):
	"""This will generate a final class-average, given a ptcl_info list as returned by class_average,
	and a final transform to be applied to each of the relative transforms in ptcl_info. ptcl_info will
	be modified in-place to contain the aggregate transformations, and the final aligned average will be returned"""
	
	if isinstance(images[0],EMData) : nimg=len(images)
	elif isinstance(images[0],str) and isinstance(images[1],int) : nimg=len(images)-1
	else : raise Exception,"Bad images list"
	
	incl=[]
	excl=[]
	avgr=Averagers.get(averager[0], averager[1])
	for i in range(nimg):
		img=get_image(images,i,normproc)
		ptcl_info[i]=(ptcl_info[i][0],xform*ptcl_info[i][1],ptcl_info[i][2])		# apply the new Transform to the existing one
		img.process_inplace("xform",{"transform":ptcl_info[i][1]})
		try: use=ptcl_info[i][2]
		except: use=1
		if use : 
			avgr.add_image(img)				# only include the particle if we've tagged it as good
			if img.has_attr("source_n") : incl.append(img["source_n"])
		elif img.has_attr("source_n") : excl.append(img["source_n"])

	avg=avgr.finish()

	# set some useful attributes
	if len(incl)>0 or len(excl)>0 :
		if len(incl)>0 : avg["class_ptcl_idxs"]=incl
		if len(excl)>0 : avg["exc_class_ptcl_idxs"]=excl
		avg["class_ptcl_src"]=img["source_path"]
		
	return avg

def class_average(images,ref=None,niter=1,normproc=("normalize.edgemean",{}),prefilt=0,align=("rotate_translate_flip",{}),
		aligncmp=("ccc",{}),ralign=None,raligncmp=None,averager=("mean",{}),scmp=("ccc",{}),keep=1.5,keepsig=1,automask=0,verbose=0,callback=None):
	"""Create a single class-average by iterative alignment and averaging.
	images - may either be a list/tuple of images OR a tuple containing a filename followed by integer image numbers
	ref - optional reference image (EMData). 
	niter - Number of alignment/averaging iterations. If 0, will align to the reference with no further iterations.
	normproc - a processor tuple, normalization applied to particles before alignments
	prefilt - boolean. If set will 'match' reference to particle before alignment
	align - aligner tuple to align particle to averaged
	aligncmp - cmp for aligner
	ralign - aligner tuple for refining alignment
	raligncmp - cmp for ralign
	averager - averager tuple to generate class-average
	scmp - cmp tuple for comparing particle to reference for purposes of discarding bad particles
	keep - 'keep' value. Meaning depends on keepsig.
	keepsig - if set, keep is a 'sigma multiplier', otherwise keep is a fractional value (ie - 0.9 means discard the worst 10% of particles)
	
	returns (average,((cmp,xform,used),(cmp,xform,used),...))
	"""

	if verbose>2 : print "class_average(",images,ref,niter,normproc,prefilt,align,aligncmp,ralign,raligncmp,averager,scmp,keep,keepsig,automask,verbose,callback,")"

	# nimg is the number of particles we have to align/average
	if isinstance(images[0],EMData) : nimg=len(images)
	elif isinstance(images[0],str) and isinstance(images[1],int) : nimg=len(images)-1
	else : raise Exception,"Bad images list (%s)"%str(images)
	
	if verbose>2 : print "Average %d images"%nimg

	# If one image and no reference, just return it
	if nimg==1 and ref==None : return (get_image(images,0,normproc),[(0,Transform(),1)])
	
	# If one particle and reference, align and return
	if nimg==1:
		if averager[0]!="mean" : raise Exception,"Cannot perform correct average of single particle"
		ali=align_one(get_image(images,0,normproc),ref,prefilt,align,aligncmp,ralign,raligncmp)
		sim=ali.cmp(scmp[0],ref,scmp[1])			# compare similarity to reference (may use a different cmp() than the aligner)
		return (ali,[(sim,ali["xform.align2d"],1)])
	
	# If we don't have a reference image, we need to make one 
	if ref==None :
		if verbose : print "Generating reference"
#		sigs=[(get_image(i)["sigma"],i) for i in range(nimg)]		# sigma for each input image, inefficient
#		ref=get_image(images,max(sigs)[1])
		ref=get_image(images,0,normproc)										# just start with the first, as EMAN1

		# now align and average the set to the gradually improving average
		for i in range(1,nimg):
			if verbose>1 :
				print ".",
				sys.stdout.flush()
			ali=align_one(get_image(images,i,normproc),ref,prefilt,align,aligncmp,ralign,raligncmp)
			ref.add(ali)
		
		# A little masking and centering
		ref.process_inplace("normalize.circlemean")
		gmw=max(5,ref["nx"]/16)		# gaussian mask width
		ref.process_inplace("mask.gaussian",{"inner_radius":ref["nx"]/2-gmw,"outer_radius":gmw/1.3})
		ref.process_inplace("normalize.circlemean")
		ref.process_inplace("xform.centerofmass",{"threshold":ref["mean"]+ref["sigma"]})						# TODO: should probably check how well this works
		ref_orient=None
	else:
		try: ref_orient=ref["xform.projection"]
		except: ref_orient=None
		
		try: ref_model=ref["model_id"]
		except: ref_model=0
		
	if verbose>1 : print ""

	init_ref=ref.copy()

	# Iterative alignment
	ptcl_info=[None]*nimg		# empty list of particle info
	
	# This is really niter+1 1/2 iterations. It gets terminated 1/2 way through the final loop
	for it in range(niter+2):
		if verbose : print "Starting iteration %d"%it
		if callback!=None : callback(int(it*100/(niter+2)))
		
		# Evaluate quality from last iteration, and set a threshold for keeping particles
		if it>0:
			# measure statistics of quality values
			mean,sigma=0,0
			for sim,xf,use in ptcl_info:
				mean+=sim
				sigma+=sim**2
			mean/=len(ptcl_info)
			sigma=sqrt(sigma/len(ptcl_info)-mean**2)
			
			# set a threshold based on statistics and options
			if keepsig:					# keep a relative fraction based on the standard deviation of the similarity values
				thresh=mean+sigma*keep
				if verbose>1 : print "mean = %f\tsigma = %f\tthresh=%f"%(mean,sigma,thresh)
			else:						# keep an absolute fraction of the total
				l=[i[0] for i in ptcl_info]
				l.sort()
				try: thresh=l[int(len(l)*keep)]
				except:
					if verbose: print "Keeping all particles"
					thresh=l[-1]+1.0
			
			if verbose:
				print "Threshold = %1.4f   Quality: min=%f max=%f mean=%f sigma=%f"%(thresh,min(ptcl_info)[0],max(ptcl_info)[0],mean,sigma)
			
			# mark the particles to keep and exclude
			nex=0
			for i,pi in enumerate(ptcl_info):
				if pi[0]>thresh :
					nex+=1
					ptcl_info[i]=(pi[0],pi[1],0)
				elif pi[2]==0:
					ptcl_info[i]=(pi[0],pi[1],1)
			
			if verbose : print "%d/%d particles excluded"%(nex,len(ptcl_info))
			
			# if all of the particles were thrown out for some reason, we keep the best one
			if nex==len(ptcl_info) :
				best=ptcl_info.index(min(ptcl_info))
				ptcl_info[best]=(ptcl_info[best][0],ptcl_info[best][1],1)
				if verbose : print "Best particle reinstated"
				
		if it==niter+1 : break		# This is where the loop actually terminates. This makes sure that inclusion/exclusion is updated at the end
		
		# Now align and average
		avgr=Averagers.get(averager[0], averager[1])
		for i in range(nimg):
			if callback!=None and nimg%10==9 : callback(int((it+i/float(nimg))*100/(niter+2.0)))
			ptcl=get_image(images,i,normproc)					# get the particle to align 
			ali=align_one(ptcl,ref,prefilt,align,aligncmp,ralign,raligncmp)  # align to reference
			sim=ali.cmp(scmp[0],ref,scmp[1])			# compare similarity to reference (may use a different cmp() than the aligner)

			try: use=ptcl_info[i][2]
			except: use=1
			if use : 
				avgr.add_image(ali)				# only include the particle if we've tagged it as good
				if verbose>1 :
					sys.stdout.write(".")
					sys.stdout.flush()
			elif verbose>1:
				sys.stdout.write("X")
				sys.stdout.flush()
			ptcl_info[i]=(sim,ali["xform.align2d"],use)
			
		if verbose>1 : print ""

		ref=avgr.finish()
		
		# A little masking before the next iteration
		ref.process_inplace("normalize.circlemean")
		gmw=max(5,ref["nx"]/12)		# gaussian mask width
		if automask :
			ref.process_inplace("mask.auto2d",{"nmaxseed":10,"nshells":gmw-2,"nshellsgauss":gmw,"sigma":0.2})
		else :
			ref.process_inplace("mask.gaussian",{"inner_radius":ref["nx"]/2-gmw,"outer_radius":gmw/1.3})
		
	if ref_orient!=None : 
		ref["xform.projection"]=ref_orient
		ref["model_id"]=ref_model
	return [ref,ptcl_info]
	
def classmx_ptcls(classmx,n):
	"""Scans a classmx file to determine which images are in a specific class. classmx may be a filename or an EMData object.
	returns a list of integers"""
	
	if isinstance(classmx,str) : classmx=EMData(classmx,0)
	
	plist=[i.y for i in classmx.find_pixels_with_value(float(n))]
	
	return plist


def check(options,verbose=0):
	error = False

	if ( options.nofilecheck == False ):
		
		if options.output == None:
			print "Error: you must specify the output file"
			error = True
		elif file_exists(options.output):
			if not options.force:
				error = True
				if (verbose):
					print "Error: output file %s exists, force not specified, will not overwrite, exiting" %options.output
		
		if options.classmx == None:
			options.bootstrap = True # turn on boot strapping
		if options.classmx != None and not file_exists(options.classmx):
			error = True
			if (verbose):
				print "Error: the file expected to contain the classification matrix (%s) was not found, cannot run e2classaverage.py" %(options.classmx)
		
		if options.input == None:
			print "Error: you must specify the input file"
			error = True
		elif not file_exists(options.input):
			error = True
			if (verbose):
				print "Error:  failed to find the particle data (%s)" %options.input
		else:
			if (options.usefilt != None):
				if not file_exists(options.usefilt):
					error = True
					if verbose: print "Error: failed to find usefilt file %s" %options.usefilt
				
				n1 = EMUtil.get_image_count(options.usefilt)
				n2 = EMUtil.get_image_count(options.input)
				if n1 != n2:
					if verbose: print "Error, the number of images in the starting particle set:",n2,"does not match the number in the usefilt set:",n1
					error = True
					
				read_header_only=True
				img1 = EMData()
				img1.read_image(options.input,0,read_header_only)
				img2 = EMData()
				img2.read_image(options.usefilt,0,read_header_only)
				
				nx1 = img1.get_attr("nx") 
				nx2 = img2.get_attr("nx") 
				
				ny1 = img1.get_attr("ny") 
				ny2 = img2.get_attr("ny") 
				
				if nx1 != nx2 or ny1 != ny2:
					error = True
					if verbose: print "Error, the dimensions of particle data (%i x %i) and the usefilt data (%i x %i) do not match" %(nx1,ny1,nx2,ny2)
				
		
		if options.classmx != None and os.path.exists(options.classmx) and os.path.exists(options.input):
			(xsize, ysize ) = gimme_image_dimensions2D(options.classmx);
			numimg = EMUtil.get_image_count(options.input)
			if ( numimg != ysize ):
				error = True
				if (verbose):
					print "Error - the number of rows (%d) in the classification matrix image %s does not match the number of images (%d) in %s" %(ysize, options.classmx,numimg,options.input)
				
			
		if options.ref != None and not file_exists(options.ref):
			print "Error: the file expected to contain the reference images (%s) does not exist" %(options.ref)
			error = True
		elif options.ref and (os.path.exists(options.input) or db_check_dict(options.input)):
			(xsize, ysize ) = gimme_image_dimensions2D(options.input);
			(pxsize, pysize ) = gimme_image_dimensions2D(options.ref);
			if ( xsize != pxsize ):
				error = True
				if (verbose):
					print "Error - the dimensions of the reference and particle images do not match"
		elif options.classmx != None and options.ref:
				# classes contains the classifications - row is particle number, column data contains class numbers (could be greater than 1)
			classes = EMData()
			classes.read_image(options.classmx, 0,True)
			class_max = int(classes["maximum"])
			num_ref= EMUtil.get_image_count(options.ref)
			if ( class_max > num_ref ):
				print "Error, the classification matrix refers to a class number (%d) that is beyond the number of images (%d) in the reference image (%s)." %(class_max,num_ref,options.ref)
			
	if (options.iter > 1 or options.bootstrap) and options.align == None:
		print "Error: you must specify the align argument"
		error = True
		
	if ( options.even and options.odd ):
		print "Error, the even and odd arguments are mutually exclusive"
		error = True
	#if ( options.keep and options.keepsig ):
		#error = True
		#if ( verbose ):
			#print "Error: --keep and --keepsig are mutually exclusive"
	
	if ( options.keep > 1 or options.keep <= 0) and not options.keepsig :
		error = True
		if (verbose):
			print "The --keep option is a percentage expressed as a fraction - it must be between 0 and 1"
	
	if ( options.iter < 0 ):
		error = True
		if (verbose):
			print "Error, --iter must be greater than or equal to 0 - you specified %d" %(options.iter)
		
	if ( check_eman2_type(options.averager,Averagers,"Averager") == False ):
		if (verbose):
			print "Unknown averager",options.averager
		error = True
	
	if ( options.iter > 0 ):
		
		if ( check_eman2_type(options.cmp,Cmps,"Comparitor") == False ):
			error = True
			
		if (options.align == None):
			print "If --classiter is greater than zero, the -align argument must be specified"
			error = True
			
		if ( check_eman2_type(options.align,Aligners,"Aligner") == False ):
			error = True

		if ( check_eman2_type(options.aligncmp,Cmps,"Comparitor") == False ):
			error = True
		
		if ( options.ralign != None ):
			if ( check_eman2_type(options.ralign,Aligners,"Aligner") == False ):
				error = True
				
			if ( check_eman2_type(options.raligncmp,Cmps,"Comparitor") == False ):
				error = True
		if ( str(options.normproc) != "None" ):
			if ( check_eman2_type(options.normproc,Processors,"Processor") == False ):
				error = True


	if hasattr(options,"parallel") and options.parallel != None:
  		if len(options.parallel) < 2:
  			print "The parallel option %s does not make sense" %options.parallel
  			error = True

	return error
	
if __name__ == "__main__":
    main()
