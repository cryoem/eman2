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

	parser.add_argument("--path", default=None, type=str,help="The name of an existing refine_xx folder, where e2refine_easy ran to completion")
	parser.add_argument("--parallel", default="thread:2", help="Standard parallelism option. Default=thread:2")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	logger=E2init(sys.argv,options.ppid)

	# Initialize parallelism
	if options.parallel :
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
		pclist=[options.input]
		etc.precache(pclist)

	# prepare tasks
	tasks=[]
		if options.oneclass==None : clslst=range(ncls)
		else : clslst=[options.oneclass]

		for cl in clslst:
			ptcls=classmx_ptcls(classmx[0],cl)
			tasks.append(ClassAvTask(options.input,ptcls,options.usefilt,options.ref,options.iter,options.normproc,options.prefilt,
			  options.align,options.aligncmp,options.ralign,options.raligncmp,options.averager,options.cmp,options.keep,options.keepsig,
			  options.automask,options.saveali,options.setsfref,options.verbose,cl,options.center))

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
						if options.decayedge:
							nx=rslt[1]["average"]["nx"]
							rslt[1]["average"].process_inplace("normalize.circlemean",{"radius":nx/2-nx/15})
							rslt[1]["average"].process_inplace("mask.gaussian",{"inner_radius":nx/2-nx/15,"outer_radius":nx/20})
							#rslt[1]["average"].process_inplace("mask.decayedge2d",{"width":nx/15})

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

			taskids=[j for i,j in enumerate(taskids) if curstat[i]!=100]

			if options.verbose and 100 in curstat :
				print "%d/%d tasks remain"%(len(taskids),len(alltaskids))
			if 100 in curstat :
				E2progress(logger,1.0-(float(len(taskids))/len(alltaskids)))

			time.sleep(3)


		if options.verbose : print "Completed all tasks"


	print "Class averaging complete"
	E2end(logger)

class ClassAvTask(JSTask):
	"""This task will create a single task-average"""

	def __init__(self,imagefile,imagenums,usefilt=None,ref=None,niter=1,normproc=("normalize.edgemean",{}),prefilt=0,align=("rotate_translate_flip",{}),
		  aligncmp=("ccc",{}),ralign=None,raligncmp=None,averager=("mean",{}),scmp=("ccc",{}),keep=1.5,keepsig=1,automask=0,saveali=0,setsfref=0,verbose=0,n=0,center="xform.center"):
		if usefilt==None : usefilt=imagefile
		self.center=center
		data={"images":["cache",imagefile,imagenums],"usefilt":["cache",usefilt,imagenums]}
		if ref!=None : data["ref"]=["cache",ref,n]
		JSTask.__init__(self,"ClassAv",data,{},"")

		self.options={"niter":niter, "normproc":normproc, "prefilt":prefilt, "align":align, "aligncmp":aligncmp,
			"ralign":ralign,"raligncmp":raligncmp,"averager":averager,"scmp":scmp,"keep":keep,"keepsig":keepsig,
			"automask":automask,"saveali":saveali,"setsfref":setsfref,"verbose":verbose,"n":n}

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
				options["automask"],options["saveali"],options["verbose"],callback)
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
			# This was commented out because the class-average doesn't have CTF parameters (or shouldn't) and
			# often will be using a comparator which makes use of CTF. Hard-coding the aligner for now
			#ali=align_one(avg,ref,True,self.options["align"],self.options["aligncmp"],self.options["ralign"],self.options["raligncmp"])
#			ali=align_one(avg,ref,True,("rotate_translate_flip_iterative",{}),("ccc",{}),("refine",{}),("ccc",{}))
			# changed to this in 3/6/14 because it was causing class-averages done without flipping to sometimes become flipped. Also not sure if I trust the _iterative aligner
			ali=align_one(avg,ref,True,("rotate_translate",{}),("ccc",{}),("refine",{}),("ccc",{}))
			fxf=ali["xform.align2d"]
			avg1=avg
			if options["verbose"]>0 : print "Final realign:",fxf
#			avg=class_average_withali([self.data["images"][1]]+self.data["images"][2],ptcl_info,Transform(),options["averager"],options["normproc"],options["verbose"])
#			avg.write_image("bdb:xf",-1)
			avg=class_average_withali([self.data["images"][1]]+self.data["images"][2],ptcl_info,fxf,ref,options["averager"],options["normproc"],options["setsfref"],options["verbose"])
#			avg.write_image("bdb:xf",-1)

			#self.data["ref"].write_image("tst.hdf",-1)
			#avg1.write_image("tst.hdf",-1)
			#avg.write_image("tst.hdf",-1)
		else :
			# Nothing to align to, so we just try to center the final average
			#gmw=max(5,avg["nx"]/16)		# gaussian mask width
			#avg.process_inplace("filter.highpass.gauss",{"cutoff_pixels":min(avg["nx"]/10,5)})	# highpass to reduce gradient issues
			#avg.process_inplace("normalize.circlemean")
			#avg.process_inplace("mask.gaussian",{"inner_radius":avg["nx"]/2-gmw,"outer_radius":gmw/1.3})
			#avg.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.07})
			#avg.process_inplace("normalize.circlemean")
			#ali=avg.process("threshold.binary",{"value":avg["mean"]+avg["sigma"]*1.5})
			#ali.process_inplace("xform.centerofmass",{"threshold":0.5})
			ali=avg.process(self.center)
			fxf=ali["xform.align2d"]
			if options["verbose"]>0 : print "Final center:",fxf.get_trans_2d()
			avg1=avg
			avg=class_average_withali([self.data["images"][1]]+self.data["images"][2],ptcl_info,fxf,None,options["averager"],options["normproc"],options["setsfref"],options["verbose"])
		try:
			avg["class_ptcl_qual"]=avg1["class_ptcl_qual"]
			avg["class_ptcl_qual_sigma"]=avg1["class_ptcl_qual_sigma"]
		except: pass

		if ref_orient!=None:
			avg["xform.projection"]=ref_orient
			avg["model_id"]=ref_model
			try: avg["projection_image_idx"]=self.data["ref"][2]
			except: pass

		return {"average":avg,"info":ptcl_info,"n":options["n"]}

jsonclasses["ClassAvTask"]=ClassAvTask.from_jsondict

def get_image(images,n,normproc=("normalize.edgemean",{})):
	"""used to get an image from a descriptor as provided to class_average function. Always a copy of the actual image."""
	if isinstance(images[0],EMData) : ret=images[n].copy()
	elif n>len(images)-2 : raise Exception, "get_image() outside range"
	else: ret=EMData(images[0],images[n+1])

	if normproc!=None : ret.process_inplace(normproc[0],normproc[1])

	return ret

def align_one(ptcl,ref,prefilt,align,aligncmp,ralign,raligncmp):
	"""Performs the multiple steps of a single particle-alignment"""

	if prefilt : ref=ref.process("filter.matchto",{"to":ptcl})

	# initial alignment
	if align!=None :
		ali=ptcl.align(align[0],ref,align[1],aligncmp[0],aligncmp[1])

	# refine alignment if requested
	if ralign!=None:
		ralign[1]["xform.align2d"] = ali.get_attr("xform.align2d")
		ali=ptcl.align(ralign[0],ref,ralign[1],raligncmp[0],raligncmp[1])

	return ali

def class_average_withali(images,ptcl_info,xform,ref,averager=("mean",{}),normproc=("normalize.edgemean",{}),setsfref=0,verbose=0):
	"""This will generate a final class-average, given a ptcl_info list as returned by class_average,
	and a final transform to be applied to each of the relative transforms in ptcl_info. ptcl_info will
	be modified in-place to contain the aggregate transformations, and the final aligned average will be returned"""

	if isinstance(images[0],EMData) : nimg=len(images)
	elif isinstance(images[0],str) and isinstance(images[1],int) : nimg=len(images)-1
	else : raise Exception,"Bad images list"

	incl=[]
	excl=[]
#	xforms=[]
	avgr=Averagers.get(averager[0], averager[1])
	for i in range(nimg):
		img=get_image(images,i,normproc)
		ptcl_info[i]=(ptcl_info[i][0],xform*ptcl_info[i][1],ptcl_info[i][2])		# apply the new Transform to the existing one
#		ptcl_info[i]=(ptcl_info[i][0],ptcl_info[i][1]*xform,ptcl_info[i][2])		# apply the new Transform to the existing one
		img.process_inplace("xform",{"transform":ptcl_info[i][1]})
		try: use=ptcl_info[i][2]
		except: use=1
		if use :
			avgr.add_image(img)				# only include the particle if we've tagged it as good
			if img.has_attr("source_n") : incl.append(img["source_n"])
#			xforms.append(ptcl_info[i][1])
		elif img.has_attr("source_n") : excl.append(img["source_n"])

	avg=avgr.finish()

	# normalize to the reference, this should make make3dpar work better as we can skip the normalization step
	if ref!=None :
		if setsfref:
			avg.process_inplace("filter.matchto",{"to":ref,"interpolate":0,"keephires":1})
			avg-=avg.get_edge_mean()
		else : avg.process_inplace("normalize.toimage",{"to":ref})

		avg["class_qual"]=avg.cmp("ccc",ref)

	# set some useful attributes
	if len(incl)>0 or len(excl)>0 :
		if len(incl)>0 : avg["class_ptcl_idxs"]=incl
		if len(excl)>0 : avg["exc_class_ptcl_idxs"]=excl
#		if len(xforms)>0: avg["class_ptcl_xforms"]=xforms
		avg["class_ptcl_src"]=img["source_path"]

	return avg

def class_average(images,ref=None,niter=1,normproc=("normalize.edgemean",{}),prefilt=0,align=("rotate_translate_flip",{}),
		aligncmp=("ccc",{}),ralign=None,raligncmp=None,averager=("mean",{}),scmp=("ccc",{}),keep=1.5,keepsig=1,automask=0,saveali=0,verbose=0,callback=None,center="xform.center"):
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
		try: ali["model_id"]=ref["model_id"]
		except: pass
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
		try:
			gmw=max(5,ref["nx"]/16)		# gaussian mask width
			#ref.process_inplace("filter.highpass.gauss",{"cutoff_pixels":min(ref["nx"]/10,5)})	# highpass to reduce gradient issues
			#ref.process_inplace("normalize.circlemean")
			#ref2=ref.process("mask.gaussian",{"inner_radius":ref["nx"]/2-gmw,"outer_radius":gmw/1.3})
			#ref2.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.07})	# highpass to reduce gradient issues
			#ref2.process_inplace("normalize.circlemean")
			#ref2.process_inplace("threshold.binary",{"value":ref["mean"]+ref["sigma"]*1.5})
			#ref2.process_inplace("xform.centerofmass",{"threshold":0.5})						# TODO: should probably check how well this works
			#fxf=ref2["xform.align2d"]
			#ref.translate(fxf.get_trans())
			ref.process_inplace(center)
			ref.process_inplace("normalize.circlemean",{"radius":ref["nx"]/2-gmw})
			ref.process_inplace("mask.gaussian",{"inner_radius":ref["nx"]/2-gmw,"outer_radius":gmw/1.3})
			ref_orient=None
		except:
			traceback.print_exc()
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

		mean,sigma=0.0,1.0		# defaults for when similarity isn't computed

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
			if saveali and it==niter : ali.write_image("aligned.hdf",-1)

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
		ref["class_ptcl_qual"]=mean
		ref["class_ptcl_qual_sigma"]=sigma

		# A little masking before the next iteration
		gmw=max(5,ref["nx"]/12)		# gaussian mask width
		ref.process_inplace("normalize.circlemean",{"radius":ref["nx"]/2-gmw})
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
		elif options.ref and os.path.exists(options.input):
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
