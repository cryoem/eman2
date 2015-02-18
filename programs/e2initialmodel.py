#!/usr/bin/env python

#
# Author: Steven Ludtke, 11/30/2008 (ludtke@bcm.edu)
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
import random
from math import *
import os
import sys
from e2simmx import cmponetomany
from EMAN2jsondb import JSTask,jsonclasses

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]
	This program will take a set of reference-free class-averages (or other projections) and generate a set of possible
	3-D initial models. It does this by heavily downsampling the data, then running a number of very fast, full iterative
	refinements, each seeded with a random starting model. The results are sorted in order of apparent agreement with the
	data, such that at the end, the first numbered model should be the best result. Ideally the top few answers will all
	qualtitatively agree on the overall structure. If they do not, the results should be thoroughly assessed manually to
	insure a sensible result. By default this routine will generate 10 initial models, but this may be fewer or more than
	is strictly necessary depending on a number of factors. If the data is highly structurally heterogeneous, particularly
	if combined with a strongly preferred orientation, a correct solution using this technique may not be possible, but
	for most situations it will work well. For other situations, single particle tomography presents a good alternative
	for generating initial models."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_header(name="initialmodelheader", help='Options below this label are specific to e2initialmodel', title="### e2initialmodel options ###", row=1, col=0, rowspan=1, colspan=3)
	parser.add_argument("--input", dest="input", default=None,type=str, help="This file should contain good class-averages to use in constructing the initial model", browser='EMBrowserWidget(withmodal=True,multiselect=False)', guitype='filebox', row=0, col=0, rowspan=1, colspan=3)
	parser.add_argument("--iter", type = int, default=8, help = "The total number of refinement iterations to perform, typically 5-10", guitype='intbox', row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--tries", type=int, default=10, help="The number of different initial models to generate in search of a good one", guitype='intbox', row=2, col=1, rowspan=1, colspan=1)
	parser.add_argument("--shrink", dest="shrink", type = int, default=0, help="Optionally shrink the input particles by an integer factor prior to reconstruction. Default=0, no shrinking", guitype='shrinkbox', row=2, col=2, rowspan=1, colspan=1)
	parser.add_argument("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos",default="c1", guitype='symbox', row=4, col=0, rowspan=1, colspan=2)
	parser.add_argument("--maskproc", default=None, type=str,help="Default=none. If specified, this mask will be performed after the built-in automask, eg - mask.soft to remove the core of a virus", )
#	parser.add_argument("--savemore",action="store_true",help="Will cause intermediate results to be written to flat files",default=False, guitype='boolbox', expert=True, row=5, col=0, rowspan=1, colspan=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--orientgen",type=str, default="eman:delta=9.0:inc_mirror=0:perturb=1",help="The type of orientation generator. Default is eman:delta=9.0:inc_mirror=0:perturb=1. See e2help.py orientgens", guitype='strbox', expert=True, row=4, col=2, rowspan=1, colspan=1)
	parser.add_argument("--parallel","-P",type=str,help="Run in parallel, specify type:<option>=<value>:<option>=<value>. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel",default="thread:1", guitype='strbox', row=6, col=0, rowspan=1, colspan=2)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	# Database Metadata storage
	#parser.add_argument("--dbls",type=str,default=None,help="data base list storage, used by the workflow. You can ignore this argument.")

	(options, args) = parser.parse_args()
	verbose=options.verbose

	try: ptcls=EMData.read_images(options.input)
	except:
		print "Error: bad input file"
		exit(1)
	apix=ptcls[0]["apix_x"]
	if options.shrink>1 : apix*=options.shrink

	for i in range(len(ptcls)):
		ptcls[i].process_inplace("normalize.edgemean",{})
		if options.shrink>1 :
			ptcls[i]=ptcls[i].process("math.meanshrink",{"n":options.shrink})
	if ptcls[0]["nx"]>160 : print "WARNING: using a large box size may be slow. Suggest trying --shrink="
	if not ptcls or len(ptcls)==0 : parser.error("Bad input file")
	boxsize=ptcls[0].get_xsize()
	if verbose>0 : print "%d particles %dx%d"%(len(ptcls),boxsize,boxsize)
	print "Models will be %1.3f A/pix"%apix

	[og_name,og_args] = parsemodopt(options.orientgen)

	try:
			sfcurve=XYData()
			sfcurve.read_file("strucfac.txt")

			sfcurve.update()
	except : sfcurve=None

	if options.maskproc!=None :
		mask2=EMData(boxsize,boxsize,boxsize)
		mask2.to_one()
		parms=parsemodopt(options.maskproc)
		if parms[0]=="mask.auto3d":
			print "Error, maskproc may not be mask.auto3d, it must be a processor that does not rely on the input map density to function"
			sys.exit(1)
		mask2.process_inplace(parms[0],parms[1])
	else: mask2=None

	# angles to use for refinement
	sym_object = parsesym(options.sym)
	orts = sym_object.gen_orientations(og_name,og_args)

	logid=E2init(sys.argv,options.ppid)
	results=[]

	try: os.mkdir("initial_models")
	except: pass
	iters=[int(i[10:12]) for i in os.listdir("initial_models") if i[:10]=="particles_"]
	try : newiter=max(iters)+1
	except : newiter=0
	results_name="initial_models/model_%02d"%newiter
	particles_name="initial_models/particles_%02d.hdf"%newiter

	# we write the pre-processed "particles" (usually class-averages) to disk, both as a record and to prevent collisions
	for i,p in enumerate(ptcls):
		p.write_image(particles_name,i)

	# parallelism
	from EMAN2PAR import EMTaskCustomer			# we need to put this here to avoid a circular reference

	etc=EMTaskCustomer(options.parallel)
	pclist=[particles_name]

	etc.precache(pclist)		# make sure the input particles are precached on the compute nodes

	tasks=[]
	for t in xrange(options.tries):
		tasks.append(InitMdlTask(particles_name,len(ptcls),orts,t,sfcurve,options.iter,options.sym,mask2,options.verbose))

	taskids=etc.send_tasks(tasks)
	alltaskids=taskids[:]			# we keep a copy for monitoring progress

	# This loop runs until all subtasks are complete (via the parallelism system
	ltime=0
	while len(taskids)>0 :
		time.sleep(0.1)
		curstat=etc.check_task(taskids)			# a list of the progress on each task
		if options.verbose>1 :
			if time.time()-ltime>1 :
				print "progress: ",curstat
				ltime=time.time()
		for i,j in enumerate(curstat):
			if j==100 :
				rslt=etc.get_results(taskids[i])		# read the results back from a completed task as a one item dict
				results.append(rslt[1]["result"])
				if options.verbose==1 : print "Task {} ({}) complete".format(i,taskids[i])

		# filter out completed tasks. We can't do this until after the previous loop completes
		taskids=[taskids[i] for i in xrange(len(taskids)) if curstat[i]!=100]


	# Write out the final results
	results.sort()
	for i,j in enumerate(results):
		out_name = results_name+"_%02d.hdf"%(i+1)
		j[1].write_image(out_name,0)
		j[4].write_image(results_name+"_%02d_init.hdf"%(i+1),0)
		print out_name,j[1]["quality"],j[0],j[1]["apix_x"]
		for k,l in enumerate(j[3]): l[0].write_image(results_name+"_%02d_proj.hdf"%(i+1),k)	# set of projection images
		for k,l in enumerate(j[2]):
			l.process("normalize").write_image(results_name+"_%02d_aptcl.hdf"%(i+1),k*2)						# set of aligned particles
			j[3][l["match_n"]][0].process("normalize").write_image(results_name+"_%02d_aptcl.hdf"%(i+1),k*2+1)	# set of projections matching aligned particles


	E2end(logid)

class InitMdlTask(JSTask):

	def __init__(self,ptclfile=None,ptcln=0,orts=[],tryid=0,strucfac=None,niter=5,sym="c1",mask2=None,verbose=0) :
		data={"images":["cache",ptclfile,(0,ptcln)],"strucfac":strucfac,"orts":orts,"mask2":mask2}
		JSTask.__init__(self,"InitMdl",data,{"tryid":tryid,"iter":niter,"sym":sym,"verbose":verbose},"")


	def execute(self,progress=None):
		sfcurve=self.data["strucfac"]
		ptcls=EMData.read_images(self.data["images"][1])
		orts=self.data["orts"]
		options=self.options
		verbose=options["verbose"]
		boxsize=ptcls[0].get_xsize()
		apix=ptcls[0]["apix_x"]
		mask2=self.data["mask2"]

		# We make one new reconstruction for each loop of t
		threed=[make_random_map(boxsize,sfcurve)]		# initial model
		apply_sym(threed[0],options["sym"])		# with the correct symmetry

		# This is the refinement loop
		for it in range(options["iter"]):
			if progress != None: progress(it*100/options["iter"])
			if verbose>0 : print "Iteration %d"%it
#			if options.savemore : threed[it].write_image("imdl.%02d.%02d.mrc"%(t,it))
			projs=[(threed[it].project("standard",ort),None) for ort in orts]		# projections
			for i in projs : i[0].process_inplace("normalize.edgemean")
			if verbose>2: print "%d projections"%len(projs)

			# determine particle orientation
			bss=0.0
			bslst=[]
			for i in range(len(ptcls)):
				sim=cmponetomany(projs,ptcls[i],align=("rotate_translate_flip",{"maxshift":boxsize/5}),alicmp=("ccc",{}),cmp=("frc",{"maxres":20}))
				bs=min(sim)
#				print bs[0]
				bss+=bs[0]
				bslst.append((bs[0],i))
				if verbose>2 : print "align %d \t(%1.3f)\t%1.3g"%(i,bs[0],bss)
				n=sim.index(bs)
				ptcls[i]["match_n"]=n
				ptcls[i]["match_qual"]=bs[0]
				ptcls[i]["xform.projection"]=orts[n]	# best orientation set in the original particle

			bslst.sort()					# sorted list of all particle qualities
			bslst.reverse()
			aptcls=[]
#			for i in range(len(ptcls)*3/4):		# We used to include 3/4 of the particles
			for i in range(len(ptcls)*7/8):
				n=ptcls[bslst[i][1]]["match_n"]
				aptcls.append(ptcls[bslst[i][1]].align("rotate_translate_flip",projs[n][0],{},"ccc",{}))
				if it<2 : aptcls[-1].process_inplace("xform.centerofmass",{})

			bss/=len(ptcls)

			# 3-D reconstruction
			pad=(boxsize*3/2)
			pad-=pad%8
			recon=Reconstructors.get("fourier", {"sym":options["sym"],"size":[pad,pad,pad]})

			# insert slices into initial volume
			recon.setup()
			for p in aptcls:
				p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
				p3=recon.preprocess_slice(p2,p["xform.projection"])
				recon.insert_slice(p3,p["xform.projection"],p.get_attr_default("ptcl_repr",1.0))

			# check particle qualities
			quals=[]
			for i,p in enumerate(aptcls):
				p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
				p3=recon.preprocess_slice(p2,p["xform.projection"])
				recon.determine_slice_agreement(p3,p["xform.projection"],p.get_attr_default("ptcl_repr",1.0),True)
				p.mult(p3["reconstruct_norm"])
				p["reconstruct_absqual"]=p3["reconstruct_absqual"]
				quals.append(p3["reconstruct_absqual"])

			quals.sort()
			qual_cutoff=quals[len(quals)/8]		#not using this right now
			qual=sum(quals)

			mdl=recon.finish(True)

			# insert normalized slices
			recon.setup()
			for p in aptcls:
				p2=p.get_clip(Region(-(pad-boxsize)/2,-(pad-boxsize)/2,pad,pad))
				p3=recon.preprocess_slice(p2,p["xform.projection"])
				recon.insert_slice(p3,p["xform.projection"],p.get_attr_default("ptcl_repr",1.0))

			threed.append(recon.finish(True))

			if verbose>1 : print "Iter %d \t %1.4g (%1.4g)"%(it,bss,qual)

			threed[-1].process_inplace("xform.centerofmass")
			threed[-1]=threed[-1].get_clip(Region((pad-boxsize)/2,(pad-boxsize)/2,(pad-boxsize)/2,boxsize,boxsize,boxsize))
			threed[-1].process_inplace("mask.gaussian",{"inner_radius":boxsize/3.0,"outer_radius":boxsize/12.0})
			threed[-1].process_inplace("filter.lowpass.gauss",{"cutoff_abs":.2})
			threed[-1].process_inplace("normalize.edgemean")
			if it>1 : threed[-1].process_inplace("mask.auto3d",{"radius":boxsize/6,"threshold":threed[-1]["sigma_nonzero"]*.85,"nmaxseed":30,"nshells":boxsize/20,"nshellsgauss":boxsize/20})
			if mask2!=None:threed[-1].mult(mask2)
			threed[-1]["apix_x"]=apix
			threed[-1]["apix_y"]=apix
			threed[-1]["apix_z"]=apix
#			threed[-1]["quality"]=bss
			threed[-1]["quality"]=qual

#			threed[-1]["quality_projmatch"]=
#			display(threed[0])
#			display(threed[-1])
			#debugging output
			#for i in range(len(aptcls)):
				#projs[aptcls[i]["match_n"]].write_image("x.%d.hed"%t,i*2)
				#aptcls[i].write_image("x.%d.hed"%t,i*2+1)
		#display(threed[-1])
		#threed[-1].write_image("x.mrc")
		if verbose>0 : print "Model %d complete. Quality = %1.4f (%1.4f)"%(options["tryid"],bss,qual)

#		if progress!=None : progress(100)		# this should be done automatically when we return
		return {"result":(bss,threed[-1],aptcls,projs,threed[0])}



jsonclasses["InitMdlTask"]=InitMdlTask.from_jsondict

def make_random_map(boxsize,sfcurve=None):
	"""This will make a map consisting of random noise, low-pass filtered and center-weighted for use
	as a random starting model in initial model generation. Note that the mask is eliptical and has random aspect."""

	ret=EMData(boxsize,boxsize,boxsize)
	ret.process_inplace("testimage.noise.gauss",{"mean":0.02,"sigma":1.0})
	#if sfcurve!=None:
		#ret.process_inplace("filter.setstrucfac",{"strucfac":sfcurve})
	ret.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
	ret.process_inplace("xform.centerofmass",{})
#	ret.process_inplace("mask.gaussian",{"inner_radius":boxsize/3.0,"outer_radius":boxsize/12.0})
	ret.process_inplace("mask.gaussian.nonuniform",{"radius_x":boxsize/random.uniform(2.0,5.0),"radius_y":boxsize/random.uniform(2.0,5.0),"radius_z":boxsize/random.uniform(2.0,5.0)})

	return ret

def apply_sym(data,sym):
	"""applies a symmetry to a 3-D volume in-place"""
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(sym)
	ref=data.copy()
	for i in range(1,nsym):
		dc=ref.copy()
		dc.transform(xf.get_sym(sym,i))
		data.add(dc)
	data.mult(1.0/nsym)


if __name__ == "__main__":
    main()

