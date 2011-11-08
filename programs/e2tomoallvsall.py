#!/usr/bin/env python

#
# Author: Jesus Galaz (with adapted code from e2classaverage3d), 07/2011
# Copyright (c) 2011 Baylor College of Medicine
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
from pprint import pprint
from EMAN2db import EMTask
from operator import itemgetter	

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <output> [options]

	STILL HEAVILY UNDER DEVELOPMENT.
	This program produces a final average of a dataset (and mutually exclusive classes of a given size [in terms of a minimum # of particles in each class]),
	where all particles have been subjected to all vs all alignments and hierarchical ascendent classification.
	
	See the e2spt Users Guide downloadable in PDF format from the EMAN2 Wiki for an explanation of this procedure.

	Three pre-processing operations are provided: mask, normproc and preprocess. They are executed in that order. Each takes
	a generic <processor>:<parm>=<value>:...  string. While you could provide any valid processor for any of these options, if
	the mask processor does not produce a valid mask, then the default normalization will fail. It is recommended that you
	specify the following, unless you really know what you're doing:
	
	--mask=mask.sharp:outer_radius=<safe radius>
	--preprocess=filter.lowpass.gauss:cutoff_freq=<1/resolution in A>
	
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")

	
	parser.add_argument("--input", type=str, help="The name of the input volume stack. MUST be HDF or BDB, since volume stack support is required.", default=None)
	parser.add_argument("--output", type=str, help="The name of the output class-average stack. MUST be HDF or BDB, since volume stack support is required.", default=None)
		
	parser.add_argument("--iter", type=int, help="The number of iterations to perform. Default is 1.", default=1)
	
	parser.add_argument("--savesteps",action="store_true", help="If set, will save the average after each iteration to class_#.hdf. Each class in a separate file. Appends to existing files.",default=False)
	parser.add_argument("--saveali",action="store_true", help="If set, will save the aligned particle volumes in class_ptcl.hdf. Overwrites existing file.",default=False)
	
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_argument("--preprocess",type=str,help="A processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.",default=None)
	parser.add_argument("--npeakstorefine", type=int, help="The number of best coarse alignments to refine in search of the best final alignment. Default=4.", default=4)
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average. Default is rotate_translate_3d:search=10:delta=15:dphi=15", default="rotate_translate_3d:search=10:delta=15:dphi=15")
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. Default is refine.3d, specify 'None' to disable", default="refine_3d")
	parser.add_argument("--raligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo")
	
	parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean")
	
	parser.add_argument("--postprocess",type=str,help="A processor to be applied to the volume after averaging the raw volumes, before subsequent iterations begin.",default=None)
	
#	parser.add_argument('--reverse_contrast', action="store_true", default=False, help=""" This multiplies the input particles by -1. Remember that EMAN2 **MUST** work with 'white protein' """)
	
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkrefine", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")
#	parser.add_argument("--automask",action="store_true",help="Applies a 3-D automask before centering. Can help with negative stain data, and other cases where centering is poor.")
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if options.align: 
		options.align=parsemodopt(options.align)
	if options.ralign: 
		options.ralign=parsemodopt(options.ralign)
	
	if options.aligncmp: 
		options.aligncmp=parsemodopt(options.aligncmp)
	if options.raligncmp: 
		options.raligncmp=parsemodopt(options.raligncmp)
	
	if options.averager: 
		options.averager=parsemodopt(options.averager)
	if options.normproc: 
		options.normproc=parsemodopt(options.normproc)
	if options.mask: 
		options.mask=parsemodopt(options.mask)
	if options.preprocess: 
		options.preprocess=parsemodopt(options.preprocess)
	if options.postprocess: 
		options.postprocess=parsemodopt(options.postprocess)
		
	if options.path and ("/" in options.path or "#" in options.path) :
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
		sys.exit(1)
		
	if options.path and options.path[:4].lower()!="bdb:": 
		options.path="bdb:"+options.path
	if not options.path: 
		options.path="bdb:"+numbered_path("spt",True)

	hdr = EMData(options.input,0,True)
	nx = hdr["nx"]
	ny = hdr["ny"]
	nz = hdr["nz"]
	if nx!=ny or ny!=nz :
		print "ERROR, input volumes are not cubes"
		sys.exit(1)

	logger = E2init(sys.argv,options.ppid)
	
	# Initialize parallelism if being used
	if options.parallel :
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
		pclist=[options.input]
	
		etc.precache(pclist)
	
	nptcl = EMUtil.get_image_count(options.input)
	if nptcl<3: 
		print "ERROR: at least 3 particles are required in the input stack for all vs all. Otherwise, to align 2 particles (one to the other) use e2classaverage3d.py"
		sys.exit(1)
	
	newptcls=[]
	allptcls=[]
	for i in range(nptcl):
		#a=EMData(options.input,i,True)				#Read the header only for each particle and set the "multiplicity" parameter to 1 if there is no such parameter
		a=EMData(options.input,i)
		if 'spt_multiplicity' not in a.get_attr_dict():
			a['spt_multiplicity']=1
			a.write_image(options.input,i)
		if 'spt_ptcl_indxs' not in a.get_attr_dict():		#Set the spt_ptcl_indxs header parameter to keep track of what particles the current particle is an average of
			a['spt_ptcl_indxs']=i				#In this case, the fresh/new stack should contain particles where this parameter is the particle number itself
			a.write_image(options.input,i)
		newptcls.append(a)
		allptcls.append(a)
		
	#allptcls = range(nptcl)
	#nptcls = range(nptcl)
	
	setOLD = set([])
	oldptcls = []

	for k in range(options.iter):
		'''
		Make ALL vs ALL comparisons among all NEW particle INDEXES.
		NOTE: In the first round all the particles are "new"
		'''
		
		nnew = len(newptcls)
		newptclsmap = list(enumerate([range(i,nnew) for i in range(1,nnew)]))
		
		#file_map = [ (newptcls[x[0]] , [newptcls[z] for z in x[1]] ) for x in enum ]
		#numbers = [i +1 for i in range(total)]
		#num_map = [ (numbers[x(0)]] , [numbers[z] for z in x[1]] ) for x in enum ]
		
		ranks = []
		tasks = []
		
		jj=0					#counter to keep track of the number of comparisons
		for ptcl1, compare in newptclsmap:
			for ptcl2 in compare:
				
				#ref = EMData(options.input,ptcl1)
				
				ref = newptcls[ptcl1]
				
				#task = Align3DTask(ref,["cache",options.input,ptcl2],jj,ptcl1,ptcl2,"Aligning particle#%d VS particle#%d in iteration %d" % (ptcl1,ptcl2,k),options.mask,options.normproc,options.preprocess,
				#options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
				
				task = Align3DTask(ref,["cache",newptcls[ptcl2]],jj,'new_' + str(ptcl1),'new_' +str(ptcl2),"Aligning particle#%d VS particle#%d of the 'NEW SET' in iteration %d" % (ptcl1,ptcl2,k),options.mask,options.normproc,options.preprocess,
				options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
				
				tasks.append(task)
				
				jj+=1
				
				#info = [name1, name2, rt2apply, score]
				#ranks.append(info)							#Add comparison to a list that will shortly be sorted by ccc score, to determine the "UNIQUE best pairs"
		
		'''
		Make comparisons for all NEW VS all OLD particles. "NEW" means particles that didn't exist in the previous round.
		There are no "new" and "old" particles in the first round; thus the loop below is needed only for k>0
		'''
				
		if k > 0:		
			if len(allptcls) == 1:
				print "The all vs all alignment has finalized and converged into one average"
				print "TERMINATING"
				sys.exit()
				
			setNEW = set(newptcls)
			setALL = set(allptcls)
			setOLD = setALL - setNEW
			
			oldptcls = list(setOLD)

			if options.verbose:
				print "The set of NEW particles has these many in it", len(setNEW)
				print "The set of ALL particles has these many in it", len(setALL)
				print "Therefore, the difference is the amount of old particles remaining", len(setOLD)

			xx=0
			for ptcl1 in newptcls:
				yy=0
				for ptcl2 in oldptcls:
					
					#ref = EMData(ptcl1)	
					#task = Align3DTask(ref,["cache",options.input,ptcl2],comparison,ptcl1,ptcl2,"Aligning particle#%d VS particle#%d in iteration %d" % (ptcl1,ptcl2,k),options.mask,options.normproc,options.preprocess,
					#options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
					
					task = Align3DTask(ptcl1,["cache",ptcl2],jj,'new_'+str(xx),'old_'+str(yy),"Aligning particle#%d of the OLD set VS particle#%d of the NEW set, in iteration %d" % (xx,yy,k),options.mask,options.normproc,options.preprocess,
					options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
					
					tasks.append(task)
										
					yy+=1
				xx+=1
				jj+=1	
		
		# start the alignments running
		tids=etc.send_tasks(tasks)
		if options.verbose: 
			print "%d tasks queued in iteration %d"%(len(tids),k) 

		# Wait for alignments to finish and get results
		results = get_results(etc,tids,options.verbose)
		
		results = sorted(results, key=itemgetter('score'))

		if options.verbose>2 : 
			print "The SORTED results are:"
			print results
			
			coeffs=[]
			for i in results:
				coeffs.append(i['score'])

			print "Therefore the ordered coefficients are", coeffs
		
		print "This many comparisons were made", len(results)
		
		tried = set()
		averaged = set()
		averages =[]
		
		mm=0
		for pair in results:
			if pair['ptcl1'] not in tried and pair['ptcl2'] not in tried:
				
				tried.add(pair['ptcl1'])		
				tried.add(pair['ptcl2'])						#Add both of the scanned particles to both groups of scanned/tried AND averaged particles, 
													#since all particle pairs that aren't already in the 'tried' group MUST be averaged
				averaged.add(pair['ptcl1'])
				averaged.add(pair['ptcl2'])
				
				avgr=Averagers.get(averager[0], averager[1])				#Call the averager
				
				ptcl1=EMAN2.EMData()
				
				if 'new' in pair['ptcl1']:
					ptcl1 = newptcls[int(pair['ptcl1'].split('_')[-1])]				
				elif 'old' in pair['ptcl1']:
					ptcl1 = oldptcls[int(pair['ptcl1'].split('_')[-1])]
				
				ptcl1 = ptcl1 * ptcl1['spt_multiplicity']				#Take the multiplicity of ptcl1 into account
					
				avgr.add_image(ptcl1)							#Add particle 1 to the average
				
				ptcl2=EMAN2.EMData()
				
				if 'new' in pair['ptcl2']:
					ptcl2 = newptcls[int(pair['ptcl2'].split('_')[-1])]				
				elif 'old' in pair['ptcl2']:
					ptcl2 = oldptcls[int(pair['ptcl2'].split('_')[-1])]
				
				ptcl2 = ptcl2 * ptcl2['spt_multiplicity']				#Take the multiplicity of ptcl1 into account				
				
				ptcl2.process_inplace("xform",{"transform":pair["xform.align3d"]})	#Apply the relative alignment between particles 1 and 2 to particle 2,
													#since particle 1 is always the "fixed" one and particle 2 the "moving" one during alignment
				
				avgr.add_image(ptcl2)							#Add the transformed (rotated and translated) particle 2 to the average
		
				avg=avgr.finish()
				avg["spt_ptcl_indxs"] = list(set(ptcl1["spt_ptcl_indxs"]).union(set(ptcl2["spt_ptcl_indxs"])))	#Keep track of what particles go into each average or "new particle"				
				avg["spt_ptcl_src"] = options.input
				
				avg['origin_x'] = 0							#The origin needs to be reset to ZERO to avoid display issues in Chimera
				avg['origin_y'] = 0
				avg['origin_z'] = 0
				
				#avg.write_image("bdb:new%02d_ptcl"%(k),mm)
				avg.write_image("%s/new%02d_ptcl"%(options.path,k),mm)			#Particles from a "new round" need to be in a "new stack" defined by counter k; the number
													#of particles in it is determined by counter mm, which increments when a pair is averaged
				averages.append(avg)	   #The list of averages will become the new set of "newptcls"
				
				mm+=1
				
			if pair['ptcl1'] not in tried:
				tried.add(pair['ptcl1'])
				
			if pair['ptcl2'] not in tried:
				tried.add(pair['ptcl2'])
						
		for particle in averaged:
			if particle in newptcls:
				newptcls.remove(particle)				
			if particle in oldptcls:
				oldptcls.remove(particle)
				
		oldptcls.append(newptcls)			#All the particles from the newptcls list that were not averaged become "old"
		newptcls=averages				#All the new averages become part of the new "newptcls" list
	
		E2end(logger)

		return()


'''
def make_averages(ptcl_file,align_parms,averager,saveali,verbose=1):
	"""Will take a set of alignments and an input particle stack filename and produce a new class-average.
	Particles may be excluded based on the keep and keepsig parameters. If keepsig is not set, then keep represents
	an absolute fraction of particles to keep (0-1). Otherwise it represents a sigma multiplier akin to e2classaverage.py"""
	
	avgr=Averagers.get(averager[0], averager[1])

	for i,ptcl_parms in enumerate(align_parms):
		ptcl=EMData(ptcl_file,i)
		ptcl.process_inplace("xform",{"transform":ptcl_parms[0]["xform.align3d"]})
		
		avgr.add_image(ptcl)
		
		if saveali:
			ptcl['origin_x'] = 0
			ptcl['origin_y'] = 0		# jesus - the origin needs to be reset to ZERO to avoid display issues in Chimera
			ptcl['origin_z'] = 0
			ptcl.write_image("bdb:class_ptcl",i)
	
	ret=avgr.finish()
	ret["class_ptcl_idxs"]=range(align_parms)
	ret["class_ptcl_src"]=ptcl_file
	
	return ret
'''





def get_results(etc,tids,verbose):
	"""This will get results for a list of submitted tasks. Won't return until it has all requested results.
	aside from the use of options["ptcl"] this is fairly generalizable code. """
	
	# wait for them to finish and get the results
	# results for each will just be a list of (qual,Transform) pairs
	results=[0]*len(tids)		# storage for results
	ncomplete=0
	tidsleft=tids[:]
	while 1:
		time.sleep(5)
		proglist=etc.check_task(tidsleft)
		nwait=0
		for i,prog in enumerate(proglist):
			if prog==-1 : nwait+=1
			if prog==100 :
				r=etc.get_results(tidsleft[i])						# results for a completed task
				#print "\n@@@@@@The results for the completed task are", r
				comparison=r[0].options["comparison"]					# get the comparison number from the task rather than trying to work back to it
				results[comparison]=r[1]["final"][0]					# this will be a list of (qual,Transform), containing the BEST peak ONLY
				
				results[comparison]['ptcl1']=r[0].options['ptcl1']			#Associate the result with the pair of particles involved
				results[comparison]['ptcl2']=r[0].options['ptcl2']

				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if verbose:
			print "  %d tasks, %d complete, %d waiting to start        \r"%(len(tids),ncomplete,nwait)
			sys.stdout.flush()
	
		if len(tidsleft)==0: 
			break
		
	return results


class Align3DTask(EMTask):
	"""This is a task object for the parallelism system. It is responsible for aligning one 3-D volume to another, with a variety of options"""

	def __init__(self,fixedimage,image,comparison,ptcl1,ptcl2,label,mask,normproc,preprocess,npeakstorefine,align,aligncmp,ralign,raligncmp,shrink,shrinkrefine,verbose):
		"""fixedimage and image may be actual EMData objects, or ["cache",path,number]
	label is a descriptive string, not actually used in processing
	ptcl is not used in executing the task, but is for reference
	other parameters match command-line options from e2classaverage3d.py
	Rather than being a string specifying an aligner, 'align' may be passed in as a Transform object, representing a starting orientation for refinement"""
		data={}
		data={"fixedimage":fixedimage,"image":image}
		EMTask.__init__(self,"ClassAv3d",data,{},"")

		self.options={"comparison":comparison,"ptcl1":ptcl1,"ptcl2":ptcl2,"label":label,"mask":mask,"normproc":normproc,"preprocess":preprocess,"npeakstorefine":npeakstorefine,"align":align,"aligncmp":aligncmp,"ralign":ralign,"raligncmp":raligncmp,"shrink":shrink,"shrinkrefine":shrinkrefine,"verbose":verbose}
	
		#self.options={"comparison":comparison,"label":label,"mask":mask,"normproc":normproc,"preprocess":preprocess,"npeakstorefine":npeakstorefine,"align":align,"aligncmp":aligncmp,"ralign":ralign,"raligncmp":raligncmp,"shrink":shrink,"shrinkrefine":shrinkrefine,"verbose":verbose}

	def execute(self,callback=None):
		"""This aligns one volume to a reference and returns the alignment parameters"""
		options=self.options
		if options["verbose"]>1: 
			print "Aligning ",options["label"]

		if isinstance(self.data["fixedimage"],EMData):
			fixedimage=self.data["fixedimage"]
		else: 
			print "You are not passing in an EMData REFERENCE!"

			#fixedimage=EMData(self.data["fixedimage"][1],self.data["fixedimage"][2])
		
		if isinstance(self.data["image"],EMData):
			image=self.data["image"]
		else: 
			print "You are not passing in an EMData PARTICLE!"
			
			#image=EMData(self.data["image"][1],self.data["image"][2])
		
		# Preprocessing applied to both volumes.
		# Make the mask first, use it to normalize (optionally), then apply it 
		mask=EMData(image["nx"],image["ny"],image["nz"])
		mask.to_one()
		if options["mask"] != None:
			print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options["mask"][0],options["mask"][1]) 
			mask.process_inplace(options["mask"][0],options["mask"][1])
		
		# normalize
		if options["normproc"] != None:
			if options["normproc"][0]=="normalize.mask" : options["normproc"][1]["mask"]=mask
			fixedimage.process_inplace(options["normproc"][0],options["normproc"][1])
			image.process_inplace(options["normproc"][0],options["normproc"][1])
		
		fixedimage.mult(mask)
		image.mult(mask)
		
		# preprocess
		if options["preprocess"] != None:
			fixedimage.process_inplace(options["preprocess"][0],options["preprocess"][1])
			image.process_inplace(options["preprocess"][0],options["preprocess"][1])
		
		# Shrinking both for initial alignment and reference
		if options["shrink"]!=None and options["shrink"]>1 :
			sfixedimage=fixedimage.process("math.meanshrink",{"n":options["shrink"]})
			simage=image.process("math.meanshrink",{"n":options["shrink"]})
		else :
			sfixedimage=fixedimage
			simage=image
			
		if options["shrinkrefine"]!=None and options["shrinkrefine"]>1 :
			if options["shrinkrefine"]==options["shrink"] :
				s2fixedimage=sfixedimage
				s2image=simage
			else :
				s2fixedimage=fixedimage.process("math.meanshrink",{"n":options["shrinkrefine"]})
				s2image=image.process("math.meanshrink",{"n":options["shrinkrefine"]})
		else :
			s2fixedimage=fixedimage
			s2image=image
			

		if options["verbose"]>1: 
			print "Align size %d,  Refine Align size %d"%(sfixedimage["nx"],s2fixedimage["nx"])

		#If a Transform was passed in, we skip coarse alignment
		if isinstance(options["align"],Transform):
			bestcoarse=[{"score":1.0,"xform.align3d":options["align"]}]
			if options["shrinkrefine"]>1: 
				bestcoarse[0]["xform.align3d"].set_trans(bestcoarse[0]["xform.align3d"].get_trans()/float(options["shrinkrefine"]))
		
		#This is the default behavior, seed orientations come from coarse alignment
		else:
			# returns an ordered vector of Dicts of length options.npeakstorefine. The Dicts in the vector have keys "score" and "xform.align3d"
			bestcoarse=simage.xform_align_nbest(options["align"][0],sfixedimage,options["align"][1],options["npeakstorefine"],options["aligncmp"][0],options["aligncmp"][1])
			scaletrans=options["shrink"]/float(options["shrinkrefine"])
			if scaletrans!=1.0:
				for c in bestcoarse:
					c["xform.align3d"].set_trans(c["xform.align3d"].get_trans()*scaletrans)

		# verbose printout
		if options["verbose"]>1 :
			for i,j in enumerate(bestcoarse): print "coarse %d. %1.5g\t%s"%(i,j["score"],str(j["xform.align3d"]))

		if options["ralign"]!=None :
			# Now loop over the individual peaks and refine each
			bestfinal=[]
			for bc in bestcoarse:
				options["ralign"][1]["xform.align3d"]=bc["xform.align3d"]
				ali=s2image.align(options["ralign"][0],s2fixedimage,options["ralign"][1],options["raligncmp"][0],options["raligncmp"][1])
				
				try: 
					bestfinal.append({"score":ali["score"],"xform.align3d":ali["xform.align3d"],"coarse":bc})
				except:
					bestfinal.append({"xform.align3d":bc["xform.align3d"],"score":1.0e10,"coarse":bc})

			if options["shrinkrefine"]>1 :
				for c in bestfinal:
					c["xform.align3d"].set_trans(c["xform.align3d"].get_trans()*float(options["shrinkrefine"]))

			# verbose printout of fine refinement
			if options["verbose"]>1 :
				for i,j in enumerate(bestfinal): 
					print "fine %d. %1.5g\t%s"%(i,j["score"],str(j["xform.align3d"]))

		else: 
			bestfinal=bestcoarse
		
		#If you just sort 'bestfinal' it will be sorted based on the 'coarse' key in the dictionaries of the list
		#because they come before the 'score' key of the dictionary (alphabetically)
		
		bestfinal = sorted(bestfinal, key=itemgetter('score'))
		
		#print "\n$$$$\n$$$$\n$$$$\n$$$$\n$$$$\n$$$$The best peaks sorted are" 
		#
		#for i in bestfinal:
		#	print bestfinal
		
		if bestfinal[0]["score"] == 1.0e10 :
			print "Error: all refine alignments failed for %s. May need to consider altering filter/shrink parameters. Using coarse alignment, but results are likely invalid."%self.options["label"]
		
		if options["verbose"]: 
			print "Best %1.5g\t %s"%(bestfinal[0]["score"],str(bestfinal[0]["xform.align3d"]))

		if options["verbose"]: 
			print "Done aligning ",options["label"]
		
		return {"final":bestfinal,"coarse":bestcoarse}

'''
def allvsall(parameters,raw_stack,stack=''):
	
	
	all_names = []
	for i in range(n):							#Writing individual raw-particles to file should only happen if a stack is provided
		num = str(i+1).zfill(len(str(n)))
		name = parameters['ID'] + '_rt' + str(0).zfill(len(str(parameters['rounds']))) + '_#' + num + '_m1.hdf'
		raw_stack[i].write_image(name,0)	
		all_names.append(name)

	new = '_rt' + str(0).zfill(len(str(parameters['rounds']))) + '_'	#Mark particles according to round they were produced in
	
	
	ranks = []

	for k in range(parameters['rounds']):
		
		if parameters['verbose']:
			print "I am on round %d of ALL vs ALL alignment", k

				
		common = parameters['ID'] 					#Mark "common" to ALL particles no matter what round they were produced in
		
		current = os.getcwd()
		findir = os.listdir(current)
		
		new_names = []
		if k == 0:
			new_names = list(all_names) 	#All the particles are "new" in the first round
			 
		else:
			for f in findir:
				if new in f:
					new_names.append(f)
			new_names.sort()
		

		
		
		#Make ALL vs ALL comparisons among all NEW particles.
		#NOTE: In the first round all the particles are "new"
		
		total = len(new_names)
		enum = list(enumerate([range(i,total) for i in range(1,total)]))
		file_map = [ (new_names[x[0]] , [new_names[z] for z in x[1]] ) for x in enum ]

		#numbers = [i +1 for i in range(total)]
		#num_map = [ (numbers[x(0)]] , [numbers[z] for z in x[1]] ) for x in enum ]

		for name1, compare in file_map:
			for name2 in compare:
				if parameters['verbose']:
						print "I will align %s VS %s" % (name1,name2)
				
				particle1 = EMData(name1,0)		#Load the particles
				particle2 = EMData(name2,0)

				particle1_ed = editor(parameters, particle1,'no')[0]	#Preprocess the particles to prepare them for alignment
				particle2_ed = editor(parameters, particle2,'no')[0]

				rtresults = tomohunt(parameters, particle1, particle2)		#Actual alignment; get back only the returned transform
				rt2apply = rtresults['xform.align3d']
				score = rtresults['score']
				
				info = [name1, name2, rt2apply, score]
				ranks.append(info)							#Add comparison to a list that will shortly be sorted by ccc score, to determine the "UNIQUE best pairs"
				
		
		#Make comparisons for all NEW VS all OLD particles
		#"NEW" means particles that didn't exist in the previous round
		
		if k > 0:
			all_names = [] 			#Reset particles to read the new set fresh on the next round
			for f in findir:
				if common in f and '.hdf' in f and "#" in f:
					all_names.append(f)
			
			if len(all_names) == 1:
				print "The all vs all alignment has finalized and converged into one average"
				print "TERMINATING"

				#launch_childprocess('mkdir ' + parameters['ID'])	
				launch_childprocess('mv *' + parameters['ID'] + '* ' + parameters['ID'])
				sys.exit()
				
			all_names.sort()

			setNEW = set(new_names)
			setALL = set(all_names)

			if parameters['verbose']:
				print "The set of NEW particles is", setNEW
				print "The set of ALL particles is", setALL
				print "What is common to all files names is", common
				print "And the label unique to this new round is", new
				old_names = setALL - setNEW
				print "Therefore, the set of all OLD particles is their difference", old_names
				old_names = sorted(old_names)

			for np in new_names:
				for op in old_names:
					if parameters['verbose']:
						print "I will align %s VS %s" % (np,op)

					
					particle1 = EMData(np,0)		#Load the particles
					particle2 = EMData(op,0)

					particle1_ed = editor(parameters, particle1,'no')[0]		#Preprocess the particles to prepare them for alignment
					particle2_ed = editor(parameters, particle2,'no')[0]

					rtresults = tomohunt(parameters, particle1, particle2)		#Actual alignment; get back only the returned transform
					rt2apply = rtresults['xform.align3d']
					score = rtresults['score']					
					info = [np, op, rt2apply, score]
					ranks.append(info)						#Add comparison to a list that will shortly be sorted by ccc score
		
		ranks_sorted = sorted(ranks, key=operator.itemgetter(-1), reverse=True)		#Sort pairs by ccc score
		
		ranks_filename = 'ranks_' + common + new.rstrip('_') + '.lst'			#Store the ordered rank of ALL PAIRS to a file
		ranks_file = open(ranks_filename,'w')
		
		#
		#Determine WHICH PAIRS TO AVERAGE
		#	
		
		best_pairs_filename = 'best_pairs_' + common + new.rstrip('_')  + '.lst'	#Store AVERAGED PAIRS in the round for fast reference/verification
		best_pairs_file = open(best_pairs_filename,'w')
		
		done_dir = 'done_' + common + new.rstrip('_')
		launch_childprocess('mkdir ' + done_dir)
					
		new = '_rt' + str(k + 1).zfill(len(str(parameters['rounds']))) + '_'	#Mark particles according to round they were produced in
		jj = 1
		tried = set()
		averaged = set()
		ranks = []
		for rs in ranks_sorted:
			rotations = rs[2].get_rotation()
			translations = rs[2].get_trans()
			az = str(int(rotations['az']*100)/100.0).zfill(6)
			alt = str(int(rotations['alt']*100)/100.0).zfill(6)
			phi = str(int(rotations['phi']*100)/100.0).zfill(6)
			x = str(int(translations[0]*100)/100.0).zfill(6)
			y = str(int(translations[1]*100)/100.0).zfill(6)
			z = str(int(translations[2]*100)/100.0).zfill(6)
			score = str(int(rs[3]*1000000)/1000000.0).zfill(6)
			
			pair = rs[0] + ' vs ' + rs[1]
			line =  pair + ' az=%s alt=%s phi=%s tx=%s ty=%s tz=%s score=%s\n' %(az,alt,phi,x,y,z,score)
			ranks_file.write(line)
			
			file1 = rs[0]
			file2 = rs[1]
			if file1 not in tried and file2 not in tried:
				a = EMData(file1,0)
				b = EMData(file2,0)
				print "\n\n\n===========File 1 and 2 are", file1, file2
				am = int(file1.split('.hdf')[0].split('_m')[1])		#To average appropriately you need to account for the "multiplicity"
				bm = int(file2.split('.hdf')[0].split('_m')[1])		#of the particles being averaged. That is, consider how many particles
				multiplicity = am + bm					#they themselves are an average of. The filenames' "m" tags determines this.

				average = (a + b)/multiplicity
				average = average.process(parameters['normalization'])
				
				avgnum = str(jj).zfill(len(str(n)))
				avgname = common + new + '#' + avgnum + '_m' + str(multiplicity) + '.hdf'
				a.write_image(avgname,0)
				jj += 1
				
				best_pairs_file.write(pair + ' ' + score + '\n')
				print "In round", k
				print "I have averaged %s and %s into %s" % (file1, file2,avgname)
				averaged.add(file1)
				averaged.add(file2)
			if file1 not in tried:
				tried.add(file1)
			if file2 not in tried:
				tried.add(file2)
				
		for qq in averaged:
			launch_childprocess('mv ' + qq + ' ' + done_dir)
		launch_childprocess('mv *.lst ' + done_dir)

		for pp in ranks_sorted:
			pps = set(pp)
			dif = pps - averaged
			if len(dif) == 4:
				ranks.append(pp)
		k += 1
	return()
'''

if __name__ == '__main__':
	main()



