#!/usr/bin/env python

#
# Author: Jesus Galaz, 07/2011
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
from optparse import OptionParser
import math
from copy import deepcopy
import os
import sys
import random
from random import choice
from pprint import pprint
from EMAN2db import EMTask

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog <output> [options]

	This program produces a final average of a dataset (and mutually exclusive classes of a given size),
	where all particles have been subjected to all vs all alignments and hierarchical ascendent classification.

	Three preprocessing operations are provided for, mask, normproc and preprocess. They are executed in that order. Each takes
	a generic <processor>:<parm>=<value>:...  string. While you could provide any valid processor for any of these options, if
	the mask processor does not produce a valid mask, then the default normalization will fail. It is recommended that you
	specify the following, unless you really know what you're doing:
	
	--mask=mask.sharp:outer_radius=<safe radius>
	--preprocess=filter.lowpass.gauss:cutoff_freq=<1/resolution in A>
	
	"""
			
	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	parser.add_option("--input", type="string", help="The name of the input volume stack. MUST be HDF or BDB, since volume stack support is required.", default=None)
	parser.add_option("--output", type="string", help="The name of the output class-average stack. MUST be HDF or BDB, since volume stack support is required.", default=None)
		
	parser.add_option("--iter", type="int", help="The number of iterations to perform. Default is 1.", default=1)
	
	parser.add_option("--savesteps",action="store_true", help="If set, will save the average after each iteration to class_#.hdf. Each class in a separate file. Appends to existing files.",default=False)
	parser.add_option("--saveali",action="store_true", help="If set, will save the aligned particle volumes in class_ptcl.hdf. Overwrites existing file.",default=False)
	
	parser.add_option("--mask",type="string",help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_option("--normproc",type="string",help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_option("--preprocess",type="string",help="A processor (as in e2proc3d.py) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.",default=None)
	parser.add_option("--ncoarse", type="int", help="Deprecated", default=None)
	parser.add_option("--npeakstorefine", type="int", help="The number of best coarse alignments to refine in search of the best final alignment. Default=4.", default=4)
	parser.add_option("--align",type="string",help="This is the aligner used to align particles to the previous class average. Default is rotate_translate_3d:search=10:delta=15:dphi=15", default="rotate_translate_3d:search=10:delta=15:dphi=15")
	parser.add_option("--aligncmp",type="string",help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. Default is refine.3d, specify 'None' to disable", default="refine_3d")
	parser.add_option("--raligncmp",type="string",help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo")
	
	parser.add_option("--averager",type="string",help="The type of averager used to produce the class average. Default=mean",default="mean")
	
	parser.add_option("--postprocess",type="string",help="A processor to be applied to the volume after averaging the raw volumes, before subsequent iterations begin.",default=None)
	
#	parser.add_option('--reverse_contrast', action="store_true", default=False, help=""" This multiplies the input particles by -1. Remember that EMAN2 **MUST** work with 'white protein' """)
	
	parser.add_option("--shrink", type="int",default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_option("--shrinkrefine", type="int",default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")
#	parser.add_option("--automask",action="store_true",help="Applies a 3-D automask before centering. Can help with negative stain data, and other cases where centering is poor.")
	parser.add_option("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n",type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")


	logger = E2init(sys.argv)
	
	# Initialize parallelism if being used
	if options.parallel :
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
		pclist=[options.input]
	
		etc.precache(pclist)
	
	nptcl = EMUtil.get_image_count(options.input)
	if nptcl<3: 
		print "ERROR: at least 3 particles are required in the input stack for all vs all. Otherwise, to align 2 particles use e2classaverage3d.py"
		sys.exit(1)
		
	allptcls = range(nptcl)
	newptcls = range(nptcl)
	for k in range(options.iter):
		'''
		Make ALL vs ALL comparisons among all NEW particles.
		NOTE: In the first round all the particles are "new"
		'''
		
		nnew = len(newptcls)
		newptclsmap = list(enumerate([range(i,nnew) for i in range(1,nnew)]))
		
		#file_map = [ (newptcls[x[0]] , [newptcls[z] for z in x[1]] ) for x in enum ]
		#numbers = [i +1 for i in range(total)]
		#num_map = [ (numbers[x(0)]] , [numbers[z] for z in x[1]] ) for x in enum ]
		
		ranks = []
		
		for ptcl1, compare in newptclsmap:
			for ptcl2 in compare:
				
				ref = EMData(options.input,ptcl1)
				task = Align3DTask(ref,["cache",options.input,ptcl2],ptcl2,"Aligning particle#%d VS particle#%d in iteration %d" % (ptcl1,ptcl2,k),options.mask,options.normproc,options.preprocess,
				options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
				tasks.append(task)
				
				#info = [name1, name2, rt2apply, score]
				#ranks.append(info)							#Add comparison to a list that will shortly be sorted by ccc score, to determine the "UNIQUE best pairs"
		
		'''
		Make comparisons for all NEW VS all OLD particles. "NEW" means particles that didn't exist in the previous round
		There are no "new" and "old" particles in the first round; thus this is needed only for k>0
		'''
		if k > 0:		
			if len(allptcls) == 1:
				print "The all vs all alignment has finalized and converged into one average"
				print "TERMINATING"
				sys.exit()
				
			setNEW = set(newptcls)
			setALL = set(allptcls)

			if options.verbose:
				print "The set of NEW particles is", setNEW
				print "The set of ALL particles is", setALL
				setOLD = setALL - setNEW
				print "Therefore, the set of all OLD particles is their difference", setOLD

			for ptcl1 in setNEW:
				for ptcl2 in setOLD:
					ref = EMData(options.input,ptcl1)
					task = Align3DTask(ref,["cache",options.input,ptcl2],ptcl2,"Aligning particle#%d VS particle#%d in iteration %d" % (ptcl1,ptcl2,k),options.mask,options.normproc,options.preprocess,
					options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
					tasks.append(task)
					
					particle1 = EMData(np,0)		#Load the particles
					particle2 = EMData(op,0)

					particle1_ed = editor(parameters, particle1,'no')[0]		#Preprocess the particles to prepare them for alignment
					particle2_ed = editor(parameters, particle2,'no')[0]

					rtresults = tomohunt(parameters, particle1, particle2)		#Actual alignment; get back only the returned transform
					rt2apply = rtresults['xform.align3d']
					score = rtresults['score']					
					info = [np, op, rt2apply, score]
					ranks.append(info)						#Add comparison to a list that will shortly be sorted by ccc score	
		
		
		# start the alignments running
		tids=etc.send_tasks(tasks)
		if options.verbose: 
			print "%d tasks queued in iteration %d"%(len(tids),k) 

		# Wait for alignments to finish and get results
		results=get_results(etc,tids,options.verbose)

		'''
		if options.verbose>2 : 
			print "Results:"
			pprint(results)

		ref=make_average(options.input,results,options.averager,options.saveali,options.keep,options.keepsig,options.verbose)		# the reference for the next iteration

		#postprocess(ref,options.mask,options.normproc,options.postprocess) #jesus
		postprocess(ref,None,options.normproc,options.postprocess) #jesus

		if options.sym!=None : 
			if options.verbose : print "Apply ",options.sym," symmetry"
			symmetrize(ref,options.sym)

		if options.savesteps :
			ref.write_image("bdb:class_%02d"%ic,it)

		if options.verbose: 
			print "Preparing final average"
		# new average
		ref=make_average(options.input,results,options.averager,options.saveali,options.keep,options.keepsig,options.verbose)		# the reference for the next iteration

		#if options.postprocess!=None : 
			#ref.process_inplace(options.postprocess[0],options.postprocess[1])     #jesus - The post process should be applied to the refinment averages. The last one is identical
												#to the output final average, so no need to apply it to ref. Plus, you ALWAYS want to have a copy
												#of the average of the raw particles, completley raw
		ref['origin_x']=0
		ref['origin_y']=0		#jesus - The origin needs to be reset to ZERO to avoid display issues in Chimera
		ref['origin_z']=0
		ref.write_image(options.output,ic)
		'''
		
		
		
		
		
		
		'''
		ranks_sorted = sorted(ranks, key=operator.itemgetter(-1), reverse=True)		#Sort pairs by ccc score
		
		ranks_filename = 'ranks_' + common + new.rstrip('_') + '.lst'			#Store the ordered rank of ALL PAIRS to a file
		ranks_file = open(ranks_filename,'w')
		
		#
		#Determine WHICH PAIRS TO AVERAGE
		#
		
		best_pairs_filename = 'best_pairs_' + common + new.rstrip('_')  + '.lst'	#Store AVERAGED PAIRS in the round for fast reference/verification
		best_pairs_file = open(best_pairs_filename,'w')
		
		done_dir = 'done_' + common + new.rstrip('_')
		os.system('mkdir ' + done_dir)
					
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
			os.system('mv ' + qq + ' ' + done_dir)
		os.system('mv *.lst ' + done_dir)

		for pp in ranks_sorted:
			pps = set(pp)
			dif = pps - averaged
			if len(dif) == 4:
				ranks.append(pp)
		k += 1
		'''
		
		
		
				
		
		
		
		E2end(logger)

		return()

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
				r=etc.get_results(tidsleft[i])			# results for a completed task
				print "\n@@@@@@The results for the completed task are", r
				ptcl=r[0].options["ptcl"]			# get the particle number from the task rather than trying to work back to it
				results[ptcl]=r[1]["final"]			# this will be a list of (qual,Transform)
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

	def __init__(self,fixedimage,image,ptcl,label,mask,normproc,preprocess,npeakstorefine,align,aligncmp,ralign,raligncmp,shrink,shrinkrefine,verbose):
		"""fixedimage and image may be actual EMData objects, or ["cache",path,number]
	label is a descriptive string, not actually used in processing
	ptcl is not used in executing the task, but is for reference
	other parameters match command-line options from e2classaverage3d.py
	Rather than being a string specifying an aligner, 'align' may be passed in as a Transform object, representing a starting orientation for refinement"""
		data={}
		data={"fixedimage":fixedimage,"image":image}
		EMTask.__init__(self,"ClassAv3d",data,{},"")

		self.options={"ptcl":ptcl,"label":label,"mask":mask,"normproc":normproc,"preprocess":preprocess,"npeakstorefine":npeakstorefine,"align":align,"aligncmp":aligncmp,"ralign":ralign,"raligncmp":raligncmp,"shrink":shrink,"shrinkrefine":shrinkrefine,"verbose":verbose}
	
	def execute(self,callback=None):
		"""This aligns one volume to a reference and returns the alignment parameters"""
		options=self.options
		if options["verbose"]: 
			print "Aligning ",options["label"]

		if isinstance(self.data["fixedimage"],EMData):
			fixedimage=self.data["fixedimage"]
		else: 
			fixedimage=EMData(self.data["fixedimage"][1],self.data["fixedimage"][2])
		
		if isinstance(self.data["image"],EMData):
			image=self.data["image"]
		else : image=EMData(self.data["image"][1],self.data["image"][2])
		
		# Preprocessing applied to both volumes.
		# Make the mask first, use it to normalize (optionally), then apply it 
		mask=EMData(image["nx"],image["ny"],image["nz"])
		mask.to_one()
		if options["mask"] != None:
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
			

		if options["verbose"]: print "Align size %d,  Refine Align size %d"%(sfixedimage["nx"],s2fixedimage["nx"])

		# If a Transform was passed in, we skip coarse alignment
		if isinstance(options["align"],Transform):
			bestcoarse=[{"score":1.0,"xform.align3d":options["align"]}]
			if options["shrinkrefine"]>1 : bestcoarse[0]["xform.align3d"].set_trans(bestcoarse[0]["xform.align3d"].get_trans()/float(options["shrinkrefine"]))
		# this is the default behavior, seed orientations come from coarse alignment
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
		
		bestfinal.sort()
		print "\n$$$$\n$$$$\n$$$$\n$$$$\n$$$$\n$$$$The best peaks sorted are", bestfinal
		
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

				#os.system('mkdir ' + parameters['ID'])	
				os.system('mv *' + parameters['ID'] + '* ' + parameters['ID'])
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
		os.system('mkdir ' + done_dir)
					
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
			os.system('mv ' + qq + ' ' + done_dir)
		os.system('mv *.lst ' + done_dir)

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



