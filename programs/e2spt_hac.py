#!/usr/bin/env python

#
# Author: Jesus Galaz, 07/2011; modified 01/March/2013
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

from EMAN2jsondb import JSTask,jsonclasses

from operator import itemgetter	

from e2spt_classaverage import sptmakepath


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <output> [options]
	
	This program depends on e2spt_classaverage.py because it imports the preprocessing and alignment functions from it.
	
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
	
	parser.add_argument("--path",type=str,default=None,help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsimjob'; for example, sptsimjob_02 will be the directory by default if 'sptsimjob_01' already exists.")
	
	parser.add_argument("--groups",type=int,default=1,help="Breaks the set into subgroups and does ALL vs ALL on the subgroups separately. Recommended when the set is > 100")
	
	parser.add_argument("--autocenter",action="store_true", help="Autocenters each averaged pair on all rounds using their center of mass",default=False)

	parser.add_argument("--input", type=str, help="The name of the input volume stack. MUST be HDF or BDB, since volume stack support is required.", default=None)
		
	parser.add_argument("--iter", type=int, help="The number of iterations to perform. Default is 1.", default=1)
	
	parser.add_argument("--exclusive_class_min", type=int, help="""The minimum multiplicity (number of particles that went into an average) to look for mutually exclusive classes/averages.
									Two classes are mutually exclusive when non of the members in one are present in the other.
									In HAC (hierarchical ascendant classification or "all vs all" alignments, classes MERGE, so a class
									from a later round will be composed of classes from earlier rounds. Some classes remain un-merged for many rounds.
									If set, this parameter will extract classes with a minimum number of particles (from whatever round/iteration they were 
									generated in) whose members are not present in any other of the extracted classes. The mutually exclusive classes
									will be put into a separate sub-directory starting with the character 'me_classes'.""", default=None)
	
	parser.add_argument("--savesteps",action="store_true", help="""If set, this will save the averages after each iteration to round#_averages.hdf. 
																There will be one .hdf stack per round, and the averages of 2 or more particles generated in that round will be images in that stack.""", default=False)
	parser.add_argument("--saveali",action="store_true", help="""If set, this will save the aligned/averaged volumes from the immediately PREVIOUS round 
																that went into the NEW particles in the "current" round, to round#_particles.hdf. 
																It will also save the latest state of alignment (for the LAST iteration only) of ALL particles provided in the input stack.
																Overwrites existing files.""",default=False)
	parser.add_argument("--saveallalign",action="store_true", help="""NOT WORKING YET: If set, will save the alignment parameters for ALL particles, 
																	aligned and unaligned (averaged and unaveraged), at each iteration""",default=False)

	#Does save ali save the stack off ALL currently UNAVERAGED particles???
	
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize")
	parser.add_argument("--preprocess",type=str,help="A processor (as in e2proc3d.py; could be masking, filtering, etc.) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.",default=None)
	
	parser.add_argument("--lowpass",type=str,help="A processor (as in e2proc3d.py; could be masking, filtering, etc.) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.",default=None)
	parser.add_argument("--highpass",type=str,help="A processor (as in e2proc3d.py; could be masking, filtering, etc.) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.",default=None)

	parser.add_argument("--npeakstorefine", type=int, help="The number of best coarse alignments to refine in search of the best final alignment. Default=4.", default=4)
	parser.add_argument("--align",type=str,help="This is the aligner use for alignments. Default is rotate_translate_3d:search=10:delta=10:dphi=10", default="rotate_translate_3d:search=10:delta=10:dphi=10")
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo")
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. Default is refine.3d, specify 'None' to disable", default="refine_3d")
	parser.add_argument("--raligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo")
	
	parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean")
	
	parser.add_argument("--postprocess",type=str,help="A processor to be applied to the volume after averaging the raw volumes, before subsequent iterations begin.",default=None)
		
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkrefine", type=int,default=1,help="Optionally shrink the input volumes by an integer amount for refine alignment.")
#	parser.add_argument("--automask",action="store_true",help="Applies a 3-D automask before centering. Can help with negative stain data, and other cases where centering is poor.")
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	parser.add_argument("--procfinelikecoarse",type=bool,default=True,help="Turn on with --procfinelikecoarse=False, and supply fine alignment parameters, such as --lowpassfine, --highpassfine, etc; to preprocess the particles for FINE alignment differently than for COARSE alignment.")
	parser.add_argument("--randomizewedge",action="store_true", help="This parameter is EXPERIMENTAL. It randomizes the position of the particles BEFORE alignment, to minimize missing wedge bias and artifacts during symmetric alignment where only a fraction of space is scanned", default=False,)


	parser.add_argument("--minscore",type=float,help="""Percent of the maximum score to use as a threshold for the minimum score to allow.
													For example, if the best pair in the first iteration yielded a score of 15.0, and you supply --minscore=0.666,
													any pair wise alignments with a score lower than 15*0.666=10 will be forbidden.""", default=0.15)

	'''
	Parameters to compensate for the missing wedge using --cpm=fsc.tomo
	'''
	parser.add_argument("--wedgeangle",type=float,help="""Missing wedge angle, calculated as 90 minus the value yo provide, times 2. 
														For example, --wedgeangle=60 will represent a wedge of size (90-60)*2=60.
														--wedgeangle=70, results in a narrower wedge of size (90-70)*2=40.
														In reality, you should enter here the range of your DATA COLLECTION.
														I.e., if you collected your tiltseries from -60 to 60, enter --wedgeangle=60.""",default=60.0)
	parser.add_argument("--wedgei",type=float,help="Missingwedge begining (in terms of its 'height' along Z. If you specify 0, the wedge will start right at the origin.", default=0.15)
	parser.add_argument("--wedgef",type=float,help="Missingwedge ending (in terms of its 'height' along Z. If you specify 1, the wedge will go all the way to the edge of the box.", default=0.9)
	parser.add_argument("--fitwedgepost", action="store_true", help="Fit the missing wedge AFTER preprocessing the subvolumes, not before, IF using the fsc.tomo comparator for --aligncmp or --raligncmp.", default=False)
	




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
	
	if options.lowpass: 
		options.lowpass=parsemodopt(options.lowpass)

	if options.highpass: 
		options.highpass=parsemodopt(options.highpass)

	if options.postprocess: 
		options.postprocess=parsemodopt(options.postprocess)

	'''
	Make the directory where to create the database where the results will be stored
	'''
		
	options = sptmakepath(options,'spt_hac')

	group_ranges=[]
	data_files = []
	
	nptcl = EMUtil.get_image_count(options.input)
	entirestack = options.input
	originalpath = options.path
	
	logger = E2init(sys.argv,options.ppid)

	for i in range(options.groups):	
		#if options.groups > 1:
		if options.groups * 3 > nptcl:
			print "ERROR: You need at least 3 particles per groups to do all vs all within each group."
			print "You asked for %d groups; thus, the stack needs to have more than %d particles, but it only has %d" % (options.groups,3*options.groups,nptcl)
			print "Reduce the number of groups requested or provide a larger stack."
			sys.exit()
		else:
			groupsize = int( nptcl/options.groups )
			bottom_range = i * groupsize
			top_range = (i+1) * groupsize
			if i == options.groups - 1:
				top_range = nptcl
			
			groupPATH = entirestack
			if options.groups > 1:
				groupDIR = originalpath + '/group' + str(i+1).zfill(len(str(options.groups)))
				groupID = 'group' + str(i+1).zfill(len(str(options.groups))) + '_raw.hdf'
				groupPATH = groupDIR + '/' + groupID
				os.system('mkdir ' + groupDIR)				

				options.path = groupDIR

				print "I will start ALL vs ALL on group number", i+1
				print "For which options.input is", options.input

			mm = 0
			for jj in xrange(bottom_range,top_range):
				print "I am rewritting the spt_ptcl_indxs header parameter for every particle in the stack"
				a = EMData(entirestack,jj)
				a['spt_ptcl_indxs'] = mm
	
				params = a.get_attr_dict()
				
				for i in params:
					if 'sptID' in i:
						print "I'm resetting THIS parameter", i
						a[str(i)] = ''
						
					if 'spt_ID' in i:
						print "I'm resetting THIS parameter", i
						a[str(i)] = ''						
						
						
				print "The type of a should be EMData", type(a)
				a.write_image(groupPATH,mm)

				mm += 1

			options.input = groupPATH
		
		#print "\nTHe len of options is ", len(options)
		print "\n\nAnd options are", options
		allvsall(options)
		if options.exclusive_class_min:
			exclusive_classes(options)
	E2end(logger)

	return()


def exclusive_classes(options):

	findir = os.listdir(options.path)
	
	minmembers = options.exclusive_class_min

	averages = []
	for i in findir:
		if '_averages.hdf' in i:
			n=EMUtil.get_image_count(options.path + '/' + i)
			for j in range(n):
				a=EMData(options.path + '/' +i,j,True)
				print "The average I have just appended has this multiplicity", a['spt_multiplicity']
				print "And comes from this file", options.path + '/' + i
				indxs = a['spt_ptcl_indxs']
				indxs.sort()
				averages.append({'file':options.path + '/' + i,'n':j,'multiplicity':a['spt_multiplicity'],'indxs':indxs})
	me_classes = []
	repeated = []
	
	averages = sorted(averages, key=itemgetter('multiplicity'))
	
	candidates = 0
	repeats=0
	for i in range(len(averages)):
		if int(averages[i]['multiplicity']) >= int(options.exclusive_class_min):	
			print "I've found a candidate with the min number of members in the average!", averages[i]['multiplicity']
			candidates+=1
			for j in xrange(i+1,len(averages)):	
				if averages[i] not in repeated:
					#print "I'm comparing these candidates", averages[i], averages[j]			
					for k in averages[i]['indxs']:
						if k in averages[j]['indxs']:
							if averages[j] not in repeated:							
								#print "I have added this candidate to repeated", averages[j]					
								print "I've found a repeated..."
								repeats+=1							
								repeated.append(averages[j])
								break
				
			if averages[i] not in me_classes and averages[i] not in repeated:					
				print "\n@@@\nI have found a ME class $$$$$$"			
				me_classes.append(averages[i])

	print "The mutually exclusive classes with more than %d particles" %options.exclusive_class_min
	print "are:"
	for i in me_classes:
		print i

	udir = options.path + '/me_classes_' + str(options.exclusive_class_min).zfill(len(str(options.exclusive_class_min)))
	os.system('mkdir ' + udir)

	for i in range(len(me_classes)):
		out = udir + '/me_class' + str(i).zfill(len(str(len(me_classes)))) + '_s' + str(me_classes[i]['multiplicity']) + '.hdf'
		cmd = 'e2proc3d.py ' + me_classes[i]['file'] + ' ' + out + ' --first=' + str(me_classes[i]['n']) + ' --last=' + str(me_classes[i]['n']) + ' --append'
		os.system(cmd)
		if options.postprocess:
			a=EMData(out,0)
			a.process_inplace(options.postprocess[0],options.postprocess[1])
			outpp = out.replace('.hdf','_postp.hdf')
			a.write_image(outpp,0)

	return()


def allvsall(options):
	print "These are path and input received in allvsall", options.path, options.input
	
	print "With these many particles in it", EMUtil.get_image_count(options.input)
	
	print "\nI will load the header of a particle."
	
	hdr = EMData(options.input,0,True)
	nx = int(hdr["nx"])
	ny = int(hdr["ny"])
	nz = int(hdr["nz"])
	if nx!=ny or ny!=nz :
		print "ERROR, input volumes are not cubes"
		sys.exit(1)
	
	print "Officialy counting number of particles"
	nptcl = EMUtil.get_image_count(options.input)
	if nptcl<3: 
		print "ERROR: at least 3 particles are required in the input stack for all vs all. Otherwise, to align 2 particles (one to the other or to a model) use e2spt_classaverage.py"
		sys.exit(1)
	
	fillfactor = len(str(nptcl))							#Calculate this based on the number of particles so that tags are adequate ("pretty") and ordered
	roundtag='round' + str(0).zfill(fillfactor)					#We need to keep track of what round we're in
	newptcls={}													#This dictionary stores 'new particles' produced in each round as { particleID : particleDATA } elements
	allptclsRound={}								#This dictionary stores all particlces in a round ("new" and "old") as 
											#{particle_id : [EMData,{index1:totalTransform1, index2:totalTransform2...}]} elements
											#The totalTransform needs to be calculated for each particle after each round, to avoid multiple interpolations
	
	print "Starting the loop"
	for i in range(nptcl):								#In the first round, all the particles in the input stack are "new" and should have an identity transform associated to them
		a=EMData(options.input,i)
		totalt=Transform()

		if 'spt_multiplicity' not in a.get_attr_dict():				#spt_multiplicity keeps track of how many particles were averaged to make any given new particle (set to 1 for the raw data)
			a['spt_multiplicity']=1
		elif not a['spt_multiplicity']:
			a['spt_multiplicity']=1
		#if 'spt_ptcl_indxs' not in a.get_attr_dict():				#spt_ptcl_indxs keeps track of what particles from the original stack went into a particular average or "new particle"
		a['spt_ptcl_indxs']=[i]							#The raw stack should contain particles where this parameter is the particle number itself
		#else:
		if type(a['spt_ptcl_indxs']) is int:
			a['spt_ptcl_indxs'] = [a['spt_ptcl_indxs']]			#The spt_ptcl_indxs parameter should be of type 'list', to easily 'append' new particle indexes
		
		particletag = roundtag + '_' + str(i).zfill(fillfactor)
		newptcls.update({particletag :a})
		
		if 'spt_ID' not in a.get_attr_dict():					#spt_multiplicity keeps track of how many particles were averaged to make any given new particle (set to 1 for the raw data)
			a['spt_ID'] = particletag
		elif not a['spt_ID']:
			a['spt_ID'] = particletag	
		
		a.write_image(options.input,i)						#Overwrite the raw stack with one that has the appropriate header parameters set to work with e2spt_hac	
		
		allptclsRound.update({particletag : [a,{i:totalt}]})			
		
	oldptcls = {}									#'Unused' particles (those that weren't part of any unique-best-pair) join the 'oldptcls' dictionary onto the next round
	surviving_results = []								#This list stores the results for previous alignment pairs that weren't used, so you don't have to recompute them
		
	allptclsMatrix = []								#Massive matrix listing all the allptclsRound dictionaries produced after each iteration
	allptclsMatrix.append(allptclsRound)
	
	
	FinalAliStack = {}								#This will keep track of the most updated alignment state of ALL particles regardless of whether they've participated 
													#in averages or not.
													
	nptcls = EMUtil.get_image_count(options.input)
	for i in range(nptcls):
		rawptcl = EMData(options.input,i)
		rawptcl['xform.alignd3d']=Transform()
		FinalAliStack.update({i:rawptcl})
	
	maxScore = 1
	for k in range(options.iter):							#Start the loop over the user-defined number of iterations
		#avgname = options.path + '/round' + str(k).zfill(fillfactor) + '_averages.hdf'
		newstack = options.path + '/round' + str(k-1).zfill(fillfactor) + '_averages.hdf'
		if k== 0:
			newstack =options.input
		
		print "\n\n\n\n\n\n\n$$$$$$$$$$$$$$$$$\n$$$$$$$$$$$$$$$$$\nStarting this iteration!", k
		print "\n\n"
		
		nnew = len(newptcls)
		if k == (int(options.iter) - 1) or (nnew + len(oldptcls) ) == 1 :
			print "This is the final round", k
			if options.saveali:
				print "You selected ; therefore, I will write the latest state of all particles in the inpust stack."
			
				for key in FinalAliStack:
					aliptcl = FinalAliStack[key]
					aliptcl.write_image(options.path + '/finalAliStack.hdf',key)
					print "Wrote this ptcl to final stack", key
		
		if nnew + len(oldptcls) == 1:						#Stop the loop if the data has converged and you're left with one final particle (an average of all)
			print "The all vs all algorithm has converged into one average"
			break
		
		allptclsRound = {}							
		
		
		print "\nInitialize parallelism"
		if options.parallel:							# Initialize parallelism if being used
			from EMAN2PAR import EMTaskCustomer
			etc=EMTaskCustomer(options.parallel)
			pclist=[options.input]
			etc.precache(pclist)
		tasks = []
		
		'''
		Make ALL vs ALL comparisons among all NEW particle INDEXES.
		NOTE: In the first round all the particles are "new"
		'''
		
		newptclsmap = list(enumerate([range(i,nnew) for i in range(1,nnew)]))
		
		jj=0									#Counter to track the number of comparisons (also the number of tasks to parallelize)
		roundtag = 'round' + str(k).zfill(fillfactor) + '_'			#The round tag needs to change as the iterations/rounds progress
		
		print "\n Start all vs all comparisons"
		for ptcl1, compare in newptclsmap:
			for ptcl2 in compare:
				
				reftag = roundtag + str(ptcl1).zfill(fillfactor)								
				particletag = roundtag + str(ptcl2).zfill(fillfactor)
				
				#if options.verbose > 2:
				print "Setting the following comparison: %s vs %s in ALL VS ALL" %(reftag,particletag)
				
				#def __init__(self,fixedimagestack,imagestack,comparison, ptcl1, ptcl2, p1n, p2n,label,options,transform):
				
				task = Align3DTaskAVSA(newstack,newstack, jj, reftag, particletag, ptcl1, ptcl2,"Aligning particle#%s VS particle#%s in iteration %d" % (reftag,particletag,k),options,k)

				#task = Align3DTaskAVSA(newstack,newstack, jj, reftag, particletag, ptcl1, ptcl2,"Aligning particle#%s VS particle#%s in iteration %d" % (reftag,particletag,k),options.mask,options.normproc,options.preprocess,options.lowpass,options.highpass,
				#options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
				
				tasks.append(task)
				
				jj+=1
		
		'''
		Make comparisons for all NEW VS all OLD particles. "NEW" means a particle that didn't exist in the previous round.
		There are no "new" and "old" particles in the first round; thus the loop below is needed only for k>0
		'''
				
		if k > 0:
		
			oldtags = {}
			nold = EMUtil.get_image_count(options.path + '/oldptclstack.hdf')
			
			for i in range(nold):
				oldtags.update({EMData(options.path + '/oldptclstack.hdf',i,True)['spt_ID'] : i})
			
			
			print "Old tagas are:\n", oldtags
			nnn = 0
			for refkey,refvalue in newptcls.iteritems():
				ptcl1 = nnn
				for particlekey,particlevalue in oldptcls.iteritems():
					
					ptcl2 = oldtags[particlekey]
					
					#if options.verbose > 2:
					print "Setting the following comparison: %s vs %s in ALL VS ALL" %(refkey,particlekey)
					
					#task = Align3DTaskAVSA( newstack, options.path + '/oldptclstack.hdf', jj , refkey, particlekey, ptcl1, ptcl2,"Aligning particle round#%d_%d VS particle#%s, in iteration %d" % (k,ptcl1,particlekey.split('_')[0] + str(ptcl2),k),options.mask,options.normproc,options.preprocess,options.lowpass,options.highpass,
					#options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
					
					task = Align3DTaskAVSA( newstack, options.path + '/oldptclstack.hdf', jj , refkey, particlekey, ptcl1, ptcl2,"Aligning particle round#%d_%d VS particle#%s, in iteration %d" % (k,ptcl1,particlekey.split('_')[0] + '_' + str(ptcl2),k),options,k)
					
					
					tasks.append(task)
										
					jj+=1
				nnn+=1	
		
		tids=etc.send_tasks(tasks)						#Start the alignments running
		#if options.verbose > 0: 
		print "%d tasks queued in iteration %d"%(len(tids),k) 
		
		results = get_results(etc,tids,options.verbose)				#Wait for alignments to finish and get results
		#results = ret[0]
		results = results + surviving_results					#The total results to process/analyze includes results (comparisons) from previous rounds that were not used
		results = sorted(results, key=itemgetter('score'))			#Sort/rank the results by score
		
		#print "results are", results
	
		if k == 0:
			plotX=[]
			plotY=[]
	
			simmxScores = EMData(nptcls,nptcls)
			simmxXs = EMData(nptcls,nptcls)
			simmxYs = EMData(nptcls,nptcls)
			simmxZs = EMData(nptcls,nptcls)
			simmxAzs = EMData(nptcls,nptcls)
			simmxAlts = EMData(nptcls,nptcls)
			simmxPhis = EMData(nptcls,nptcls)
			simmxScales = EMData(nptcls,nptcls)
			
			simmxScores.to_zero()
			simmxXs.to_zero()
			simmxYs.to_zero()
			simmxZs.to_zero()
			simmxAzs.to_zero()
			simmxAlts.to_zero()
			simmxPhis.to_zero()
			simmxScales.to_one()
		
		compNum=0
		for i in results:
			if options.verbose > 0:
				print "In iteration %d the SORTED results are:", k	
				print "%s VS %s , score=%f, transform=%s" %(i['ptclA'], i['ptclB'], i['score'], i['xform.align3d'] )
			
			if k == 0:
				if compNum == 0:
					maxScore = i['score']
				plotX.append( compNum )
				plotY.append( i['score'] )
				
				indxA = int( i['ptclA'].split('_')[-1] )
				indxB = int( i['ptclB'].split('_')[-1] )
				simmxScores.set_value_at(indxA,indxB,i['score'])
				
				t= i['xform.align3d']
				trans=t.get_trans()
				rots=t.get_rotation()
				
				simmxXs.set_value_at(indxA,indxB,float(trans[0]))
				simmxYs.set_value_at(indxA,indxB,float(trans[1]))
				simmxZs.set_value_at(indxA,indxB,float(trans[2]))
				simmxAzs.set_value_at(indxA,indxB,float(rots['az']))
				simmxAlts.set_value_at(indxA,indxB,float(rots['alt']))
				simmxPhis.set_value_at(indxA,indxB,float(rots['phi']))
				simmxScales.set_value_at(indxA,indxB,float(i['score']))
			
				compNum+=1
		
		if k == 0:
			simmxFile = 'simmx_' + str( k ).zfill( len (str (options.iter))) + '.hdf'
			simmxScores.write_image(simmxFile,0)
			simmxXs.write_image(simmxFile,1)
			simmxYs.write_image(simmxFile,2)
			simmxZs.write_image(simmxFile,3)
			simmxAzs.write_image(simmxFile,4)
			simmxAlts.write_image(simmxFile,5)
			simmxPhis.write_image(simmxFile,6)	
			simmxScales.write_image(simmxFile,7)
			
			from e2figureplot import plotter
			plotter (plotX, plotY, options, 'scatter')
			
			plotName = simmxFile.replace('.hdf','_PLOT.png')
			from e2figureplot import textwriter
			textwriter(plotX, plotY, options, plotName)
	
		print "\n\n\n\nIn iteration %d, the total number of comparisons in the ranking list, either new or old that survived, is %d" % (k, len(results))
		
		tried = set()											#Tracks what particles have "appeared" on the list, whether averaged or not
		averages = {}											#Stores the new averages; you need a dict different from newptcls because you'll need to 'fetch' data from the previous version of newptcls
		used = set()											#Tracks what particles WERE actually averaged
		
		mm=0												#Counter to track new particles/averages produced and write them to output
		print "I'm in the averager!!!!!!!!!!!!"
		
		for z in range(len(results)):
			if options.minscore:
				score = results[z]['score']			
				if score > maxScore * float( options.minscore ):
					print "Breaking loop because the next comparison score", score
					print "Is worse (larger, more positive, in EMAN2) than the specified percentage", options.minscore
					print "Of the maximum score from the initial simmx matrix", maxScore
					break
			
			if results[z]['ptclA'] not in tried and results[z]['ptclB'] not in tried:
				tried.add(results[z]['ptclA'])							#If the two particles in the pair have not been tried, and they're the next "best pair", they MUST be averaged
				tried.add(results[z]['ptclB'])							#Add both to "tried" AND "used" 
				used.add(results[z]['ptclA'])		
				used.add(results[z]['ptclB'])
													
				avgr = Averagers.get(options.averager[0], options.averager[1])			#Call the averager
				
				avg_ptcls = []
								
				ptcl1 = allptclsMatrix[k][results[z]['ptclA']][0]				
														#You always add all the past particles that went into a particular new particle (for all the particles being
														#averaged in the current round) freshly from the raw stack to the averager (with
														#the appropriate transforms they've undergone, of course. Thus, YOU DON'T have to
														#worry about "multiplicity", since it takes care of itself by doing this.
				indx_trans_pairs = {}

				#print "\n\n\nThe indexes in particle 1 are", ptcl1['spt_ptcl_indxs']
				
				ptcl1info = allptclsMatrix[k][results[z]['ptclA']]
				#print "\n\nptcl1 info attached is", ptcl1info
				
				ptcl1_indxs_transforms = ptcl1info[-1]
				#print "\n\nptcl1_indexes_transforms is", ptcl1_indxs_transforms
				#print "\n\n"
								
				for p in ptcl1['spt_ptcl_indxs']:											
					pastt = ptcl1_indxs_transforms[p]
					
					subp1 = EMData(options.input,p)
					
					#subp1.process_inplace("xform",{"transform":pastt})					#RADICAL CHANGE ****************************************					
					subp1.transform(pastt)
					
					subp1['xform.align3d']=pastt
					
					avgr.add_image(subp1)
					
					if options.saveali:
						avg_ptcls.append(subp1)
				
						subp1_forFinalAliStack =subp1.copy()							#It would be insane to write out ALL the particles in each iteration in the final aligned state
						subp1_forFinalAliStack['xform.align3d']=pastt
						
						FinalAliStack.update({int(p):subp1_forFinalAliStack})		#But let's keep track of them, and write them out only when the LAST iteration has been reached
						print "I have CHANGED a particle1 in finalAliStack to have this transform", totalt

					indx_trans_pairs.update({p:pastt})
					
						
				#avgr.add_image(ptcl1)								#Add particle 1 to the average
				
				ptcl2 = allptclsMatrix[k][results[z]['ptclB']][0]
						
				resultingt = results[z]["xform.align3d"]					
					
				print "\n\n\nThe indexes in particle 2 are", ptcl2['spt_ptcl_indxs']
				
				ptcl2info = allptclsMatrix[k][results[z]['ptclB']]
				
				print "\n\nptcl2 info attached is", ptcl2info
					
				ptcl2_indxs_transforms = ptcl2info[-1]
				
				print "\n\nptcl2_indexes_transforms is", ptcl2_indxs_transforms
				print "\n\n"
				
				#ptcl2avgr = Averagers.get(options.averager[0], options.averager[1])		#You need to recompute ptcl2 "fresh" from the raw data to avoid multiple interpolations
				
				for p in ptcl2['spt_ptcl_indxs']:						#All the particles in ptcl2's history need to undergo the new transformation before averaging
					print "I'm fixing the transform for this index in the new average", p	#(multiplied by any old transforms, all in one step, to avoid multiple interpolations)
					
					pastt = ptcl2_indxs_transforms[p]
					
					print "\n\n\n@@@@@@@@@@@@@@@@@"
					print "So the past transform for ptcl2 was", pastt
					print "But the resulting transform is", resultingt
					
					totalt = resultingt * pastt
					print "Which means their product is", totalt

					subp2 = EMData(options.input,p)
					#subp2.process_inplace("xform",{"transform":totalt})					#RADICAL CHANGE **********************************
					subp2.transform(totalt)
					subp2['xform.align3d']=totalt
					
					print "Which should coincide with the xform.align3d parameter on the rotated particle's header", subp2['xform.align3d']
					
					avgr.add_image(subp2)
					
					
					if options.saveali:
						avg_ptcls.append(subp2)

						subp2_forFinalAliStack =subp2.copy()							#It would be insane to write out ALL the particles in each iteration in the final aligned state
						subp2_forFinalAliStack['xform.align3d']=totalt
						
						FinalAliStack.update({int(p):subp2_forFinalAliStack})			#But let's keep track of them, and write them out only when the LAST iteration has been reached			
						print "I have CHANGED a particle2 in finalAliStack to have this transform", totalt

					indx_trans_pairs.update({p:totalt})
								
				avg=avgr.finish()
				
				
				
				
				print "THe average was successfully finished"

				if options.autocenter:
					print "\n\n\n\nYou have selected to autocenter!\n"
					avg = avg.process('xform.centerofmass')
					tcenter = avg['xform.align3d']
					print "Thus the average HAS BEEN be translated like this", tcenter
					
					if options.saveali:
						avg_ptcls = []
					
					for p in ptcl1['spt_ptcl_indxs']:
						pastt = ptcl1_indxs_transforms[p]
						totalt = tcenter * pastt
						indx_trans_pairs.update({p:totalt})
						
						
						subp1 = EMData(options.input,p)
						subp1.process_inplace("xform",{"transform":totalt})
						avg_ptcls.append(subp1)
						
						if options.saveali:	
							subp1_forFinalAliStack =subp1.copy()							#It would be insane to write out ALL the particles in each iteration in the final aligned state
							subp1_forFinalAliStack['xform.align3d']=totalt
							
							FinalAliStack.update({int(p):subp1_forFinalAliStack})		#But let's keep track of them, and write them out only when the LAST iteration has been reached
							print "After AUTOCENTER have CHANGED a particle1 in finalAliStack to have this transform", totalt

					
						
					for p in ptcl2['spt_ptcl_indxs']:
						pastt = ptcl2_indxs_transforms[p]
						totalt = tcenter * resultingt * pastt
						indx_trans_pairs.update({p:totalt})
						
						subp2 = EMData(options.input,p)
						subp2.process_inplace("xform",{"transform":totalt})
						avg_ptcls.append(subp2)
						
						if options.saveali:
							subp2_forFinalAliStack =subp2.copy()							#It would be insane to write out ALL the particles in each iteration in the final aligned state
							subp2_forFinalAliStack['xform.align3d']=totalt
							
							FinalAliStack.update({int(p):subp2_forFinalAliStack})		#But let's keep track of them, and write them out only when the LAST iteration has been reached
							print "After AUTOCENTER I have CHANGED a particle2 in finalAliStack to have this transform", totalt
				
				
				
				
				
				
				
				
				print "I will set the multiplicity of the average"
				avgmultiplicity = ptcl1['spt_multiplicity'] + ptcl2['spt_multiplicity']		#Define and set the multiplicity of the average
				avg['spt_multiplicity'] = avgmultiplicity
				
				print avgmultiplicity

				indexes1 = ptcl1["spt_ptcl_indxs"]
				indexes2 = ptcl2["spt_ptcl_indxs"]				
				avgindexes = indexes1 + indexes2

				print "I will sort the indexes in the average"				
				avgindexes.sort()
				print avgindexes

				avg["spt_ptcl_indxs"] = avgindexes						#Keep track of what particles go into each average or "new particle"				
				
				avg["spt_ptcl_src"] = options.input
				
				avg['origin_x'] = 0								#The origin needs to be set to ZERO to avoid display issues in Chimera
				avg['origin_y'] = 0
				avg['origin_z'] = 0
							
				newroundtag = 'round' + str(k+1).zfill(fillfactor) + '_'
				avgtag = newroundtag + str(mm).zfill(fillfactor)
				
				avg['spt_ID'] = avgtag
				
				print "I will normalize the average"
				
				# Make the mask first, use it to normalize (optionally), then apply it 
				#mask=EMData(avg['nx'],avg['ny'],avg['nz'])
				#mask.to_one()
				
				#if options.mask:
				#	#print "This is the mask I will apply: mask.process_inplace(%s,%s)" %(options.mask[0],options.mask[1]) 
				#	mask.process_inplace(options.mask[0],options.mask[1])
		
				# normalize
				#if options.normproc:
				#	if options.normproc[0]=="normalize.mask": 
				#		options.normproc[1]["mask"]=mask

				#	avg.process_inplace(options.normproc[0],options.normproc[1])

				#Mask after normalizing with the mask you just made, which is just a box full of 1s if not mask is specified
				#avg.mult(mask)

				#If normalizing, it's best to do normalize-mask-normalize-mask
				#if options.normproc:
				#	if options.normproc[0]=="normalize.mask": 
				#		options.normproc[1]["mask"]=mask

				#	avg.process_inplace(options.normproc[0],options.normproc[1])

				#	avg.mult(mask)
				
				avg.process_inplace(options.normproc[0],options.normproc[1])
				
				avg.write_image(options.path + '/round' + str(k).zfill(fillfactor) + '_averages.hdf',mm)
				
				if options.saveali:
					for oo in range(len(avg_ptcls)):
						avg_ptcls[oo].write_image(options.path + '/round' + str(k).zfill(fillfactor) + '_average' + str(mm).zfill(fillfactor)  + '_ptcls.hdf',oo)
				
				if options.postprocess: 
					avgp=avg.process(options.postprocess[0],options.postprocess[1])
					avgp.write_image(options.path + '/round' + str(k).zfill(fillfactor) + '_averages_postp.hdf',mm)
										
				averages.update({avgtag:avg})	   						#The list of averages will become the new set of "newptcls"
				allptclsRound.update({avgtag : [avg,indx_trans_pairs]})
				
				mm+=1
				
			if results[z]['ptclA'] not in tried:						#If a particle appeared in the ranking list but its pair was already taken, the particle must be classified as "tried"
				tried.add(results[z]['ptclA'])						#because you don't want to average it with any other available particle lower down the list that is available
													#We only average "UNIQUE BEST PAIRS" (the first occurance in the ranking list of BOTH particles in a pair).
			if results[z]['ptclB'] not in tried:
				tried.add(results[z]['ptclB'])
		
		
		surviving_results = []
		for z in range(len(results)):
			if results[z]['ptclA'] not in used and results[z]['ptclB'] not in used:
				surviving_results.append(results[z])			
		
		surviving_newptcls = {}
		surviving_oldptcls = {}		
		
		if options.verbose > 2:
			print "These were the particles in iteration", k
		
		for particlekey,particlevalue in newptcls.iteritems():
			
			if options.verbose > 2:
				print particlekey
			
			if particlekey not in used:
				surviving_newptcls.update({particlekey:particlevalue})

			else:
				if options.verbose > 1:
					print "This particle from newptcls was averaged", particlekey
		
		for particlekey,particlevalue in oldptcls.iteritems():
			if particlekey not in used:
				surviving_oldptcls.update({particlekey:particlevalue})
			else:
				if options.verbose > 1:
					print "This particle from oldptcls was averaged", particlekey
						
		if options.verbose > 0:
			print "At the end of iteration", k
			print "There were these many old ptcls NOT averaged", len(surviving_oldptcls)
			print "And these many 'new ptcls' not averaged that need to become old", len(surviving_newptcls)
		
		oldptcls = {}
		oldptcls.update(surviving_oldptcls)  
		oldptcls.update(surviving_newptcls)					#All the particles from the newptcls list that were not averaged become "old"
				
		newptcls = averages							#All the new averages become part of the new "newptcls" list
		
		fs=os.listdir( options.path)
		
		if 'oldptclstack.hdf' in fs:						#This seems not to be needed in the last round. Sometimes the program fails with an error saying there's no old stack in the directory.
			os.system('rm ' + options.path + '/oldptclstack.hdf')
		#a=EMData(nx,ny,nz)
		#a.write_image(oldptclstack.hdf,0)					#
		
		gg=0
		for particlekey,particlevalue in oldptcls.iteritems():
			allptclsRound.update({ particlekey: [particlevalue,allptclsMatrix[k][particlekey][-1]]})
			particlevalue.write_image(options.path + '/oldptclstack.hdf',gg)
			gg+=1
			
		allptclsMatrix.append(allptclsRound)
					
		
		print "And these many new averages", len(newptcls), len(averages)
		
		print "So there are these many old particles for the next round", len(oldptcls)
		print "And these many new-new ones", len(newptcls)
		
		#if k>0:
		#	os.system('rm ' + newstack)

	
	return()


def get_results(etc,tids,verbose):
	"""This will get results for a list of submitted tasks. Won't return until it has all requested results.
	aside from the use of options["ptcl"] this is fairly generalizable code. """
	
	# wait for them to finish and get the results
	# results for each will just be a list of (qual,Transform) pairs
	results=[0]*len(tids)		# storage for results
	ncomplete=0
	
	tidsleft=tids[:]
	
	numtides = len(tids)
	while 1:
		time.sleep(5)
		proglist=etc.check_task(tidsleft)
		nwait=0
		for i,prog in enumerate(proglist):
			if prog==-1: 
				nwait+=1
			
			if prog==100:
				r=etc.get_results(tidsleft[i])						# results for a completed task
				comparison=r[0].classoptions["comparison"]					# get the comparison number from the task rather than trying to work back to it
				
				results[comparison]=r[1]["final"][0]					# this will be a list of (qual,Transform), containing the BEST peak ONLY
				
				results[comparison]['ptclA']=r[0].classoptions['ptclA']			#Associate the result with the pair of particles involved
				results[comparison]['ptclB']=r[0].classoptions['ptclB']

				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		
		print "  %d tasks, %d complete, %d waiting to start        \r"%(len(tids),ncomplete,nwait)

		if verbose:
			sys.stdout.flush()
	
		if len(tidsleft)==0 or ncomplete == numtides: 
			break		
	return results


class Align3DTaskAVSA(JSTask):
	"""This is a task object for the parallelism system. It is responsible for aligning one 3-D volume to another, with a variety of options"""
	
	
	#def __init__(self,fixedimagestack,imagestack,comparison,ptcl1,ptcl2,p1n,p2n,label,mask,normproc,preprocess,lowpass,highpass,npeakstorefine,align,aligncmp,ralign,raligncmp,shrink,shrinkrefine,verbose):

	def __init__(self,fixedimagestack,imagestack,comparison, ptclA, ptclB, pAn, pBn,label,classoptions, round):
		
		data={}
		data={"fixedimage":fixedimagestack,"image":imagestack}
		JSTask.__init__(self,"SptHac",data,{},"")

		#self.options={"comparison":comparison,"ptcl1":ptcl1,"ptcl2":ptcl2,"p1number":p1n,"p2number":p2n,"label":label,"mask":mask,"normproc":normproc,"preprocess":preprocess,"lowpass":lowpass,"highpass":highpass,"npeakstorefine":npeakstorefine,"align":align,"aligncmp":aligncmp,"ralign":ralign,"raligncmp":raligncmp,"shrink":shrink,"shrinkrefine":shrinkrefine,"verbose":verbose}
		self.classoptions={"comparison":comparison,"ptclA":ptclA,"ptclB":ptclB,"pAn":pAn,"pBn":pBn,"label":label,"classoptions":classoptions, 'round':round}
	
	def execute(self,callback=None):
		
		"""This aligns one volume to a reference and returns the alignment parameters"""
		#classoptions=self.classoptions
		options=self.classoptions
		
		"""
		CALL the alignment function, which is imported from e2spt_classaverage
		"""
		
		print "Will import alignment"
		from e2spt_classaverage import alignment
				
		print "I have imported alignment and will call it"
		
		#def alignment(fixedimage,image,ptcl,label,classoptions,transform):
		
		fixedimage = EMData( self.data["fixedimage"], options['pAn'] )
		image = EMData( self.data["image"], options['pBn'] )
		
		nptcls = EMUtil.get_image_count( options['classoptions'].input )
		
		if options['classoptions'].groups:
			nptcls = ( nptcls / int(options['classoptions'].groups) ) + nptcls % int(options['classoptions'].groups)
		
		potentialcomps = ( nptcls * (nptcls - 1) )/ 2
		
		xformslabel = 'round' + str(options['round']).zfill( len( str(options['classoptions'].iter))) + '_comparison' + str(options['comparison']).zfill( len( str(potentialcomps) ) ) + '_ptclA' + str(options['pAn']).zfill( len(str(nptcls))) + '_ptclB' + str(options['pBn']).zfill( len(str(nptcls)))
		ret=alignment( fixedimage, image, options['label'], options['classoptions'],xformslabel,None,'e2spt_hac')
		
		bestfinal=ret[0]
		bestcoarse=ret[1]
		
		return {"final":bestfinal,"coarse":bestcoarse}
		

jsonclasses["Align3DTaskAVSA"]=Align3DTaskAVSA.from_jsondict

if __name__ == '__main__':
	main()



