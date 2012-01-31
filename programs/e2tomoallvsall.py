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
	
	parser.add_argument("--groups",type=int,default=1,help="Breaks the set into subgroups and does ALL vs ALL on the subgroups separately. Recommended when the set is > 100")
	
	parser.add_argument("--autocenter",action="store_true", help="Autocenters each averaged pair on all rounds using their center of mass",default=False)

	parser.add_argument("--input", type=str, help="The name of the input volume stack. MUST be HDF or BDB, since volume stack support is required.", default=None)
		
	parser.add_argument("--iter", type=int, help="The number of iterations to perform. Default is 1.", default=1)
	
	parser.add_argument("--savesteps",action="store_true", help="If set, will save the averages after each iteration to round#_averages.hdf. There will be one .hdf stack per round, and the averages of 2 or more particles generated in that round will be images in that stack",default=False)
	parser.add_argument("--saveali",action="store_true", help="If set, will save the aligned particle volumes in round#_particles.hdf. Overwrites existing file.",default=False)
	
	#Does save ali save the stack off ALL currently UNAVERAGED particles???
	
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", default="mask.sharp:outer_radius=-2")
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	parser.add_argument("--preprocess",type=str,help="A processor (as in e2proc3d.py; could be masking, filtering, etc.) to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.",default=None)
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
		
	#if options.path and options.path[:4].lower()!="bdb:": 
	#	options.path="bdb:"+options.path
	
	files=os.listdir(os.getcwd())

	if not options.path: 
		#options.path="bdb:"+numbered_path("sptavsa",True)
		options.path = "sptavsa_01"
	else:
		while options.path in files:
			if '_' not in options.path:
				options.path = options.path + '_00'
				
			options.path = options.path.split('_')[0] + '_' + str(int(options.path.split('_')[-1]) + 1).zfill(2)
			print "The new options.path is", options.path

	if options.path not in files:
		os.system('mkdir ' + options.path)
	
	group_ranges=[]
	data_files = []
	
	nptcl = EMUtil.get_image_count(options.input)
	groupsize = nptcl
	entirestack = options.input
	originalpath = options.path
	
	for i in range(options.groups):	
		if options.groups > 1:
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

				groupDIR = originalpath + '/group' + str(i+1).zfill(len(str(options.groups)))
				groupID = 'group' + str(i+1).zfill(len(str(options.groups))) + '_raw.hdf'
				groupPATH = groupDIR + '/' + groupID
				os.system('mkdir ' + groupDIR)
				#os.system('e2proc3d.py ' + entirestack + ' ' + groupPATH + ' --first=' + str(bottom_range) + ' --last=' + str(top_range)) 	
				
				mm = 0
				for jj in xrange(bottom_range,top_range):
					a = EMData(entirestack,jj)
					a['spt_ptcl_indxs'] = mm
					a.write_image(groupPATH,mm)
					mm += 1
				
				options.input = groupPATH
				options.path = groupDIR

		print "I will start ALL vs ALL on group number", i+1
		print "For which options.input is", options.input
		allvsall(options)
	return()
	
	
def allvsall(options):
	
	print "Therefore the received options.input inside ALL VS ALL is", options.input
	print "With these many particles in it", EMUtil.get_image_count(options.input)
	hdr = EMData(options.input,0,True)
	nx = hdr["nx"]
	ny = hdr["ny"]
	nz = hdr["nz"]
	if nx!=ny or ny!=nz :
		print "ERROR, input volumes are not cubes"
		sys.exit(1)
	
	nptcl = EMUtil.get_image_count(options.input)
	if nptcl<3: 
		print "ERROR: at least 3 particles are required in the input stack for all vs all. Otherwise, to align 2 particles (one to the other or to a model) use e2classaverage3d.py"
		sys.exit(1)
	
	fillfactor = len(str(nptcl))							#Calculate this based on the number of particles so that tags are adequate ("pretty") and ordered
	roundtag='round' + str(0).zfill(fillfactor)					#We need to keep track of what round we're in
	newptcls={}									#This dictionary stores 'new particles' produced in each round as { particleID : particleDATA } elements
	allptclsRound={}								#This dictionary stores all particlces in a round ("new" and "old") as 
											#{particle_id : [EMData,{index1:totalTransform1, index2:totalTransform2...}]} elements
											#The totalTransform needs to be calculated for each particle after each round, to avoid multiple interpolations
	
	for i in range(nptcl):								#In the first round, all the particles in the input stack are "new" and should have an identity transform associated to them
		a=EMData(options.input,i)
		totalt=Transform()

		if 'spt_multiplicity' not in a.get_attr_dict():				#spt_multiplicity keeps track of how many particles were averaged to make any given new particle (set to 1 for the raw data)
			a['spt_multiplicity']=1
		if 'spt_ptcl_indxs' not in a.get_attr_dict():				#spt_ptcl_indxs keeps track of what particles from the original stack went into a particular average or "new particle"
			a['spt_ptcl_indxs']=[i]						#The raw stack should contain particles where this parameter is the particle number itself
		else:
			if type(a['spt_ptcl_indxs']) is int:
				a['spt_ptcl_indxs'] = [a['spt_ptcl_indxs']]		#The spt_ptcl_indxs parameter should be of type 'list', to easily 'append' new particle indexes
		
		a.write_image(options.input,i)						#Overwrite the raw stack with one that has the appropriate header parameters set to work with e2sptallvsall
		
		particletag = roundtag + '_' + str(i).zfill(fillfactor)
		newptcls.update({particletag :a})
		
		if 'sptID' not in a.get_attr_dict():					#spt_multiplicity keeps track of how many particles were averaged to make any given new particle (set to 1 for the raw data)
			a['sptID'] = particletag
		allptclsRound.update({particletag : [a,{i:totalt}]})			
		
	oldptcls = {}									#'Unused' particles (those that weren't part of any unique-best-pair) join the 'oldptcls' dictionary onto the next round
	surviving_results = []								#This list stores the results for previous alignment pairs that weren't used, so you don't have to recompute them
		
	allptclsMatrix = []								#Massive matrix listing all the allptclsRound dictionaries produced after each iteration
	allptclsMatrix.append(allptclsRound)

	for k in range(options.iter):							#Start the loop over the user-defined number of iterations
		#avgname = options.path + '/round' + str(k).zfill(fillfactor) + '_averages.hdf'
		newstack = options.path + '/round' + str(k-1).zfill(fillfactor) + '_averages.hdf'
		if k== 0:
			newstack =options.input
		
		nnew = len(newptcls)
		if nnew + len(oldptcls) == 1:						#Stop the loop if the data has converged and you're left with one final particle (an average of all)
				print "The all vs all algorithm has converged into one average"
				break
		
		allptclsRound = {}							
		
		logger = E2init(sys.argv,options.ppid)
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
		
		for ptcl1, compare in newptclsmap:
			for ptcl2 in compare:
				
				reftag = roundtag + str(ptcl1).zfill(fillfactor)								
				particletag = roundtag + str(ptcl2).zfill(fillfactor)
				
				#if options.verbose > 2:
				print "Setting the following comparison: %s vs %s in ALL VS ALL" %(reftag,particletag)
				
				task = Align3DTaskAVSA(newstack,newstack, jj, reftag, particletag, ptcl1, ptcl2,"Aligning particle#%s VS particle#%s in iteration %d" % (reftag,particletag,k),options.mask,options.normproc,options.preprocess,
				options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
				
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
				oldtags.update({EMData(options.path + '/oldptclstack.hdf',i,True)['sptID'] : i})
			
			
			print "Old tagas are:\n", oldtags
			nnn = 0
			for refkey,refvalue in newptcls.iteritems():
				ptcl1 = nnn
				for particlekey,particlevalue in oldptcls.iteritems():
					
					ptcl2 = oldtags[particlekey]
					
					#if options.verbose > 2:
					print "Setting the following comparison: %s vs %s in ALL VS ALL" %(refkey,particlekey)
					
					task = Align3DTaskAVSA(newstack,options.path + '/oldptclstack.hdf',jj , refkey, particlekey, ptcl1, ptcl2,"Aligning particle round#%d_%d VS particle#%s, in iteration %d" % (k,ptcl1,particlekey.split('_')[0] + str(ptcl2),k),options.mask,options.normproc,options.preprocess,
					options.npeakstorefine,options.align,options.aligncmp,options.ralign,options.raligncmp,options.shrink,options.shrinkrefine,options.verbose-1)
					
					tasks.append(task)
										
					jj+=1
				nnn+=1	
		
		tids=etc.send_tasks(tasks)						#Start the alignments running
		#if options.verbose > 0: 
		print "%d tasks queued in iteration %d"%(len(tids),k) 
		
		results = get_results(etc,tids,options.verbose)				#Wait for alignments to finish and get results
		results = results + surviving_results					#The total results to process/analyze includes results (comparisons) from previous rounds that were not used
		results = sorted(results, key=itemgetter('score'))			#Sort/rank the results by score
		
		if options.verbose > 0:
			print "In iteration %d the SORTED results are:", k
			for i in results:
				print "%s VS %s , score=%f" %(i['ptcl1'], i['ptcl2'], i['score'])
		
		print "\n\n\n\nIn iteration %d, the total number of comparisons in the ranking list, either new or old that survived, is %d" % (k, len(results))
		
		tried = set()											#Tracks what particles have "appeared" on the list, whether averaged or not
		averages = {}											#Stores the new averages; you need a dict different from newptcls because you'll need to 'fetch' data from the previous version of newptcls
		used = set()											#Tracks what particles WERE actually averaged
		
		mm=0												#Counter to track new particles/averages produced and write them to output
		for z in range(len(results)):
			if results[z]['ptcl1'] not in tried and results[z]['ptcl2'] not in tried:
				tried.add(results[z]['ptcl1'])							#If the two particles in the pair have not been tried, and they're the next "best pair", they MUST be averaged
				tried.add(results[z]['ptcl2'])							#Add both to "tried" AND "averages" 
				used.add(results[z]['ptcl1'])		
				used.add(results[z]['ptcl2'])
													
				avgr = Averagers.get(options.averager[0], options.averager[1])			#Call the averager
				
				avg_ptcls = []
								
				ptcl1 = allptclsMatrix[k][results[z]['ptcl1']][0]
							
														#You always add all the past particles that went into a particular particle (being
														#averaged in the current round) freshly from the raw stack to the averager (with
														#the appropriate transforms they've undergone, of course. Thus, YOU DON'T have to
														#worry about "multiplicity", since it takes care of itself by doing this.
				indx_trans_pairs = {}

				print "The indexes in particle 1 are", ptcl1['spt_ptcl_indxs']
				
				ptcl1info = allptclsMatrix[k][results[z]['ptcl1']]
				
				ptcl1_indxs_transforms = ptcl1info[-1]
								
				for p in ptcl1['spt_ptcl_indxs']:											
					pastt = ptcl1_indxs_transforms[p]
					
					subp1 = EMData(options.input,p)
					subp1.process_inplace("xform",{"transform":pastt})
					
					avgr.add_image(subp1)
					
					if options.saveali:
						avg_ptcls.append(subp1)
					
					indx_trans_pairs.update({p:pastt})
					
						
				#avgr.add_image(ptcl1)								#Add particle 1 to the average
				
				ptcl2 = allptclsMatrix[k][results[z]['ptcl2']][0]
						
				resultingt = results[z]["xform.align3d"]					
					
				print "The indexes in particle 2 are", ptcl2['spt_ptcl_indxs']
				
				ptcl2info = allptclsMatrix[k][results[z]['ptcl2']]
					
				ptcl2_indxs_transforms = ptcl2info[-1]
								
				#ptcl2avgr = Averagers.get(options.averager[0], options.averager[1])		#You need to recompute ptcl2 "fresh" from the raw data to avoid multiple interpolations
				
				for p in ptcl2['spt_ptcl_indxs']:						#All the particles in ptcl2's history need to undergo the new transformation before averaging
					#print "I'm fixing the transform for this index", p			#(multiplied by any old transforms, all in one step, to avoid multiple interpolations)
					
					pastt = ptcl2_indxs_transforms[p]
					
					totalt = resultingt * pastt
					subp2 = EMData(options.input,p)
					subp2.process_inplace("xform",{"transform":totalt})
					
					avgr.add_image(subp2)
					
					if options.saveali:
						avg_ptcls.append(subp2)					
					
					indx_trans_pairs.update({p:totalt})
					
							
				avg=avgr.finish()
				
				if options.autocenter:
					print "\n\n\n\nYou have selected to autocenter!\n"
					avg = avg.process('xform.centerofmass')
					tcenter = avg['xform.align3d']
					print "Thus the average will be translated like this", tcenter
					
					if options.saveali:
						avg_ptcls = []
					
					for p in ptcl1['spt_ptcl_indxs']:
						pastt = ptcl1_indxs_transforms[p]
						totalt = tcenter * pastt
						indx_trans_pairs.update({p:totalt})
						
						if options.saveali:
							subp1 = EMData(options.input,p)
							subp1.process_inplace("xform",{"transform":totalt})
							avg_ptcls.append(subp1)
						
					for p in ptcl2['spt_ptcl_indxs']:
						pastt = ptcl2_indxs_transforms[p]
						totalt = tcenter * pastt
						indx_trans_pairs.update({p:totalt})
						
						if options.saveali:
							subp2 = EMData(options.input,p)
							subp2.process_inplace("xform",{"transform":totalt})
							avg_ptcls.append(subp2)
						
				
				avgmultiplicity = ptcl1['spt_multiplicity'] + ptcl2['spt_multiplicity']		#Define and set the multiplicity of the average
				avg['spt_multiplicity'] = avgmultiplicity
				
				indexes1 = ptcl1["spt_ptcl_indxs"]
				indexes2 = ptcl2["spt_ptcl_indxs"]				
				
				avg["spt_ptcl_indxs"] = indexes1 + indexes2					#Keep track of what particles go into each average or "new particle"				
				
				avg["spt_ptcl_src"] = options.input
				
				avg['origin_x'] = 0								#The origin needs to be set to ZERO to avoid display issues in Chimera
				avg['origin_y'] = 0
				avg['origin_z'] = 0
								
				newroundtag = 'round' + str(k+1).zfill(fillfactor) + '_'
				avgtag = newroundtag + str(mm).zfill(fillfactor)
				
				avg['sptID'] = avgtag
				
				avg.process_inplace("normalize.edgemean")
				
				avg.write_image(options.path + '/round' + str(k).zfill(fillfactor) + '_averages.hdf',mm)
				
				if options.saveali:
					for oo in range(len(avg_ptcls)):
						avg_ptcls[oo].write_image(options.path + '/round' + str(k).zfill(fillfactor) + '_average' + str(mm).zfill(2)  + '_ptcls.hdf',oo)
				
				if options.postprocess!=None : 
					avgp=avg.process(options.postprocess[0],options.postprocess[1])
					avgp.write_image(options.path + '/round' + str(k).zfill(fillfactor) + '_averages_postp.hdf',mm)
										
				averages.update({avgtag:avg})	   						#The list of averages will become the new set of "newptcls"
				allptclsRound.update({avgtag : [avg,indx_trans_pairs]})
				
				mm+=1
				
			if results[z]['ptcl1'] not in tried:						#If a particle appeared in the ranking list but its pair was already taken, the particle must be classified as "tried"
				tried.add(results[z]['ptcl1'])						#because you don't want to average it with any other available particle lower down the list that is available
													#We only average "UNIQUE BEST PAIRS" (the first occurance in the ranking list of BOTH particles in a pair).
			if results[z]['ptcl2'] not in tried:
				tried.add(results[z]['ptcl2'])
		
		surviving_results = []
		for z in range(len(results)):
			if results[z]['ptcl1'] not in used and results[z]['ptcl2'] not in used:
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
		
		os.system('rm ' + options.path + '/oldptclstack.hdf')
		
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
				comparison=r[0].options["comparison"]					# get the comparison number from the task rather than trying to work back to it
				
				results[comparison]=r[1]["final"][0]					# this will be a list of (qual,Transform), containing the BEST peak ONLY
				
				results[comparison]['ptcl1']=r[0].options['ptcl1']			#Associate the result with the pair of particles involved
				results[comparison]['ptcl2']=r[0].options['ptcl2']

				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		
		print "  %d tasks, %d complete, %d waiting to start        \r"%(len(tids),ncomplete,nwait)

		if verbose:
			sys.stdout.flush()
	
		if len(tidsleft)==0 or ncomplete == numtides: 
			break		
	return results


class Align3DTaskAVSA(EMTask):
	"""This is a task object for the parallelism system. It is responsible for aligning one 3-D volume to another, with a variety of options"""

	def __init__(self,fixedimagestack,imagestack,comparison,ptcl1,ptcl2,p1n,p2n,label,mask,normproc,preprocess,npeakstorefine,align,aligncmp,ralign,raligncmp,shrink,shrinkrefine,verbose):
		"""fixedimage and image may be actual EMData objects, or ["cache",path,number]
	label is a descriptive string, not actually used in processing
	ptcl is not used in executing the task, but is for reference
	other parameters match command-line options from e2classaverage3d.py
	Rather than being a string specifying an aligner, 'align' may be passed in as a Transform object, representing a starting orientation for refinement"""
		data={}
		data={"fixedimage":fixedimagestack,"image":imagestack}
		EMTask.__init__(self,"ClassAv3d",data,{},"")

		self.options={"comparison":comparison,"ptcl1":ptcl1,"ptcl2":ptcl2,"p1number":p1n,"p2number":p2n,"label":label,"mask":mask,"normproc":normproc,"preprocess":preprocess,"npeakstorefine":npeakstorefine,"align":align,"aligncmp":aligncmp,"ralign":ralign,"raligncmp":raligncmp,"shrink":shrink,"shrinkrefine":shrinkrefine,"verbose":verbose}
	
	def execute(self,callback=None):
		"""This aligns one volume to a reference and returns the alignment parameters"""
		options=self.options
		if options["verbose"]>1: 
			print "Aligning ",options["label"]
		
		fixedimage=EMData(self.data["fixedimage"],options['p1number'])
		fn=EMUtil.get_image_count(self.data["fixedimage"])
		
		#if type(self.data) != libpyEMData2.EMData:
		#	print 
		
		image=EMData(self.data["image"],options['p2number'])
		iin=EMUtil.get_image_count(self.data["image"])
		
		mask=EMData(int(image['nx']),int(image['ny']),int(image['nz']))
		mask.to_one()
		
		if options["mask"] != None:
			mask.process_inplace(options["mask"][0],options["mask"][1])
		
		# normalize
		if options["normproc"] != None:
			if options["normproc"][0]=="normalize.mask": 
				options["normproc"][1]["mask"]=mask
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
						 
		if options["verbose"] >2:
			print "Because it was greater than 2 or not integer, I will exit"
			sys.exit()  
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
		
		if bestfinal[0]["score"] == 1.0e10 :
			print "Error: all refine alignments failed for %s. May need to consider altering filter/shrink parameters. Using coarse alignment, but results are likely invalid."%self.options["label"]
		
		if options["verbose"]>1: 
			print "Best %1.5g\t %s"%(bestfinal[0]["score"],str(bestfinal[0]["xform.align3d"])) 
			print "Done aligning ",options["label"]
		
		return {"final":bestfinal,"coarse":bestcoarse}

if __name__ == '__main__':
	main()



