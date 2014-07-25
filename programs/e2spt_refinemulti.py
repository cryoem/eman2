#!/usr/bin/env python

#
# Author: Jesus Galaz-Montoya  November/21/2013
# Copyright (c) 2011- Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA

from sys import argv
import os
from EMAN2 import *

from e2spt_classaverage import sptmakepath

import subprocess

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <stack>

	WARNING: This program still EXPERIMENTAL (It's under heavy development)	
	
	Refinement of a 3D-volume stack against multiple models. The initial models can be provided in a stack, OR generated from the data itself.
	
	When no reference is provided, you define a number of models greater than 1 to use (for example, 2, 3, 4 or more).
	[If you want to refine the data against ONE model, use e2spt_refine.py]
	The data set is divded into that specified number of groups.
	
	An initial model will be generated with the particles assigned to each group.
	Then, the entire data set will be refined against all initial models.
	
	You can increase the number of references used for each iteration by specifying the 
	--addmodel parameter.
	This will take the "best initial model" (the one that most particles preferred) and include it as an initial model for the next round of refinement.
	For exampe, if you start with two references A and B, 
	two averages will come out of aligning the data against them, A' and B'.
	So if --addmodel is on, instead of only using A' and B' as references for the next 
	refinement round, the best of A and B will also be used,
	which means you will refine the data against 3 models in the next round, not just 2.
	
	If you supply a single reference/model then --addmodel will be ASSUMED to be True; 
	otherwise, to refine a data set against a single model use
	e2spt_refine.py
	 """

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	#parser.add_argument("--ncls", type=int, help="...", default=2)
	#parser.add_argument("--nbasis", type=int, help="Basis vectors to use", default=3)

	parser.add_argument("--refs", type=str, help="""This can either be an HDF stack, where each image will be treated as a separate model/reference, 
													or a comma separatted list of individual images; e.g. --refs=ref1.hdf,ref2.hdf,ref3.hdf.
													""", default='')
	
	parser.add_argument("--nrefs", type=int, help="""Number of references to generate from the data for reference-free alignment. Default=1""", default=1)
	parser.add_argument("--refsgenmethod", type=str, help="""Method for generating the initial reference(s). Options are 'binarytree' and 'hac'. Default=binarytree""", default='binarytree') 
	#parser.add_argument("--refpreprocess",action="store_true",default=False,help="""This 
	#	will preprocess the reference identically to the particles. It is off by default, but it is internally turned on when no reference is supplied.""")
	
	
	'''
	PARAMETERS TO BE PASSED ON TO e2spt_classaverage.py
	'''
	
	parser.add_header(name="caheader", help='Options below this label are specific to e2spt_classaverage', title="### e2spt_classaverage options ###", default=None, guitype='filebox', row=3, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--path",type=str,default='',help="""Directory to store results in. 
		The default is a numbered series of directories containing the prefix 'spt_refinemulti'; 
		for example, spt_refinemulti02 will be the directory by default if 'spt_refinemulti01' 
		already exists.""")
	
	parser.add_argument("--syms", type=str, help="""List comma-separated symmetries to apply
		separately on the different references. For example, if you provide --syms=d8,d7
		and provide 2 references via --nrefs=2 or supply two references via --refs=r1.hdf,r2.hdf, 
		d8 symmetry will be applied to the first reference and d7 to the second after each iteration
		of refinement (the final average in one iteration becomes a reference for the next).""", default='')
	
	parser.add_argument("--input", type=str, help="The name of the input volume stack. MUST be HDF or BDB, since volume stack support is required.", default=None, guitype='filebox', browser='EMSubTomosTable(withmodal=True,multiselect=False)', row=0, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--output", type=str, help="The name of the output class-average stack. MUST be in  .hdf format, since volume stack support is required.", default=None, guitype='strbox', row=2, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#parser.add_argument("--oneclass", type=int, help="Create only a single class-average. Specify the class number.",default=None)
	#parser.add_argument("--classmx", type=str, help="The name of the classification matrix specifying how particles in 'input' should be grouped. If omitted, all particles will be averaged.", default='')
	
	#parser.add_argument("--ref", type=str, help="Reference image(s). Used as an initial alignment reference and for final orientation adjustment if present. This is typically the projections that were used for classification.", default=None, guitype='filebox', browser='EMBrowserWidget(withmodal=True,multiselect=True)', filecheck=False, row=1, col=0, rowspan=1, colspan=3, mode='alignment')
	
	#parser.add_argument("--resultmx",type=str,help="Specify an output image to store the result matrix. This is in the same format as the classification matrix. http://blake.bcm.edu/emanwiki/EMAN2/ClassmxFiles", default=None)
	
	#parser.add_argument("--refinemultireftag", type=str, help="DO NOT USE THIS PARAMETER. It is passed on from e2spt_refinemulti.py if needed.", default='')
	
	
	parser.add_argument("--radius", type=float, help="""Will make --align and --falign None. Hydrodynamic radius of the particle in Angstroms. 
													This will be used to automatically calculate the angular steps to use in search of the best alignment.
													Make sure the apix is correct on the particles' headers, sine the radius will be converted from Angstroms to pixels.
													Then, the fine angular step is equal to 360/(2*pi*radius), and the coarse angular step 4 times that""", default=0)
	
	parser.add_argument("--iter", type=int, help="The number of iterations to perform. Default is 1.", default=1, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, nosharedb=True, mode='alignment,breaksym')
	parser.add_argument("--savesteps",action="store_true", help="If set, will save the average after each iteration to class_#.hdf. Each class in a separate file. Appends to existing files.",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='alignment,breaksym')
	parser.add_argument("--saveali",action="store_true", help="If set, will save the aligned particle volumes in class_ptcl.hdf. Overwrites existing file.",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')
	parser.add_argument("--saveallalign",action="store_true", help="If set, will save the alignment parameters after each iteration",default=False, guitype='boolbox', row=4, col=2, rowspan=1, colspan=1, mode='alignment,breaksym')
	parser.add_argument("--sym", dest = "sym", default=None, help = "Symmetry to impose - choices are: c<n>, d<n>, h<n>, tet, oct, icos", guitype='symbox', row=9, col=1, rowspan=1, colspan=2, mode='alignment,breaksym')
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", returnNone=True, default="mask.sharp:outer_radius=-2", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'mask\')', row=11, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize")
	
	
	parser.add_argument("--threshold",type=str,help="""A threshold applied to the subvolumes after normalization. 
													For example, --threshold=threshold.belowtozero:minval=0 makes all negative pixels equal 0, so that they do not contribute to the correlation score.""", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	
	
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--preprocessfine",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.", default=None)
	
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=17, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--lowpassfine",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=18, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--highpassfine",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--postprocess",type=str,help="A processor to be applied to the FINAL volume after averaging the raw volumes in their FINAL orientations, after all iterations are done.",default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=16, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--procfinelikecoarse",type=bool,default=True,help="Turn on with --procfinelikecoarse=False, and supply fine alignment parameters, such as --lowpassfine, --highpassfine, etc; to preprocess the particles for FINE alignment differently than for COARSE alignment.")
	
	parser.add_argument("--ncoarse", type=int, help="Deprecated. Use --npeakstorefine instead.", default=None)
	parser.add_argument("--npeakstorefine", type=int, help="The number of best coarse alignments to refine in search of the best final alignment. Default=4.", default=4, guitype='intbox', row=9, col=0, rowspan=1, colspan=1, nosharedb=True, mode='alignment,breaksym[1]')
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average. Default is rotate_translate_3d:search=10:delta=15:dphi=15, specify 'None' to disable", returnNone=True, default="rotate_translate_3d:search=10:delta=15:dphi=15", guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'3d\')', row=12, col=0, rowspan=1, colspan=3, nosharedb=True, mode="alignment,breaksym['rotate_symmetry_3d']")
	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo", guitype='comboparambox',choicelist='re_filter_list(dump_cmps_list(),\'tomo\')', row=13, col=0, rowspan=1, colspan=3,mode="alignment,breaksym")
	parser.add_argument("--falign",type=str,help="This is the second stage aligner used to refine the first alignment. Default is refine.3d, specify 'None' to disable", default="refine_3d", returnNone=True, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine.*3d\')', row=14, col=0, rowspan=1, colspan=3, nosharedb=True, mode='alignment,breaksym[None]')
	parser.add_argument("--faligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo", guitype='comboparambox',choicelist='re_filter_list(dump_cmps_list(),\'tomo\')', row=15, col=0, rowspan=1, colspan=3,mode="alignment,breaksym")
	parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean")
	parser.add_argument("--keep",type=float,help="The fraction of particles to keep in each class.",default=1.0, guitype='floatbox', row=6, col=0, rowspan=1, colspan=1, mode='alignment,breaksym')

	parser.add_argument("--inixforms",type=str,help="directory containing a dict of transform to apply before reference generation", default="", guitype='dirbox', dirbasename='spt_|sptsym_', row=7, col=0,rowspan=1, colspan=2, nosharedb=True, mode='breaksym')

	parser.add_argument("--breaksym",action="store_true", help="Break symmetry. Do not apply symmetrization after averaging", default=False, guitype='boolbox', row=7, col=2, rowspan=1, colspan=1, nosharedb=True, mode=',breaksym[True]')
	
#	parser.add_argument("--groups",type=int,help="WARNING: This parameter is EXPERIMENTAL, and will only work if --iter=1. It's the number of final averages you want from the set after ONE iteration of alignment. Particles will be separated in groups based on their correlation to the reference",default=0)

	parser.add_argument("--randomizewedge",action="store_true", help="This parameter is EXPERIMENTAL. It randomizes the position of the particles BEFORE alignment, to minimize missing wedge bias and artifacts during symmetric alignment where only a fraction of space is scanned", default=False)
	parser.add_argument("--savepreprocessed",action="store_true", help="Will save stacks of preprocessed particles (one for coarse alignment and one for fine alignment if preprocessing options are different).", default=False)
	parser.add_argument("--keepsig", action="store_true", help="Causes the keep argument to be interpreted in standard deviations.",default=False, guitype='boolbox', row=6, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')
	parser.add_argument("--nocenterofmass", action="store_true", help="Disable Centering of mass of the subtomogram every iteration.", default=False, guitype='boolbox', row=6, col=2, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	#parser.add_argument('--reverse_contrast', action="store_true", default=False, help=""" This multiplies the input particles by -1. Remember that EMAN2 **MUST** work with 'white protein' """)
	
	parser.add_argument("--shrink", type=int,default=0,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	parser.add_argument("--shrinkfine", type=int,default=0,help="Optionally shrink the input volumes by an integer amount for fine alignment.")
	
	#parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default='', guitype='strbox', row=19, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')

	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1", guitype='strbox', row=19, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#parser.add_argument("--automask",action="store_true",help="Applies a 3-D automask before centering. Can help with negative stain data, and other cases where centering is poor.")
	#parser.add_argument("--resample",action="store_true",help="If set, will perform bootstrap resampling on the particle data for use in making variance maps.",default=False)
	#parser.add_argument("--odd", default=False, help="Used by EMAN2 when running eotests. Includes only odd numbered particles in class averages.", action="store_true")
	#parser.add_argument("--even", default=False, help="Used by EMAN2 when running eotests. Includes only even numbered particles in class averages.", action="store_true")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	parser.add_argument("--plotccc", action='store_true', help="""Turn this option on to generate
		a plot of the ccc scores both during model generation with e2spt_classaverage.py or
		e2spt_hac.py and for refinement results of e2spt_refinemulti.py.
		Running on a cluster or via ssh remotely might not support plotting.""",default=False)
	
	(options, args) = parser.parse_args()

	try:
		hdr = EMData(options.input,0,True) #This is done just to test whether the input file exists where it should
		boxsize = hdr['nx']
	except:
		print "ERROR: Can't find the file provided through --input"

	if options.radius:
		options.align = None
		options.falign = None
	
		
	#print "\n\n\n(e2spt_refinemulti.py) options.refpreprocess is", options.refpreprocess
	#print "\n\n\n"
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	
	print "(e2spt_refinemulti.py) BEFORE sptmakepath, otions.path is", options.path
	
	options = sptmakepath(options,'spt_refinemulti')

	print "(e2spt_refinemulti.py) AFTER sptmakepath, otions.path is", options.path

	rootpath = os.getcwd()
	
	originalCompletePath = rootpath + '/' + options.path
		
	if '/' not in options.input:
		#relativeInput = '../' + options.input
		#absoluteInput = rootpath + '/' + options.input
		#options.input = absoluteInput
		originalCompleteStack = rootpath + '/' + options.input
		options.input = originalCompleteStack
	
	avgs={}
	finalize=0

	
	
	'''
	Store parameters in parameters.txt file inside --path
	'''
	from e2spt_classaverage import writeParameters
	writeParameters(options,'e2spt_refinemulti.py', 'refinemulti')
	
	'''
	Determine how many references there are and put them into one file classAvg.hdf, or 
	classAvg_iterXX.hdf, if they come from separate files
	'''
	nrefs=0
	
	refsfiles = set([])
	
	if not options.refs:
		'''
		If no references are provided, the program has to generate them, using the method
		specified through --refsgenmethod
		'''
		
		print "(e2spt_refinemulti.py) No references provided; therefore, genrefs function will be called"
		ret = genrefs( options, originalCompletePath)
		
		print "\n(e2spt_refinemulti.py) Back to main, genrefs has returned", ret
		
		for rf in ret:
			actualFyle = os.path.basename(rf)
			actualFyleFixed = actualFyle.replace('group','ref').replace('avg','')
			rf2add = originalCompletePath + '/' + actualFyleFixed
			os.system('cp ' + rf + ' ' + rf2add)
			
			refsfiles.update( [rf2add] )
		
		
		
		
	else:
		refsFylesOrig = options.refs.split(',')
		
		#if ',' in options.refs:
		#refsorig = options.refs.split(',')
		
		if len(refsFylesOrig) > 1:
			nrefs = len(refsFylesOrig)
			
			i=0
			for refFyle in refsFylesOrig:
				outref = originalCompletePath + '/' + os.path.basename(refFyle).replace('.hdf','_ref' + str(i).zfill(len(str(nrefs))) + '.hdf' )
				os.system('cp ' + refFyle + ' ' + outref)		
				#nrefs = len(refsorig)
		
				refsfiles.update( [ outref ] )
				i+=1
		else:
			nrefs = EMUtil.get_image_count( options.refs )
			for i in range(nrefs):
				outref = originalCompletePath + '/' + options.refs.replace('.hdf','_ref' + str(i).zfill(len(str(nrefs))) + '.hdf' )
				os.system('e2proc3d.py ' + options.refs + ' ' + outref + ' --first=' + str(i) + ' --last=' + str(i) )
				refsfiles.update( [ outref ] )


	#print "ERROR: You must provide at least one reference through --ref, or specify the number of references to generate from the data through --nrefs."
		
		
	'''
	Generate the commands to refine the data against each reference
	'''		
	
	for it in range( options.iter ):
		print "\n\nIteration", it
		print "\n\n"
				#The program goes into --path to execute the alignment command; therefore, --input will be one level furtherback
	
		if it > 0:
			newrefsfiles=set()				
			#avgsName =  'classAvgs.hdf'
			#if options.savesteps and int( options.iter ) > 1 :
			#	avgsName = 'classAvgs_iter' + str().zfill( len( str( options.iter ))) + '.hdf'
	
			newreffile = ''
			for reftag in avgs:
				if avgs[ reftag ] and avgs[reftag] != None and avgs[reftag]!='None':
					newref = avgs[ reftag ]
					
					for ref in refsfiles:
						if reftag in ref:
							if '_iter' in ref:
								newreffile = ref.split('_iter')[0] + '_iter' + str(it).zfill( len( str( options.iter ))) + '.hdf'
							else:
								newreffile = ref.replace('.hdf','_iter' + str(it).zfill( len( str( options.iter ))) + '.hdf')
						
							newref.write_image( newreffile, 0)
				#else:
				#	for ref in refsfiles:
				#		if reftag in ref:					
				#			newreffile = ref
				
				newrefsfiles.update( [newreffile] )
			
			refsfiles = newrefsfiles
			
			if len(refsfiles) < 2:
				finalize = 1
				print "(e2spt_refinemulti.py, line 268) All particles preferred one average and therefore multirefine has failed/converged"
				sys.exit()
				
		
		#print "\n\n\n\nRRRRRRRRRRRRRR \n BEFORE alignment loop for all refts, len refsfiles is and refsfiles are", len(refsfiles), refsfiles
		#print "\nRRRRRRRRRRRRRRRRR \n\n\n\n\n"
		
		k=0
		reftags = []
		masterInfo = {}
		
		if options.syms and it == 0:
			options.syms = options.syms.split(',')
		
			if len(options.syms) != len(refsfiles):
				if len(options.syms) > len(refsfiles):
					options.syms = options.syms[0,len(refsfiles)]
				elif len(options.syms) < len(resfiles):
				
					howMany = len(refsfiles) - len(options.syms)
				
					for pi in range(howMany):
						options.syms.append('c1')
			
		for ref in refsfiles:
			#print "\n\n\n\n\n\n\n\n\nAAAAAAAAAAAAAAAAAAAAAAAA\nAligning data to ref,refnumber", ref,k
			#print "\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n\n\n\n\n"
			#ref = ref
			
			#thisRefinementPath = ref.split('/')[-1].replace('.hdf','')
			
			sym=''
			if options.syms:
				sym = options.syms[k]
			
			alicmd ='cd ' + originalCompletePath + ' && e2spt_classaverage.py --ref=' + ref + ' --path=sptTMP' #+ ' --input=' + relativeInput
	
			names = dir(options)
			for name in names:
				if getattr(options,name) and 'refs' not in name and 'syms' not in name and "__" not in name and "_" not in name and 'path' not in name and str(getattr(options,name)) != 'True' and 'iter' not in name:	
					#if "__" not in name and "_" not in name and str(getattr(options,name)) and 'path' not in name and str(getattr(options,name)) != 'False' and str(getattr(options,name)) != 'True' and str(getattr(options,name)) != 'None':			
					alicmd += ' --' + name + '=' + str(getattr(options,name))
			alicmd += ' --donotaverage'
			
			if options.sym:
				alicmd += ' --sym=' + sym
			
			#if options.refpreprocess:
			alicmd += ' --refpreprocess'
				
			#tag = str(k).zfill( len( str ( nrefs )))
			
			tag = os.path.basename(ref).split('ref')[-1].split('_')[0].replace('.hdf','').zfill( len( str ( len(refsfiles) )))
			
			reftag = 'ref' + tag
			reftags.append(reftag)
		
			jsAliParamsPathActual = originalCompletePath + '/tomo_xforms.json'
			
			jsAliParamsPathNew = jsAliParamsPathActual.replace('.json','_' + reftag + '.json')
			jsAliParamsPathNew = jsAliParamsPathNew.replace('.json', '_it' + str(it).zfill( len(str(options.iter))) + '.json')
			
			#print "(e2spt_refinemulti.py) The actual .json file with ali params will be", jsAliParamsPathNew
			
			alicmd += ' --refinemultireftag=' + tag	+ ' && mv ' + originalCompletePath + '/sptTMP/* ' + originalCompletePath + '/ && rm -r ' + originalCompletePath +'/sptTMP* && cp ' + jsAliParamsPathActual + ' ' + jsAliParamsPathNew
		
		
			
			#alicmd += ' --refinemultireftag=' + tag	+ ' && cp ' + options.path + '/sptTMP/* ' + options.path + '/ && rm -r ' + options.path +'/sptTMP* && mv ' + jsAliParamsPathActual + ' ' + jsAliParamsPathNew

		
			#print "Command is", alicmd
		
			'''
			Make sure the supbrocess that executes e2spt_classaverage.py ends before script continues
			'''
			
			#print "The command to execute is", alicmd
			
			p=subprocess.Popen( alicmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		
			for line in iter(p.stdout.readline, ''):
				print line.replace('\n','')
		
			#'''
			#Option B
			#'''
			#while True:
			#	line = p.stdout.readline()
			#	if not line: 
			#		break
			#	else:
			#		print line.replace('\n','')
		
			#returnedtxt = p.communicate()
			text=p.communicate()	
			p.stdout.close()
		
		
			print "\n\n\nThis is the alicmd", alicmd
			print "\ntomo_xforms.json should be in direct, see", os.listdir(originalCompletePath)
			print"\n\n\n"
			
			print "Feedback from p was", text
		
		
			'''
			Open alignment results file for current reference
			'''
			#aliParamsFile = options.path + '/' + thisRefinementPath + '/tomo_xforms' + tag + '.json'
			#aliParamsFile = options.path + '/tomo_xforms' + tag + '.json'
			aliParams = js_open_dict(jsAliParamsPathNew)
			nparams = len(aliParams)
		
			#print "(e2spt_refinemulti.py) The aliParams file to read is", jsAliParamsPathNew
			#print "(e2spt_refinemulti.py) I read these many params", nparams
			#print "Which are", aliParams
		
			'''
			Add info per particle for results from all references to a master dictionary, 'masterInfo'
			'''
			for i in range(nparams):
				ptclID = "tomo_" + str(i).zfill( len(str( nparams )) )
				ptclScore = float( aliParams[ptclID][-1] )
				ptclAliParams = aliParams[ptclID]
			
				infolist = [ ptclScore, ptclAliParams, reftag]
			
				if k==0:
					masterInfo.update({ ptclID: [] })
			
				#print "\n\nptclID to update is", ptclID
				#print "infolist to append is", infolist
				#print "BEFORE appending, masterInfo[ ptclID ] is", masterInfo[ ptclID ]
				#print "Of type", type(masterInfo[ ptclID ])
			
				value = masterInfo[ ptclID ]
				value.append(infolist)
			
				#print "Therfore value is", value
				#print "Of type", type(value)
				#print "\n\n"
			
				masterInfo.update({ ptclID: value })
				
				
				'''
				Open scores file for current reference
				'''
				scoresFile = originalCompletePath +'/subtomo_FinalScores_' + reftag + '.json'
				jsScores = js_open_dict(scoresFile)
				
				jsScores.setval(ptclID, 'score='+str(ptclScore) + ' vs ref=' + ref  )
				
				jsScores.close()
				
			k+=1
		
		'''
		Analyze all results and classify particles based on them
		'''			
		from operator import itemgetter						
	
		print "I've aligned all the particles in the data set to all the references for iter %d and will now classify them from the masterInfo dict" %(it), masterInfo
	
		classes = {}
		for reftag in reftags:
			print "\n\n\n\n\n\n\n\n\nRRRRRRRRRR\nreftag is", reftag
			classes.update({ reftag : [] })
	
		for ele in masterInfo:
			sortedPtclInfo = sorted( masterInfo[ele], key=itemgetter(0))	#Sorted works because you want the scores from SMALLEST to BIGGEST. Remember, the MORE NEGATIVE (smaller) the better score in EMAN2
			bestPtclInfo = sortedPtclInfo[0]
		
			bestreftag = bestPtclInfo[-1]	
			bestAliParams = bestPtclInfo[1]
			bestScore = bestPtclInfo[0]
		
			ptclIndx = int( ele.split('_')[-1] )
			value = classes[ bestreftag ]
			value.append( [ ptclIndx, bestAliParams, bestScore] )
			classes.update({ bestreftag : value })
			
	
		#klassIndx = 0
		klassesLen = len(classes)
		#print "(e2spt_refinemulti.py) THERE ARE THESE MANY surviving CLASSESS with particles",len(classes)
		for klass in classes:
			#print "\n\nThe particles and their aliparams, for this class", klass
			#print "are:", classes[ klass ]
		
			klassIndx = int( klass.replace('ref','') )
		
			#ptclsFinal = {}
			#for key in klass:
				#ptclnum = int)
				#ptclTransform = klass[ key ]
				#ptclsFinal.update({ ptclnum : ptclTransform })
		
			if classes[ klass ]:
				print "\n\nWill average. There are particles in the klass %d see" %( klassIndx ) 
				print classes[klass]
				print "\n\n"
				ret = makeAverage( options, classes[klass], klassIndx, klassesLen, it, finalize, originalCompletePath)
				
				if ret:
					avgsName =  originalCompletePath + '/classAvgs.hdf'
					if options.savesteps and int( options.iter ) > 1 :
						avgsName = originalCompletePath + '/classAvgs_iter' + str( it ).zfill( len( str( options.iter ))) + '.hdf'						
					ret.write_image(avgsName,-1)
			else:
				print "The klass %d was empty (no particles were assgined to it). You might have too many classes." % ( klassIndx )	
				#dummyClassAvg=EMData(boxsize,boxsize,boxsize)
				#dummyClassAvg.to_zero()
				ret = None
			
			avgs.update({ klass : ret })
		
			#klassIndx += 1				
			#os.system(cmd)
		
		nonNulls=0
		for avg in avgs:
			if avgs[ avg ] and avgs[ avg ] != None:
				nonNulls +=1
		
		if nonNulls < 2:
			print "e2spt_refinemulti.py has allocated all particles to one average and therefore has failed/converged. EXITING."
			sys.exit()
			
			
		
		print "\n\n\n(e2spt_refinemulti.py) Final Averages are", avgs
		print "\n\n\n"
		
		logger = E2init(sys.argv,options.ppid)	
		E2end(logger)
	
	return()


def genrefs( options, originalCompletePath ):
	
	nptcls = EMUtil.get_image_count( options.input )
	print "\n(e2spt_refinemulti.py) Inside genrefs, originalPath is", originalCompletePath
	
	refsFyles = []
	nrefs = int(options.nrefs)
	
	groupsize = nptcls
	if nrefs > 1:
		groupsize = int( int(nptcls)/int(options.nrefs) )
	
	print "\nTherefore, groupsize is", groupsize
	
	if options.refsgenmethod == 'binarytree' or 'binary' in options.refsgenmethod:
				
		print "\n(e2spt_refinemulti.py) (genrefs) refsgenmethod and nrefs are", options.refsgenmethod, nrefs	

			
		for i in range(nrefs):
			print "\nIterating over nrefs; i is", i
						
			bottom_range = i * groupsize
			top_range = (i+1) * groupsize - 1	#Since e2proc3d.py includes the top range
												#for a set of 16 particles, for example
												#the top range would be index 15, becuase
												#numeration starts at 0. Therefore, you need to
												#subtract 1 to "top_range" if separating the particles
												#using e2proc3d.py
			if i == options.nrefs - 1:
				top_range = nptcls - 1
			print "\nbottom and top ranges are", bottom_range, top_range
			
			#groupPATH = options.input
			if nrefs > 1:
				groupID = 'group' + str(i+1).zfill(len(str(options.nrefs)))
				
				#groupDIR = originalPath + '/' + options.path + '/' + groupID
				
				groupStack = originalCompletePath + '/' + os.path.basename(options.input).replace('.hdf','group' + str(i+1).zfill(len(str(options.nrefs) ) ) + 'ptcls.hdf')
				
				divisioncmd = 'e2proc3d.py ' + options.input + ' ' + groupStack + ' --append --first=' + str(bottom_range) + ' --last=' + str(top_range)
				
			
				alicmd ='cd ' + originalCompletePath + ' && e2spt_classaverage.py --path=' + groupID + ' --input=' + groupStack

				names = dir(options)
				for name in names:
					if getattr(options,name) and 'refs' not in name and 'syms' not in name and 'iter' not in name and 'output' not in name and 'input' not in name and "__" not in name and "_" not in name and 'path' not in name and str(getattr(options,name)) != 'True' and 'iter' not in name:	
						#if "__" not in name and "_" not in name and str(getattr(options,name)) and 'path' not in name and str(getattr(options,name)) != 'False' and str(getattr(options,name)) != 'True' and str(getattr(options,name)) != 'None':			
						alicmd += ' --' + name + '=' + str(getattr(options,name))
				
				alicmd += ' --output=' + groupID + 'avg.hdf --iter=1'
			
				#alicmd += ' && mv ' + groupStack + ' ' + groupID

				finalBinTreeCmd = divisioncmd + ' && ' + alicmd
				
				print "\nfinalBinTreeCmd is", finalBinTreeCmd
				
				p=subprocess.Popen( finalBinTreeCmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				text=p.communicate()	
				p.stdout.close()
			
				refFyle = originalCompletePath + '/' + groupID + '/' + groupID + 'avg.hdf'
				refsFyles.append( refFyle )
		print "\nI have finished building references with the binary tree method"
		
	elif options.refsgenmethod == 'hac':
		avsacmd ='cd ' +  originalCompletePath + ' && e2spt_hac.py --path=sptTMP'
		print "(e2spt_refinemulti.py) (genrefs) refsgenmethod and nrefs are", options.refsgenmethod, nrefs
	
		names = dir(options)
		for name in names:
			if getattr(options,name) and 'keep' not in name and 'syms' not in name and 'refs' not in name and 'iter' not in name and "__" not in name and "_" not in name and 'path' not in name and str(getattr(options,name)) != 'True' and 'iter' not in name:	
				#if "__" not in name and "_" not in name and str(getattr(options,name)) and 'path' not in name and str(getattr(options,name)) != 'False' and str(getattr(options,name)) != 'True' and str(getattr(options,name)) != 'None':			
				avsacmd += ' --' + name + '=' + str(getattr(options,name))
		#avsacmd += ' --autocenter --groups=' + str(nrefs) + ' --iter=' + str(groupsize) + ' && mv ' + originalCompletePath + '/sptTMP/* ' + originalCompletePath + ' && rm -r sptTMP'

		avsacmd += ' --autocenter --groups=' + str(nrefs) + ' --iter=' + str(groupsize)  + ' && mv ' + originalCompletePath + '/sptTMP/* ' + originalCompletePath + ' && rm -r sptTMP'

		print "\n\n\n\n\n(e2spt_refinemulti.py) (genrefs) the command for hac ref generation is", avsacmd
		
		p=subprocess.Popen( avsacmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
		
		for i in range(options.nrefs):
			groupID = 'group' + str(i+1).zfill(len(str(options.nrefs)))
			os.system('mv ' + originalCompletePath + '/' + groupID + '/finalAvg.hdf ' +  originalCompletePath + '/' + groupID +  '/' + groupID + 'avg.hdf')
		
				
			refFyle = originalCompletePath + '/' + groupID + '/' + groupID + 'avg.hdf'
			refsFyles.append( refFyle )
		print "I have finished building references with the hac method"
		
	print "\nDone building references. The returned references are", refsFyles

	return refsFyles
	
	
def makeAverage(options, klass, klassIndx, klassesLen, iterNum, finalize, originalCompletePath):
	"""Will take a set of alignments and an input particle stack filename and produce a new class-average.
	Particles may be excluded based on the keep and keepsig parameters. If keepsig is not set, then keep represents
	an absolute fraction of particles to keep (0-1). Otherwise it represents a sigma multiplier akin to e2classaverage.py"""

	if options.averager: 
		parsedAverager=parsemodopt(options.averager)
		
	print "\n\n\nTHe parsed averager is!!!", parsedAverager
	print "\n"
			
	if options.keepsig:
		# inefficient memory-wise
		val = sum([ score[-1] for score in klass ])
		val2 = sum([ score[-1]**2 for score in klass ])

		mean = val/len( klass )
		sig = sqrt( val2 /len(klass) - mean*mean )
		thresh = mean + sig*options.keep
		if options.verbose: 
			print "Keep threshold : %f (mean=%f  sigma=%f)"%(thresh,mean,sig)

	if options.keep:
		#print "Len of align_parms is", len(klass)
		
		for score in klass:
			if score[-1]:
				pass
				#print "\nscore!"
			else:
				print "\n(e2pt_refinemulti.py) (makeAverage) the score was 0, see", score[-1] 
				#print "see, p[0] is", p[0]
				#sys.exit()
		
		val = [ score[-1] for score in klass]
		val.sort()
		print "The len of val is", len(val)
		print "these are the vals", val
		print "Which shuld be the same as len(klass) see", len(klass)
		print "The clossest position to threshold value based on keep", options.keep
		threshIndx =  int (options.keep * len(klass) ) - 1 
		print "is", threshIndx
		
		thresh = val[ threshIndx ]
		if options.verbose: 
			print "Keep threshold : %f (min=%f  max=%f)"%(thresh,val[0],val[-1])

	'''
	# Make variance image if available
	variance = EMData( ptcl_file, 0 ).copy_head()
	if options.averager[0] == 'mean':
		options.averager[1]['sigma'] = variance
	'''
	
	avgr = Averagers.get(parsedAverager[0], parsedAverager[1])
	included = []
	
	#print "The path to save the class average is", options.path
			
	#jsdict = path + '/tomo_xforms.json'
	#js = js_open_dict(jsdict)
			
	ptclsAdded = 0
	for k in klass:
		print "klass is", klass
		
		if klass and len(klass) > 0:
			#print "\n\nk in klass is", k
			#print "\nThe index of the particle to add is",k[0]
			#print "\nAnd this its transform", k[1][0]
			ptcl = EMData(options.input,k[0])
		
			#print "\n\n\n(e2spt_refinemuti.py) in makeAverage, ptcl and its type are",ptcl,type(ptcl)
		
			ptclTransform =k[1][0]
			#print "And the ptcl transform is", ptclTransform
			ptcl.process_inplace("xform",{"transform" : ptclTransform})
		
			#print "I've applied the transform"
			#print "I have applied this transform before averaging", ptcl_parms[0]["xform.align3d"]			
		
			if k[-1] <= thresh: 
				avgr.add_image(ptcl)
				included.append(k[0])
				ptclsAdded += 1

			#js["tomo_%04d"%i] = ptcl_parms[0]['xform.align3d']
		
			if int(options.iter) -1 == iterNum:
				finalize = 1
		
			if options.saveali and finalize:
				ptcl['origin_x'] = 0
				ptcl['origin_y'] = 0		#The origin needs to be reset to ZERO to avoid display issues in Chimera
				ptcl['origin_z'] = 0
				ptcl['spt_score'] = k[-1]
			
				#print "\nThe score is", ptcl_parms[0]['score']
				#print "Because the zero element is", ptcl_parms[0]
			
				ptcl['xform.align3d'] = Transform()
				#ptcl['spt_ali_param'] = ptcl_parms[0]['xform.align3d']
				ptcl['xform.align3d'] = ptclTransform
			
				print "\n\nFinal iteration, and options.saveali on, so saving class_ptcls\n\n"
			
				classStack = originalCompletePath + "/class" + str( klassIndx ).zfill( len( str (klassesLen))) + "_ptcl.hdf"
				#print "The class name is", classname
				#sys.exit()
				ptcl.write_image(classStack,-1)	
	#js.close()
	
	if options.verbose: 
		print "Kept %d / %d particles in average"%(len(included),len(klass))

	avg=avgr.finish()
	#if options.symmetry and not options.breaksym:
	#	avg=avg.process('xform.applysym',{'sym':options.symmetry})
	
	avg["class_ptcl_idxs"] = included
	avg["class_ptcl_src"] = options.input
	
	#if options.averager[0] == 'mean' and variance:
	#	variance.write_image(path+"/class_varmap.hdf",it)
				
	if not options.nocenterofmass:
		avg.process_inplace("xform.centerofmass")
	
	
	
	if ptclsAdded > 0:
		#avg.write_image(avgsName,klassIndx)
		
		return avg
	else:
		return 0



	
if __name__ == '__main__':
	main()

