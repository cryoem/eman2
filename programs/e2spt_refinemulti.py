#!/usr/bin/env python

#
# Author: Jesus Galaz  July/16/2013
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
	[If you want to refine the data against ONE model, use e2spt_classaverage.py]
	The data set is divded into that same number of groups (4).
	 
	An initial model will be generated with the data in each group.
	Then, the entire data set will be refined against all 4 initial models.
	
	You can increase the number of references used for each iteration by specifying the --addmodel parameter.
	This will take the "best initial model" (the one that most particles prefered) and include it as an initial model for the next round of refinement.
	For exampe, if you start with two references A and B, two averages will come out of aligning the data against them, avgA and avgB.
	So if --addmodel is on, instead of only using avgA and avgB as references for the next refinement round, the best of A and B will also be used,
	which means you will refine the data against 3 models in the next round, not just 2.
	
	If you supply a single reference/model then --addmodel MUST be supplied too; otherwise, to refine a data set against a single model use
	e2spt_classaverage.py
	 """

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	#parser.add_argument("--ncls", type=int, help="...", default=2)
	#parser.add_argument("--nbasis", type=int, help="Basis vectors to use", default=3)

	parser.add_argument("--refs", type=str, help="""This can either be an HDF stack, where each image will be treated as a separate model/reference, 
													or a comma separatted list of individual images; e.g. --refs=ref1.hdf,ref2.hdf,ref3.hdf""", default=None)
	
	'''
	PARAMETERS TO BE PASSED ON TO e2spt_classaverage.py
	'''
	
	parser.add_header(name="caheader", help='Options below this label are specific to e2spt_classaverage', title="### e2spt_classaverage options ###", default=None, guitype='filebox', row=3, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--path",type=str,default='spt_refinemulti',help="Directory to store results in. The default is a numbered series of directories containing the prefix 'spt_refinemulti'; for example, spt_refinemulti02 will be the directory by default if 'spt_refinemulti01' already exists.")
	parser.add_argument("--input", type=str, help="The name of the input volume stack. MUST be HDF or BDB, since volume stack support is required.", default=None, guitype='filebox', browser='EMSubTomosTable(withmodal=True,multiselect=False)', row=0, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--output", type=str, help="The name of the output class-average stack. MUST be HDF or BDB, since volume stack support is required.", default=None, guitype='strbox', row=2, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#parser.add_argument("--oneclass", type=int, help="Create only a single class-average. Specify the class number.",default=None)
	#parser.add_argument("--classmx", type=str, help="The name of the classification matrix specifying how particles in 'input' should be grouped. If omitted, all particles will be averaged.", default='')
	
	#parser.add_argument("--ref", type=str, help="Reference image(s). Used as an initial alignment reference and for final orientation adjustment if present. This is typically the projections that were used for classification.", default=None, guitype='filebox', browser='EMBrowserWidget(withmodal=True,multiselect=True)', filecheck=False, row=1, col=0, rowspan=1, colspan=3, mode='alignment')
	
	#parser.add_argument("--resultmx",type=str,help="Specify an output image to store the result matrix. This is in the same format as the classification matrix. http://blake.bcm.edu/emanwiki/EMAN2/ClassmxFiles", default=None)
	
	#parser.add_argument("--refinemultireftag", type=str, help="DO NOT USE THIS PARAMETER. It is passed on from e2spt_refinemulti.py if needed.", default='')
	
	parser.add_argument("--iter", type=int, help="The number of iterations to perform. Default is 1.", default=1, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, nosharedb=True, mode='alignment,breaksym')
	parser.add_argument("--savesteps",action="store_true", help="If set, will save the average after each iteration to class_#.hdf. Each class in a separate file. Appends to existing files.",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='alignment,breaksym')
	parser.add_argument("--saveali",action="store_true", help="If set, will save the aligned particle volumes in class_ptcl.hdf. Overwrites existing file.",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')
	parser.add_argument("--saveallalign",action="store_true", help="If set, will save the alignment parameters after each iteration",default=False, guitype='boolbox', row=4, col=2, rowspan=1, colspan=1, mode='alignment,breaksym')
	parser.add_argument("--sym", dest = "sym", default=None, help = "Symmetry to impose - choices are: c<n>, d<n>, h<n>, tet, oct, icos", guitype='symbox', row=9, col=1, rowspan=1, colspan=2, mode='alignment,breaksym')
	parser.add_argument("--mask",type=str,help="Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2", returnNone=True, default="mask.sharp:outer_radius=-2", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'mask\')', row=11, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	parser.add_argument("--normproc",type=str,help="Normalization processor applied to particles before alignment. Default is to use normalize.mask. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'", default="normalize.mask")
	
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
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. Default is refine.3d, specify 'None' to disable", default="refine_3d", returnNone=True, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine.*3d\')', row=14, col=0, rowspan=1, colspan=3, nosharedb=True, mode='alignment,breaksym[None]')
	parser.add_argument("--raligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo", guitype='comboparambox',choicelist='re_filter_list(dump_cmps_list(),\'tomo\')', row=15, col=0, rowspan=1, colspan=3,mode="alignment,breaksym")
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
	
	parser.add_argument("--shrink", type=int,default=0,help="Optionally shrink the input volumes by an integer amount for coarse alignment.", guitype='shrinkbox', row=5, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')
	parser.add_argument("--shrinkrefine", type=int,default=0,help="Optionally shrink the input volumes by an integer amount for refine alignment.", guitype='intbox', row=5, col=2, rowspan=1, colspan=1, mode='alignment')
	
	#parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default='', guitype='strbox', row=19, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')

	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1", guitype='strbox', row=19, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#parser.add_argument("--automask",action="store_true",help="Applies a 3-D automask before centering. Can help with negative stain data, and other cases where centering is poor.")
	#parser.add_argument("--resample",action="store_true",help="If set, will perform bootstrap resampling on the particle data for use in making variance maps.",default=False)
	#parser.add_argument("--odd", default=False, help="Used by EMAN2 when running eotests. Includes only odd numbered particles in class averages.", action="store_true")
	#parser.add_argument("--even", default=False, help="Used by EMAN2 when running eotests. Includes only even numbered particles in class averages.", action="store_true")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	'''
	Parameters to compensate for the missing wedge using --cpm=fsc.tomo
	'''
	parser.add_argument("--wedgeangle",type=float,help="""Missing wedge angle, calculated as 90 minus the value yo provide, times 2. 
														For example, --wedgeangle=60 will represent a wedge of size (90-60)*2=60.
														--wedgeangle=70, results in a narrower wedge of size (90-70)*2=40.
														In reality, you should enter here the range of your DATA COLLECTION.
														I.e., if you collected your tiltseries from -60 to 60, enter --wedgeangle=60.""",default=60.0)
	parser.add_argument("--wedgei",type=float,help="Missingwedge begining (in terms of its 'height' along Z. If you specify 0, the wedge will start right at the origin.", default=0.10)
	parser.add_argument("--wedgef",type=float,help="Missingwedge ending (in terms of its 'height' along Z. If you specify 1, the wedge will go all the way to the edge of the box.", default=0.9)
	parser.add_argument("--fitwedgepost", action="store_true", help="Fit the missing wedge AFTER preprocessing the subvolumes, not before, IF using the fsc.tomo comparator for --aligncmp or --raligncmp.", default=False)
	parser.add_argument("--writewedge", action="store_true", help="Write a subvolume with the shape of the fitted missing wedge if --raligncmp or --aligncmp are fsc.tomo. Default is 'True'. To turn on supply --writewedge", default=False)		
	
	(options, args) = parser.parse_args()
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	options = sptmakepath(options,'spt_refinemulti')

	rootpath = os.getcwd()
	
	try:
		hdr = EMData(options.input,0,True) #This is done just to test whether the input file exists where it should
	except:
		print "ERROR: Can't find the file provided through --input"
	
	'''
	Determine how many references there are and separate them if bundled up in one HDF file
	'''
	nrefs=0
	
	refsfiles = []

	if ',' in options.refs:
		refsorig = options.refs.split(',')
		nrefs = len(refsorig)
		for i in range(nrefs):
			outref = options.path + '/' + refsorig[i].replace('.hdf','_ref' + str(i).zfill(len(str(nrefs))) + '.hdf' )
			os.system('cp ' + refsorig[i] + ' ' + outref)		
			#nrefs = len(refsorig)
			refsfiles.append(outref)
			
			
	else:
		nrefs = EMUtil.get_image_count(options.refs)
		for i in range(nrefs):
			outref = options.path + '/' + options.refs.replace('.hdf','_ref' + str(i).zfill(len(str(nrefs))) + '.hdf' )
			os.system('e2proc3d.py ' + options.refs + ' ' + outref + ' --first=' + str(i) + ' --last=' + str(i) )
			refsfiles.append(outref)
	
	'''
	Generate the commands to refine the data against each reference
	'''
	
	#filesindir = os.listdir(rootpath)

	options.input = '../' + options.input 		#The program goes into --path to execute the alignment command; therefore, --input will be one level furtherback
	
	k=0
	reftags = []
	masterInfo = {}
	for ref in refsfiles:
		print "Aligning data to ref number", k
		ref = rootpath + '/' + ref
		
		thisRefinementPath = ref.split('/')[-1].replace('.hdf','')
		
		alicmd = 'cd ' + options.path + ' && e2spt_classaverage.py --ref=' + ref + ' --path=' + thisRefinementPath
	
		names = dir(options)
		for name in names:
			if getattr(options,name) and 'refs' not in name and "__" not in name and "_" not in name and 'path' not in name and str(getattr(options,name)) != 'True':	
				#if "__" not in name and "_" not in name and str(getattr(options,name)) and 'path' not in name and str(getattr(options,name)) != 'False' and str(getattr(options,name)) != 'True' and str(getattr(options,name)) != 'None':			
				alicmd += ' --' + name + '=' + str(getattr(options,name))
		alicmd += ' --donotaverage'
		
		tag = str(k).zfill( len( str ( nrefs )))
		reftag = 'ref' + tag
		reftags.append(reftag)
	
		alicmd += ' --refinemultireftag=' + tag		
		
		#print "Command is", alicmd
		
		p=subprocess.Popen( alicmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		'''
		Option A
		'''
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
		p.communicate()	
		p.stdout.close()
			
		scoresFile = options.path + '/' + thisRefinementPath + '/subtomo_scores' + tag + '.json'
		scores = js_open_dict(scoresFile)
		nscores = len(scores)
		
		print "\n\nThe scores file to read is", scoresFile
		print "I read these many scores", nscores
		print "scores are", scores
		
		#for ele in scores:
		#	print "one score element is", ele
		#	#print "And therefore score is", scores[ele]
			
		print "\n\n"
		
		aliParamsFile = options.path + '/' + thisRefinementPath + '/tomo_xforms' + tag + '.json'
		aliParams = js_open_dict(aliParamsFile)
		nparams = len(aliParams)
		
		#print "The aliParams file to read is", aliParamsFile
		#print "I read these many params", nparams
		
		if nparams != nscores:
			print "nscores is", nscores
			print "nparams is", nparams
			print "WARNING! They should be the same."	
		
		for i in range(nparams):
			ptclID = "tomo_" + str(i).zfill( len(str( nparams )) )
			ptclScore = float( scores[ptclID] )
			ptclAliParams = aliParams[ptclID]
			
			infolist = [ptclScore,ptclAliParams,reftag]
			
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
			
			#print "masterInfo has been updated and now is",masterInfo
		k+=1
		
				
	from operator import itemgetter						
	
	print "I've aligned all the particles in the data set to all the references and will now classify them from the masterInfo dict", masterInfo
	
	classes = {}
	for reftag in reftags:
		classes.update({ reftag : [] })
	
	for ele in masterInfo:
		sortedPtclInfo = sorted( masterInfo[ele], key=itemgetter(0))	#Sorted works because you want the scores from SMALLEST to BIGGEST. Remember, the MORE NEGATIVE (smaller) the better score in EMAN2
		bestPtclInfo = sortedPtclInfo[0]
		bestreftag = bestPtclInfo[-1]	
		bestAliParams = bestPtclInfo[1]
		
		print "\n\nFor particle", ele
		print "The sorted data is", sortedPtclInfo
		print "\n\n"
		
		value = classes[ bestreftag ]
		value.append( [ele,bestAliParams] )
		classes.update({ bestreftag : value })
		
	for klass in classes:
		print "\n\nThe particles and their aliparams, for this class", klass
		print "are:", classes[ klass ]
		
				
		
		#os.system(cmd)
		
		
		
		
		
	logger = E2init(sys.argv,options.ppid)
	

	
	E2end(logger)
	
	return()
	
	
if __name__ == '__main__':
	main()

