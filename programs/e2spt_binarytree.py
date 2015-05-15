#!/usr/bin/env python

#
# Author: Jesus Galaz-Montoya 01/Nov/2014
# Last modification: 19/Feb/2015
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

from EMAN2 import *
import os
import sys
from EMAN2jsondb import JSTask,jsonclasses
from pprint import pprint


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <output> [options]
	Program to build an initial subtomogram average by averaging pairs from the largest subset
	in --input that is a power of 2. For example, if you supply an input stack with 100 subtomograms,
	this program will build an initial reference using 64, since 64 is the largest power of 2 contained in 100.
	In the first iteration, particle 1 will be averaged with 2, 3 with 4, 5 with 6... etc.
	32 new averages (each an average of 2 subtomograms) will be used for the second iteration.
	Again, 1 will be averaged with 2, 3 with 4, etc... yielding 16 new averages.
	The algorithm continues until the entire subset (64) has been merged into 1 average.
	
	This program depends on e2spt_classaverage.py because it imports the preprocessing 
	and alignment functions from it.
	
	--mask=mask.sharp:outer_radius=<safe radius>
	--preprocess=filter.lowpass.gauss:cutoff_freq=<1/resolution in A>
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_header(name="sptbtheader", help="""Options below this label are specific to 
		sptbinarytree""", title="### sptbinarytree options ###", row=6, col=0, rowspan=1, colspan=3,mode="align")
	
	parser.add_header(name="caheader", help="""Options below this label are specific to sptclassaverage""", title="### sptclassaverage options ###", row=3, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--path",type=str,default='spt',help="""Default=spt. Directory to store results in. The default is a numbered series of directories containing the prefix 'spt'; for example, spt_02 will be the directory by default if 'spt_01' already exists.""")
	
	parser.add_argument("--input", type=str, default='',help="""Default=None. The name of the input volume stack. MUST be HDF since volume stack support is required.""", guitype='filebox', browser='EMSubTomosTable(withmodal=True,multiselect=False)', row=0, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--output", type=str, default='avg.hdf', help="""Default=avg.hdf. The name of the output class-average stack. MUST be HDF since volume stack support is required.""", guitype='strbox', row=2, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	#parser.add_argument("--classmx", type=str, default='', help="""Default=None. The name of the classification matrix specifying how particles in 'input' should be grouped. If omitted, all particles will be averaged.""")
	
	#parser.add_argument("--ref", type=str, default='', help="""Default=None. Reference image(s). Used as an initial alignment reference and for final orientation adjustment if present. This is typically the projections that were used for classification.""", guitype='filebox', browser='EMBrowserWidget(withmodal=True,multiselect=True)', filecheck=False, row=1, col=0, rowspan=1, colspan=3, mode='alignment')
	
	#parser.add_argument("--refpreprocess",action="store_true",default=False,help="""Default=False. This will preprocess the reference identically to the particles. It is off by default, but it is internally turned on when no reference is supplied.""")
	
	#parser.add_argument("--resultmx",type=str,default=None,help="""Default=Npone. Specify an output image to store the result matrix. This is in the same format as the classification matrix. http://blake.bcm.edu/emanwiki/EMAN2/ClassmxFiles""")
	
	#parser.add_argument("--refinemultireftag", type=str, default='', help="""Default=''. DO NOT USE THIS PARAMETER. It is passed on from e2spt_refinemulti.py if needed.""")

	parser.add_argument("--radius", type=float, default=0, help="""Default=0 (which means it's not used by default). Hydrodynamic radius of the particle in Angstroms. This will be used to automatically calculate the angular steps to use in search of the best alignment. Make sure the apix is correct on the particles' headers, sine the radius will be converted from Angstroms to pixels. Then, the fine angular step is equal to 360/(2*pi*radius), and the coarse angular step 4 times that.""")
	
	parser.add_argument("--precision",type=float,default=1.0,help="""Default=1.0. Precision in pixels to use when figuring out alignment parameters automatically using --radius. Precision would be the number of pixels that the the edge of the specimen is moved (rotationally) during the finest sampling, --falign. If precision is 1, then the precision of alignment will be that of the sampling (apix of your images) times the --shrinkfine factor specified.""")
	
	parser.add_argument("--search", type=int,default=8,help=""""Default=8. During COARSE alignment translational search in X, Y and Z, in pixels. Default=8. This WILL overwrite any search: provided through --align, EXCEPT if you provide --search=8, which is the default. In general, just avoid providing search twice (through here and through the aligner, --align). If you do, just be careful to make them consistent to minimize misinterpretation and error.""")
	
	parser.add_argument("--searchfine", type=int,default=2,help=""""Default=2. During FINE alignment translational search in X, Y and Z, in pixels. Default=2. This WILL overwrite any search: provided through --falign, EXCEPT if you provide --searchfine=2, which is the default. In general, just avoid providing search twice (through here and through the fine aligner --falign). If you do, just be careful to make them consistent to minimize misinterpretation and error.""")
	
	#parser.add_argument("--donotaverage",action="store_true", help="""If e2spt_refinemulti.py is calling e2spt_classaverage.py, the latter need not average any particles, but rather only yield the alignment results.""", default=False)
	
	parser.add_argument("--iterstop", type=int, default=0, help="""Default=0. (Not used). The program is called to convergence by default (all particles merge into one final average). To stop at an intermediate iteration, provide this parameter. For example, --iterstop=1, will only allow the algorithm to complete 1 iteration; --iterstop=2 will allow it to go through 2, etc.""")
	
	parser.add_argument("--savesteps",action="store_true", default=False, help="""Default=False. If set, will save the average after each iteration to class_#.hdf. Each class in a separate file. Appends to existing files.""", guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	parser.add_argument("--saveali",action="store_true", default=False, help="""Default=False. If set, will save the aligned particle volumes in class_ptcl.hdf. Overwrites existing file.""", guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	parser.add_argument("--saveallalign",action="store_true", default=False, help="""Default=False. If set, will save the alignment parameters after each iteration""", guitype='boolbox', row=4, col=2, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	parser.add_argument("--sym", dest = "sym", default='', help = """Default=None (equivalent to c1). Symmetry to impose -choices are: c<n>, d<n>, h<n>, tet, oct, icos""", guitype='symbox', row=9, col=1, rowspan=1, colspan=2, mode='alignment,breaksym')
	
	parser.add_argument("--mask",type=str,default="mask.sharp:outer_radius=-2", help="""Default is mask.sharp:outer_radius=-2. Masking processor applied to particles before alignment. IF using --clipali, make sure to express outer mask radii as negative pixels from the edge.""", returnNone=True, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'mask\')', row=11, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--maskfile",type=str,default='',help="""Default=None. Mask file (3D IMAGE) applied to particles before alignment. Must be in HDF format. Default is None.""")
	
	parser.add_argument("--normproc",type=str, default='normalize.edgemean',help="""Default is 'normalize.edgemean' (see 'e2help.py processors -v 10' at the command line). Normalization processor applied to particles before alignment. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'""")
	
	parser.add_argument("--threshold",type=str,default='',help="""Default=None. A threshold applied to the subvolumes after normalization. For example, --threshold=threshold.belowtozero:minval=0 makes all negative pixels equal 0, so that they do not contribute to the correlation score.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--preprocess",type=str,default='',help="""Any processor (see 'e2help.py processors -v 10' at the command line) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--preprocessfine",type=str,default='',help="""Any processor (see 'e2help.py processors -v 10' at the command line) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.""")
	
	parser.add_argument("--lowpass",type=str,default='',help="""Default=None. A lowpass filtering processor (see 'e2help.py processors -v 10' at the command line) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=17, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--lowpassfine",type=str,default='',help="""Default=None. A lowpass filtering processor (see 'e2help.py processors -v 10' at the command line) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.""")

	parser.add_argument("--highpass",type=str,default='',help="""Default=None. A highpass filtering processor (see 'e2help.py processors -v 10' at the command line) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.""", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=18, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--highpassfine",type=str,default='',help="""Default=None. A highpass filtering processor (see 'e2help.py processors -v 10' at the command line) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.""")

	parser.add_argument("--shrink", type=int,default=1,help="""Default=1 (no shrinking). Optionally shrink the input volumes by an integer amount for coarse alignment.""", guitype='shrinkbox', row=5, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	parser.add_argument("--shrinkfine", type=int,default=1,help="""Default=1 (no shrinking). Optionally shrink the input volumes by an integer amount for refine alignment.""", guitype='intbox', row=5, col=2, rowspan=1, colspan=1, mode='alignment')
	
	parser.add_argument("--clipali",type=int,default=0,help="""Default=0 (which means it's not used). Boxsize to clip particles as part of preprocessing to speed up alignment. For example, the boxsize of the particles might be 100 pixels, but the particles are only 50 pixels in diameter. Aliasing effects are not always as deleterious for all specimens, and sometimes 2x padding isn't necessary; still, there are some benefits from 'oversampling' the data during averaging; so you might still want an average of size 2x, but perhaps particles in a box of 1.5x are sufficiently good for alignment. In this case, you would supply --clipali=75""")
	
	parser.add_argument("--postprocess",type=str,default='',help="""A processor to be applied to the FINAL volume after averaging the raw volumes in their FINAL orientations, after all iterations are done.""",guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=16, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--procfinelikecoarse",action='store_true',default=False,help="""If you supply this parameters, particles for fine alignment will be preprocessed identically to particles for coarse alignment by default. If you supply this, but want specific parameters for preprocessing particles for also supply: fine alignment, nd supply fine alignment parameters, such as --lowpassfine, --highpassfine, etc; to preprocess the particles for FINE alignment differently than for COARSE alignment.""")
	
	parser.add_argument("--npeakstorefine", type=int, help="""Default=1. The number of best coarse alignments to refine in search of the best final alignment. Default=1.""", default=4, guitype='intbox', row=9, col=0, rowspan=1, colspan=1, nosharedb=True, mode='alignment,breaksym[1]')
	
	parser.add_argument("--align",type=str,default="rotate_translate_3d:search=8:delta=12:dphi=12",help="""This is the aligner used to align particles to the previous class average. Default is rotate_translate_3d:search=8:delta=12:dphi=12, specify 'None' (with capital N) to disable.""", returnNone=True,guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'3d\')', row=12, col=0, rowspan=1, colspan=3, nosharedb=True, mode="alignment,breaksym['rotate_symmetry_3d']")
	
	parser.add_argument("--aligncmp",type=str,default="ccc.tomo",help="""Default=ccc.tomo. The comparator used for the --align aligner. Do not specify unless you need to use anotherspecific aligner.""",guitype='comboparambox',choicelist='re_filter_list(dump_cmps_list(),\'tomo\')', row=13, col=0, rowspan=1, colspan=3,mode="alignment,breaksym")
	
	parser.add_argument("--falign",type=str,default="refine_3d_grid:delta=3:range=15:search=2",help="""Default="refine_3d_grid:delta=3:range=15:search=2". This is the second stage aligner used to fine-tune the first alignment. Specify 'None' to disable.""", returnNone=True, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine.*3d\')', row=14, col=0, rowspan=1, colspan=3, nosharedb=True, mode='alignment,breaksym[None]')
		
	parser.add_argument("--faligncmp",type=str,default="ccc.tomo",help="""Default=ccc.tomo. The comparator used by the second stage aligner.""", guitype='comboparambox', choicelist='re_filter_list(dump_cmps_list(),\'tomo\')', row=15, col=0, rowspan=1, colspan=3,mode="alignment,breaksym")		
		
	parser.add_argument("--averager",type=str,default="mean.tomo",help="""Default=mean.tomo. The type of averager used to produce the class average. Default=mean.tomo.""")
		
	#parser.add_argument("--keep",type=float,default=1.0,help="""Default=1.0 (all particles kept). The fraction of particles to keep in each class.""", guitype='floatbox', row=6, col=0, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	#parser.add_argument("--keepsig", action="store_true", default=False,help="""Default=False. Causes the keep argument to be interpreted in standard deviations.""", guitype='boolbox', row=6, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')

	#parser.add_argument("--inixforms",type=str,default="",help="""Default=None. .json file containing a dict of transforms to apply to 'pre-align' the particles.""", guitype='dirbox', dirbasename='spt_|sptsym_', row=7, col=0,rowspan=1, colspan=2, nosharedb=True, mode='breaksym')
	
	parser.add_argument("--breaksym",action="store_true", default=False,help="""Default=False. Break symmetry. Do not apply symmetrization after averaging, even if searching the asymmetric unit provided through --sym only for alignment. Default=False""", guitype='boolbox', row=7, col=2, rowspan=1, colspan=1, nosharedb=True, mode=',breaksym[True]')
	
	#parser.add_argument("--groups",type=int,default=0,help="""Default=0 (not used; data not split). This parameter will split the data into a user defined number of groups. For purposes of gold-standard FSC computation later, select --group=2.""")
		
	parser.add_argument("--randomizewedge",action="store_true",  default=False,help="""Default=False. This parameter is EXPERIMENTAL. It randomizes the position of the particles BEFORE alignment, to minimize missing wedge bias and artifacts during symmetric alignment where only a fraction of space is scanned""")
	
	parser.add_argument("--savepreprocessed",action="store_true",  default=False,help="""Default=False. Will save stacks of preprocessed particles (one for coarse alignment and one for fine alignment if preprocessing options are different).""")
	
	parser.add_argument("--autocenter",type=str, default='',help="""Default=None. Autocenters each averaged pair during initial average generation with --btref and --hacref. Will also autocenter the average of all particles after each iteration of iterative refinement. Options are --autocenter=xform.centerofmass (self descriptive), or --autocenter=xform.centeracf, which applies auto-convolution on the average.""")
	
	parser.add_argument("--autocentermask",type=str, default='',help="""Default=None. Masking processor to apply before autocentering. See 'e2help.py processors -v 10' at the command line.""")
	
	parser.add_argument("--autocenterpreprocess",action='store_true', default=False,help="""Default=False. This will apply a highpass filter at a frequency of half the box size times the apix, shrink by 2, and apply a low pass filter at half nyquist frequency to any computed average for autocentering purposes if --autocenter is provided. Default=False.""")
	
	parser.add_argument("--parallel",default="thread:1",help="""default=thread:1. Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel""", guitype='strbox', row=19, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--ppid", type=int, help="""Default=-1. Set the PID of the parent process, used for cross platform PPID""",default=-1)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="""Default=0. Verbose level [0-9], higner number means higher level of verboseness""")
		
	#parser.add_argument("--resume",type=str,default='',help="""(Not working currently). tomo_fxorms.json file that contains alignment information for the particles in the set. If the information is incomplete (i.e., there are less elements in the file than particles in the stack), on the first iteration the program will complete the file by working ONLY on particle indexes that are missing. For subsequent iterations, all the particles will be used.""")
															
	parser.add_argument("--plots", action='store_true', default=False,help="""Default=False. Turn this option on to generatea plot of the ccc scores during each iteration. Running on a cluster or via ssh remotely might not support plotting.""")

	parser.add_argument("--subset",type=int,default=0,help="""Default=0 (not used). Refine only this substet of particles from the stack provided through --input""")
	
	parser.add_argument("--notmatchimgs",action='store_true',default=False,help="""Default=True. This option prevents applying filter.match.to to one image so that it matches the other's spectral profile during preprocessing for alignment purposes.""")

	parser.add_argument("--preavgproc1",type=str,default='',help="""Default=None. A processor (see 'e2help.py processors -v 10' at the command line) to be applied to the raw particle after alignment but before averaging (for example, a threshold to exclude extreme values, or a highphass filter if you have phaseplate data.)""")
	
	parser.add_argument("--preavgproc2",type=str,default='',help="""Default=None. A processor (see 'e2help.py processors -v 10' at the command line) to be applied to the raw particle after alignment but before averaging (for example, a threshold to exclude extreme values, or a highphass filter if you have phaseplate data.)""")

	parser.add_argument("--weighbytiltaxis",type=str,default='',help="""Default=None. A,B, where A is an integer number and B a decimal. A represents the location of the tilt axis in the tomogram in pixels (eg.g, for a 4096x4096xZ tomogram, this value should be 2048), and B is the weight of the particles furthest from the tomogram. For example, --weighbytiltaxis=2048,0.5 means that praticles at the tilt axis (with an x coordinate of 2048) will have a weight of 1.0 during averaging, while the distance in the x coordinates of particles not-on the tilt axis will be used to weigh their contribution to the average, with particles at the edge(0+radius or 4096-radius) weighing 0.5, as specified by the value provided for B.""")
	
	parser.add_argument("--weighbyscore",action='store_true',default=False,help="""Default=False. This option will weigh the contribution of each subtomogram to the average by score/bestscore.""")


	'''
	BT SPECIFIC PARAMETERS
	'''
	
		
	parser.add_argument("--nseedlimit",type=int,default=0,help="""Maximum number of particles
		to use. For example, if you supply a stack with 150 subtomograms, the program will
		automatically select 128 as the limit to use because it's the largest power of 2 that is
		smaller than 150. But if you provide, say --nseedlimit=100, then the number of particles
		used will be 64, because it's the largest power of 2 that is still smaller than 100.""")
	
	

	(options, args) = parser.parse_args()
	
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'spt_bt')
	
	rootpath = os.getcwd()
	if rootpath not in options.path:
		options.path = rootpath + '/' + options.path
	
	
	if not options.input:
		parser.print_help()
		exit(0)
	elif options.subset:
		subsetStack = options.path + '/subset' + str( options.subset ).zfill( len( str( options.subset))) + '.hdf' 
		print "\nSubset to be written to", subsetStack
		
		subsetcmd = 'e2proc3d.py ' + options.input + ' ' + subsetStack + ' --first=0 --last=' + str(options.subset-1) 
		print "Subset cmd is", subsetcmd
		
		p=subprocess.Popen( subsetcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE )
		text=p.communicate()	
		p.stdout.close()
		
		options.input = subsetStack
		
	from e2spt_classaverage import sptParseAligner
	options = sptParseAligner( options )

	'''
	If --radius of the particle is provided, we calculate the optimal alignment steps for 
	coarse and fine alignment rounds using --shrink and --shrinkfine options and apix info
	'''
	
	if options.shrink < options.shrinkfine:
		options.shrink = options.shrinkfine
		print "It makes no sense for shrinkfine to be larger than shrink; therefore, shrink will be made to match shrinkfine"
	
	if options.radius:
		from e2spt_classaverage import calcAliStep
		options = calcAliStep(options)
	
	'''
	Parse parameters such that "None" or "none" are adequately interpreted to turn of an option
	'''
	
	from e2spt_classaverage import sptOptionsParser
	options = sptOptionsParser( options )
	
	from e2spt_classaverage import writeParameters
	writeParameters(options,'e2spt_binarytree.py', 'bt')
	
					
	hdr = EMData(options.input,0,True)
	nx = hdr["nx"]
	ny = hdr["ny"]
	nz = hdr["nz"]
	if nx!=ny or ny!=nz :
		print "ERROR, input volumes are not cubes"
		sys.exit(1)
		
	logger = E2init(sys.argv, options.ppid)
	
	
	'''
	Initialize parallelism if being used
	'''
	
	if options.parallel :
	
		if options.parallel == 'none' or options.parallel == 'None' or options.parallel == 'NONE':
			options.parallel = ''
			etc = ''
		
		else:
			print "\n\n(e2spt_classaverage.py) INITIALIZING PARALLELISM!"
			print "\n\n"

			from EMAN2PAR import EMTaskCustomer
			etc=EMTaskCustomer(options.parallel)

			pclist=[options.input]

			etc.precache(pclist)
		
	else:
		etc=''
	
	nptcl=EMUtil.get_image_count(options.input)
	if nptcl < 1: 
		print "ERROR : at least 2 particles required in input stack"
		sys.exit(1)
	
	ptclnums=range(nptcl)
	nptclForRef = len(ptclnums)
	
	nseed=2**int(floor(log(len(ptclnums),2)))	# we stick with powers of 2 for this to make the tree easier to collapse
	
	if options.nseedlimit:
		nseed=2**int(floor(log( options.nseedlimit , 2)))
		
	binaryTreeRef(options,nptclForRef,nseed,-1,etc)

		
	print "Will end logger"	
	E2end(logger)
	
	print "logger ended"
	sys.stdout.flush()
	
	return
	

def binaryTreeRef(options,nptclForRef,nseed,ic,etc):
	
	from e2spt_classaverage import Align3DTask,align3Dfunc,get_results
	
	if nptclForRef == 1: 
		print "Error: More than 1 particle required to build a reference."
		sys.exit(1)
			
	# we need to make an initial reference. Due to the parallelism scheme we're using in 3-D and the slow speed of the
	# individual alignments we use a slightly different strategy than in 2-D. We make a binary tree from the first 2^n particles and
	# compute pairwise alignments until we get an average out. 

	
	
	#if nseed>64 : 
	#	nseed=64
	#	print "Limiting seeding to the first 64 images"

	nseediter=int(log(nseed,2))			# number of iterations we'll need
	
	if options.iterstop:
		nseediter=options.iterstop
	
	if options.verbose: 
		print "Seedtree to produce initial reference. Using %d particles in a %d level tree"%(nseed,nseediter)
	
	# We copy the particles for this class into bdb:seedtree_0
	'''
	for i,j in enumerate(ptclnums[:nseed]):
		emdata = EMData(options.input,j)
		
		#if options.inixforms:
		#	emdata.process_inplace("xform",{"transform":js["tomo_%04d"%i]})
		#	emdata.set_attr("test_xfm",js["tomo_%04d"%i])
		
		seedfile = "%s/seedtree_0_cl_%d.hdf" % (options.path,ic)
		emdata.write_image(seedfile, i)
		
		print "Creating this seed file for this class", seedfile, ic
	'''
	
	#for i in range( ptclnums[:nseed] ):
	
	ii=0
	
	seedfile = options.path + '/seedtree_0_cl_' + str(ic) + '.hdf'
	
	if ic < 0:
		seedfile = options.path + '/seedtree_0.hdf'	
	
	#for j in ptclnums[:nseed]:
	for j in range(nseed):
		emdata = EMData(options.input,j)
		emdata.write_image(seedfile,ii)
		print "have taken particle %d and written it into index %d of the seedfile" %(j,ii)
		ii+=1
		if ic >= 0:
			print "Creating this seed file for this class", seedfile, ic
	
	
	'''
	#Outer loop covering levels in the converging binary tree
	'''
	for i in range( nseediter ):
		infile="%s/seedtree_%d_cl_%d.hdf"%(options.path,i,ic)
		if ic < 0:
			infile="%s/seedtree_%d.hdf"%(options.path,i)
			
		print "Infile will be", infile
		
		outfile="%s/seedtree_%d_cl_%d.hdf"%(options.path,i+1,ic)
		if ic < 0:
			outfile="%s/seedtree_%d.hdf"%(options.path,i+1)
	
		if i == nseediter-1:
			outfile = options.path + '/final_avg.hdf'
		
		print "Outfile will be", outfile
	
		tasks=[]
		results=[]
		transform = None
		# loop over volumes in the current level
		
		for j in range(0,nseed/(2**i),2):

			#Unfortunately this tree structure limits the parallelism to the number of pairs at the current level :^(
			if options.parallel:
				#task=Align3DTask(["cache",infile,j],["cache",infile,j+1],j/2,"Seed Tree pair %d at level %d"%(j/2,i),options.mask,options.normproc,options.preprocess,options.lowpass,options.highpass,
				#	options.npeakstorefine,options.align,options.aligncmp,options.falign,options.faligncmp,options.shrink,options.shrinkfine,transform,options.verbose-1,options.randomizewedge,options.wedgeangle,options.wedgei,options.wedgef)
				
				task=Align3DTask(["cache",infile,j],["cache",infile,j+1],j/2,"Seed Tree pair #%d at level %d"%(j/2,i),options,transform,0)
				tasks.append(task)
			else:
				#print "No parallelism specified"
				result=align3Dfunc(["cache",infile,j],["cache",infile,j+1],j/2,"Seed Tree pair #%d at level %d"%(j/2,i),options,transform,0)
				results.append(result['final'])
		'''		
		#Start the alignments for this level
		'''
		if options.parallel:
			tids=etc.send_tasks(tasks)
			if options.verbose: 
				print "%d tasks queued in seedtree level %d"%(len(tids),i) 

			"""Wait for alignments to finish and get results"""
			results=get_results(etc,tids,options.verbose,nseed,'binarytree')

			#results=get_results(etc,tids,options.verbose,{},len(ptclnums),0,'binarytree')
			#results=get_results(etc,tids,options.verbose,{},nptclForRef,0)
			#def get_results(etc,tids,verbose,nptcls,refmethod=''):



			if options.verbose>2 : 
				print "Results:"
				pprint(results)
		else:
			#print "No parallelism specified"
			#results=tasks
			if options.verbose>2 : 
				print "Results:" 
				pprint(results)
						
		makeAveragePairs(options,infile,outfile,results)
		
	#ref = EMData( outfile, 0 )		# result of the last iteration
	
	#return ref
	return
	

def makeAveragePairs(options,ptcl_file,outfile, results):
	"""Will take a set of alignments and an input particle stack filename and produce a new set of class-averages over pairs"""
	
	current = os.getcwd()
	print "\n(e2spt_classaverage.py) (make_average_pairs) current directory is", current
	findir = os.listdir(current)
	print "\noptions.path is", options.path
	findirpath = os.listdir(options.path)
	print "\nThe particle file where the particles ought to be read from is", ptcl_file
	print "\nLets see if ptcl_file is in path. Files in path are", findirpath
	
	print "\nresults are", results
	print "\nTheir len", len(results)
	
	#for i,ptcl_parms in enumerate(align_parms):
	ii=0
	for r in results:
		print "r is", r
		print "\nr[0] is", r[0]
		print "\nr[0][0] is", r[0][0]
		print "\nr[0][0]['xform.align3d'] is", r[0][0]["xform.align3d"]
		print "\nr[0][-1]", r[0][-1]
		
				
		#if ptcl_parms:
		ptcl0=EMData(ptcl_file,ii*2)
		
		ptcl1=EMData(ptcl_file,ii*2+1)
		ptcl1.process_inplace("xform",{"transform":r[0][0]["xform.align3d"]})
	
		#ptcl1.process_inplace("xform",{"transform":align_parms[0]["xform.align3d"]})
	
		# While this is only 2 images, we still use the averager in case something clever is going on
		print "averager is", options.averager
		avgr = Averagers.get(options.averager[0], options.averager[1])
		avgr.add_image(ptcl0)
		avgr.add_image(ptcl1)
	
		avg=avgr.finish()
		#postprocess(avg,optmask,optnormproc,optpostprocess)		#There should be NO postprocessing of the intermediate averages
	
		if options.autocenter:
			print "\n\n\n\nYou have selected to autocenter!\n", options.autocenter
			
			avgac = avg.copy()
			if options.autocentermask:
				avgac.process_inplace( options.autocentermask[0],options.autocentermask[1] )
				
			if options.autocenterpreprocess:
				apix = avg['apix_x']
				halfnyquist = apix*4
				highpassf = apix*a['nx']/2.0
				
				avgac.process_inplace( 'filter.highpass.gauss',{'cutoff_freq':highpassf,'apix':apix})
				avgac.process_inplace( 'filter.lowpass.gauss',{'cutoff_freq':halfnyquist,'apix':apix})
				avgac.process_inplace( 'math.meanshrink',{'n':2})
				
			avgac.process_inplace(options.autocenter[0],options.autocenter[1])
			
			tcenter = avgac['xform.align3d']
			print "Thus the average HAS BEEN be translated like this", tcenter
	
		avg['origin_x']=0
		avg['origin_y']=0
		avg['origin_z']=0
	
		avg.write_image(outfile,ii)
		
		ii+=1
		
	return


if __name__ == '__main__':
	main()

