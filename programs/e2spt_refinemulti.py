#!/usr/bin/env python

#
# Author: Jesus Galaz-Montoya  November/21/2013
# Last modification: 19/Feb/2015
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
from EMAN2jsondb import JSTask,jsonclasses

from e2spt_classaverage import Align3DTask, align3Dfunc, get_results

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
	
	parser.add_header(name="sptrefinemultiheader", help="""Options below this label are specific to sptrefinemulti.""", title="### sptrefinemulti options ###", row=5, col=0, rowspan=1, colspan=3,mode="align")
	
	#parser.add_argument("--ncls", type=int, help="...", default=2)
	#parser.add_argument("--nbasis", type=int, help="Basis vectors to use", default=3)

	parser.add_argument("--input", type=str, help="""The name of the input volume stack. MUST be HDF since volume stack support is required.""", default=None, guitype='filebox', browser='EMSubTomosTable(withmodal=True,multiselect=False)', row=0, col=0, rowspan=1, colspan=3, mode="align")
	
	parser.add_argument("--ref", type=str, help="""Comma separated list of individual images; e.g. --refs=ref1.hdf,ref2.hdf,ref3.hdf. If a single image is provided, several copies will be made based on the number of references specified through --nref.""", default='')
	
	parser.add_argument("--nref", type=int,  default=2, help="""Default=2. (For single reference refinement use e2spt_classaverage.py). Number of references to generate from a single image provided through --ref (random-phase filtered differently), or number of different initial references to generate from scratch from the data set (--input). Default=2""")
	
	parser.add_argument("--refgenmethod", type=str, help="""Method for generating the initial reference(s). Options are 'bt', for binary tree (see e2spt_binarytree.py), 'hac', for hierarchical ascendant classification (see e2spt_hac.py), or 'ssa' for self-symmetry alignment (see e2symsearch3d.py). Default=bt""", default='bt') 
	
	parser.add_argument("--subset4ref",type=int, help=""" Size of the subset of particles to use for generating each reference. Default=0, which means all particles in each subgroup will be used (for example, if --input has 100 particles and --nref is 10, 10 references will be generated using 10 particles for each). If --ref not provided, the program generates an --nref number of references from --input.""")
	
	parser.add_argument("--subset",type=int,default=0,help=""" WARNING: NOT IMPLEMENTED HERE YET. Default=0 (not used). Refine only this substet of particles from the stack provided through --input""")

	
	#parser.add_argument("--refpreprocess",action="store_true",default=False,help="""This 
	#	will preprocess the reference identically to the particles. It is off by default, but it is internally turned on when no reference is supplied.""")
	
	'''
	PARAMETERS TO BE PASSED ON TO e2spt_classaverage.py
	'''
	
	parser.add_argument("--apix",type=float,default=0.0,help="""Default=0.0 (not used). Use this apix value where relevant instead of whatever is in the header of the reference and the particles.""")
	
	parser.add_header(name="caheader", help="""Options below this label are specific to e2spt_classaverage""", title="### e2spt_classaverage options ###", default=None, guitype='filebox', row=3, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--path",type=str,default='',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'spt_refinemulti'; for example, spt_refinemulti02 will be the directory by default if 'spt_refinemulti01' already exists.""")
	
	parser.add_argument("--syms", type=str, help="""List comma-separated symmetries to apply separately on the different references. For example, if you provide --syms=d8,d7 and provide 2 references via --nref=2 or supply two references via --refs=r1.hdf,r2.hdf, d8 symmetry will be applied to the first reference and d7 to the second after each iteration of refinement (the final average in one iteration becomes a reference for the next).""", default='')
	
	parser.add_argument("--output", type=str, help="The name of the output class-average stack. MUST be in  .hdf format, since volume stack support is required.", default=None, guitype='strbox', row=2, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#parser.add_argument("--oneclass", type=int, help="Create only a single class-average. Specify the class number.",default=None)
	#parser.add_argument("--classmx", type=str, help="The name of the classification matrix specifying how particles in 'input' should be grouped. If omitted, all particles will be averaged.", default='')
	
	parser.add_argument("--classmx", type=str, default='', help="""Default=None. WARNING: Not implemented yet. The name of the classification matrix specifying how particles in 'input' should be grouped to generate initial averages and seed refinement.""")

	
	#parser.add_argument("--ref", type=str, help="Reference image(s). Used as an initial alignment reference and for final orientation adjustment if present. This is typically the projections that were used for classification.", default=None, guitype='filebox', browser='EMBrowserWidget(withmodal=True,multiselect=True)', filecheck=False, row=1, col=0, rowspan=1, colspan=3, mode='alignment')
	
	#parser.add_argument("--resultmx",type=str,help="Specify an output image to store the result matrix. This is in the same format as the classification matrix. http://blake.bcm.edu/emanwiki/EMAN2/ClassmxFiles", default=None)
	
	#parser.add_argument("--refinemultireftag", type=str, help="DO NOT USE THIS PARAMETER. It is passed on from e2spt_refinemulti.py if needed.", default='')
	
	parser.add_argument("--refpreprocess",action="store_true",default=False,help="""Default=False. This will preprocess the reference identically to the particles. It is off by default, but it is internally turned on when no reference is supplied. It should probably be off when using a crystal structure (with all positive densities) turned to EM density as an initial model, but it should be on when using an EM map.""")
	
	parser.add_argument("--refrandphase", type=float, default=0, help="""Default=0. Resolution to phase-randomize the reference to.""")
	
	parser.add_argument("--hacref",type=int,default=0,help="""Default=0 (not used by default). Size of the SUBSET of particles to use to build an initial reference by calling e2spt_hac.py which does Hierarchical Ascendant Classification (HAC) or 'all vs all' alignments.""") 
		
	parser.add_argument("--ssaref",type=int,default=0,help="""Default=0 (not used by default). Size of the SUBSET of particles to use to build an initial reference by calling e2symsearch3d.py, which does self-symmetry alignments. You must provide --sym different than c1 for this to make any sense.""")
		
	parser.add_argument("--btref",type=int,default=0,help="""Default=0 (internally turned on and set to 64). Size of the SUBSET of particles to use to build an initial reference by calling e2spt_binarytree.py. By default, the largest power of two smaller than the number of particles in --input will be used. For example, if you supply a stack with 150 subtomograms, the program will automatically select 128 as the limit to use because it's the largest power of 2 that is smaller than 150. But if you provide, say --btref=100, then the number of particles used will be 64, because it's the largest power of 2 that is still smaller than 100.""")
	
	parser.add_argument("--radius", type=float, help="""Will make --align and --falign None. Hydrodynamic radius of the particle in Angstroms. This will be used to automatically calculate the angular steps to use in search of the best alignment. Make sure the apix is correct on the particles' headers, sine the radius will be converted from Angstroms to pixels. Then, the fine angular step is equal to 360/(2*pi*radius), and the coarse angular step 4 times that""", default=0)
	
	parser.add_argument("--iter", type=int, help="""The number of iterations to perform. Default is 1.""", default=1, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, nosharedb=True, mode='alignment,breaksym')
	
	parser.add_argument("--savesteps",action="store_true", help="""If set, will save the average after each iteration to class_#.hdf. Each class in a separate file. Appends to existing files.""",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	parser.add_argument("--saveali",action="store_true", help="""If set, will save the aligned particle volumes in class_ptcl.hdf. Overwrites existing file.""",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	parser.add_argument("--saveallalign",action="store_true", help="""If set, will save the alignment parameters after each iteration""",default=False, guitype='boolbox', row=4, col=2, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	parser.add_argument("--sym", dest = "sym", default=None, help="""Symmetry to impose -choices are: c<n>, d<n>, h<n>, tet, oct, icos""", guitype='symbox', row=9, col=1, rowspan=1, colspan=2, mode='alignment,breaksym')
	
	parser.add_argument("--mask",type=str,help="""Mask processor applied to particles before alignment. Default is mask.sharp:outer_radius=-2""", returnNone=True, default="mask.sharp:outer_radius=-2", guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'mask\')', row=11, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--normproc",type=str,help="""Normalization processor applied to particles before alignment. Default is to use normalize. If normalize.mask is used, results of the mask option will be passed in automatically. If you want to turn this option off specify \'None\'""", default="normalize")
	
	parser.add_argument("--threshold",type=str,help="""A threshold applied to the subvolumes after normalization. For example, --threshold=threshold.belowtozero:minval=0 makes all negative pixels equal 0, so that they do not contribute to the correlation score.""", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--preprocess",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=10, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')

	parser.add_argument("--preprocessfine",type=str,help="Any processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.", default=None)
	
	parser.add_argument("--lowpass",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=17, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')

	parser.add_argument("--lowpassfine",type=str,help="A lowpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--highpass",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to COARSE alignment. Not applied to aligned particles before averaging.", default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=18, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')

	parser.add_argument("--highpassfine",type=str,help="A highpass filtering processor (as in e2proc3d.py) to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.", default=None)

	parser.add_argument("--postprocess",type=str,help="A processor to be applied to the FINAL volume after averaging the raw volumes in their FINAL orientations, after all iterations are done.",default=None, guitype='comboparambox', choicelist='re_filter_list(dump_processors_list(),\'filter\')', row=16, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	
	parser.add_argument("--procfinelikecoarse",type=bool,default=True,help="Turn on with --procfinelikecoarse=False, and supply fine alignment parameters, such as --lowpassfine, --highpassfine, etc; to preprocess the particles for FINE alignment differently than for COARSE alignment.")
	
#	parser.add_argument("--ncoarse", type=int, help="Deprecated. Use --npeakstorefine instead.", default=None)

	parser.add_argument("--npeakstorefine", type=int, default=1, help="Default=1. The number of best coarse alignments to refine in search of the best final alignment.",guitype='intbox', row=9, col=0, rowspan=1, colspan=1, nosharedb=True, mode='alignment,breaksym[1]')

	parser.add_argument("--align",type=str,default="rotate_translate_3d:search=8:delta=12:dphi=12",help="""This is the aligner used to align particles to the previous class average. Default is rotate_translate_3d:search=8:delta=12:dphi=12, specify 'None' (with capital N) to disable.""", returnNone=True,guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'3d\')', row=12, col=0, rowspan=1, colspan=3, nosharedb=True, mode="alignment,breaksym['rotate_symmetry_3d']")

	parser.add_argument("--aligncmp",type=str,help="The comparator used for the --align aligner. Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.",default="ccc.tomo", guitype='comboparambox',choicelist='re_filter_list(dump_cmps_list(),\'tomo\')', row=13, col=0, rowspan=1, colspan=3,mode="alignment,breaksym")

	parser.add_argument("--falign",type=str,default="refine_3d_grid:delta=3:range=15:search=2",help="""Default="refine_3d_grid:delta=3:range=15:search=2". This is the second stage aligner used to fine-tune the first alignment. Specify 'None' to disable.""", returnNone=True, guitype='comboparambox', choicelist='re_filter_list(dump_aligners_list(),\'refine.*3d\')', row=14, col=0, rowspan=1, colspan=3, nosharedb=True, mode='alignment,breaksym[None]')

	parser.add_argument("--faligncmp",type=str,help="The comparator used by the second stage aligner. Default is the internal tomographic ccc",default="ccc.tomo", guitype='comboparambox',choicelist='re_filter_list(dump_cmps_list(),\'tomo\')', row=15, col=0, rowspan=1, colspan=3,mode="alignment,breaksym")

	parser.add_argument("--averager",type=str,help="The type of averager used to produce the class average. Default=mean",default="mean.tomo")

	parser.add_argument("--keep",type=float,help="The fraction of particles to keep in each class.",default=1.0, guitype='floatbox', row=6, col=0, rowspan=1, colspan=1, mode='alignment,breaksym')

	parser.add_argument("--inixforms",type=str,help="directory containing a dict of transform to apply before reference generation", default="", guitype='dirbox', dirbasename='spt_|sptsym_', row=7, col=0,rowspan=1, colspan=2, nosharedb=True, mode='breaksym')

	parser.add_argument("--breaksym",action="store_true", help="Break symmetry. Do not apply symmetrization after averaging", default=False, guitype='boolbox', row=7, col=2, rowspan=1, colspan=1, nosharedb=True, mode=',breaksym[True]')
	
#	parser.add_argument("--groups",type=int,help="WARNING: This parameter is EXPERIMENTAL, and will only work if --iter=1. It's the number of final averages you want from the set after ONE iteration of alignment. Particles will be separated in groups based on their correlation to the reference",default=0)

	parser.add_argument("--randomizewedge",action="store_true", help="This parameter is EXPERIMENTAL. It randomizes the position of the particles BEFORE alignment, to minimize missing wedge bias and artifacts during symmetric alignment where only a fraction of space is scanned", default=False)
	
	parser.add_argument("--savepreprocessed",action="store_true", help="Will save stacks of preprocessed particles (one for coarse alignment and one for fine alignment if preprocessing options are different).", default=False)
	
	parser.add_argument("--keepsig", action="store_true", help="Causes the keep argument to be interpreted in standard deviations.",default=False, guitype='boolbox', row=6, col=1, rowspan=1, colspan=1, mode='alignment,breaksym')
	
	parser.add_argument("--autocenter",type=str, default='',help="""WARNING: Experimental. Default=None. Autocenters each averaged pair during initial average generation with --btref and --hacref. Will also autocenter the average of all particles after each iteration of iterative refinement. Options are --autocenter=xform.centerofmass (self descriptive), or --autocenter=xform.centeracf, which applies auto-convolution on the average.""")
	
	parser.add_argument("--autocentermask",type=str, default='',help="""WARNING: Experimental. Requires --autocenter. Default=None. Masking processor to apply before autocentering. See 'e2help.py processors -v 10' at the command line.""")
	
	parser.add_argument("--autocenterpreprocess",action='store_true', default=False,help="""WARNING: Experimental. Requires --autocenter. Default=False. This will apply a highpass filter at a frequency of half the box size times the apix, shrink by 2, and apply a low pass filter at half nyquist frequency to any computed average for autocentering purposes if --autocenter is provided. Default=False.""")
	
	#parser.add_argument('--reverse_contrast', action="store_true", default=False, help=""" This multiplies the input particles by -1. Remember that EMAN2 **MUST** work with 'white protein' """)
	
	parser.add_argument("--shrink", type=int,default=0,help="Optionally shrink the input volumes by an integer amount for coarse alignment.")
	
	parser.add_argument("--shrinkfine", type=int,default=0,help="Optionally shrink the input volumes by an integer amount for fine alignment.")
	
	parser.add_argument("--search", type=int,default=8,help=""""During COARSE alignment translational search in X, Y and Z, in pixels. Default=8. This WILL overwrite any search: provided through --align, EXCEPT if you provide --search=8, which is the default. In general, just avoid providing search twice (through here and through the aligner, --align). If you do, just be careful to make them consistent to minimize misinterpretation and error.""")
	
	parser.add_argument("--searchfine", type=int,default=2,help=""""During FINE alignment translational search in X, Y and Z, in pixels. Default=2. This WILL overwrite any search: provided through --falign, EXCEPT if you provide --searchfine=2, which is the default. In general, just avoid providing search twice (through here and through the fine aligner --falign). If you do, just be careful to make them consistent to minimize misinterpretation and error.""")
	
	#parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default='', guitype='strbox', row=19, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')

	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default="thread:1", guitype='strbox', row=19, col=0, rowspan=1, colspan=3, mode='alignment,breaksym')
	#parser.add_argument("--automask",action="store_true",help="Applies a 3-D automask before centering. Can help with negative stain data, and other cases where centering is poor.")
	#parser.add_argument("--resample",action="store_true",help="If set, will perform bootstrap resampling on the particle data for use in making variance maps.",default=False)
	#parser.add_argument("--odd", default=False, help="Used by EMAN2 when running eotests. Includes only odd numbered particles in class averages.", action="store_true")
	#parser.add_argument("--even", default=False, help="Used by EMAN2 when running eotests. Includes only even numbered particles in class averages.", action="store_true")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	parser.add_argument("--plots", action='store_true', help="""Turn this option on to generate a plot of the ccc scores both during model generation with e2spt_classaverage.py or e2spt_hac.py and for refinement results of e2spt_refinemulti.py. Running on a cluster or via ssh remotely might not support plotting.""",default=False)
	
	parser.add_argument("--notmatchimgs",action='store_true',default=False,help="""Default=True. This option prevents applying filter.match.to to one image so that it matches the other's spectral profile during preprocessing for alignment purposes.""")
	
	parser.add_argument("--preavgproc1",type=str,default='',help="""Default=None. A processor (see 'e2help.py processors -v 10' at the command line) to be applied to the raw particle after alignment but before averaging (for example, a threshold to exclude extreme values, or a highphass filter if you have phaseplate data.)""")
	
	parser.add_argument("--preavgproc2",type=str,default='',help="""Default=None. A processor (see 'e2help.py processors -v 10' at the command line) to be applied to the raw particle after alignment but before averaging (for example, a threshold to exclude extreme values, or a highphass filter if you have phaseplate data.)""")
	
	parser.add_argument("--weighbytiltaxis",type=str,default='',help="""Default=None. A,B, where A is an integer number and B a decimal. A represents the location of the tilt axis in the tomogram in pixels (eg.g, for a 4096x4096xZ tomogram, this value should be 2048), and B is the weight of the particles furthest from the tomogram. For example, --weighbytiltaxis=2048,0.5 means that praticles at the tilt axis (with an x coordinate of 2048) will have a weight of 1.0 during averaging, while the distance in the x coordinates of particles not-on the tilt axis will be used to weigh their contribution to the average, with particles at the edge(0+radius or 4096-radius) weighing 0.5, as specified by the value provided for B.""")
	
	parser.add_argument("--weighbyscore",action='store_true',default=False,help="""Default=False. This option will weigh the contribution of each subtomogram to the average by score/bestscore.""")
	
	parser.add_argument("--clipali",type=int,default=0,help="""Default=0 (which means it's not used). Boxsize to clip particles as part of preprocessing to speed up alignment. For example, the boxsize of the particles might be 100 pixels, but the particles are only 50 pixels in diameter. Aliasing effects are not always as deleterious for all specimens, and sometimes 2x padding isn't necessary; still, there are some benefits from 'oversampling' the data during averaging; so you might still want an average of size 2x, but perhaps particles in a box of 1.5x are sufficiently good for alignment. In this case, you would supply --clipali=75""")

	parser.add_argument("--precision",type=float,default=1.0,help="""Default=1.0. Precision in pixels to use when figuring out alignment parameters automatically using --radius. Precision would be the number of pixels that the the edge of the specimen is moved (rotationally) during the finest sampling, --falign. If precision is 1, then the precision of alignment will be that of the sampling (apix of your images) times the --shrinkfine factor specified.""")

	parser.add_argument("--tweak",action='store_true',default=False,help="""WARNING: BUGGY. [NOT IMPLMEMENTED HERE YET]. This will perform a final alignment with no downsampling [without using --shrink or --shrinkfine] if --shrinkfine > 1.""")

	(options, args) = parser.parse_args()

	try:
		hdr = EMData(options.input,0,True) #This is done just to test whether the input file exists where it should
		boxsize = hdr['nx']
	except:
		print """ERROR: Can't find the file provided through --input""", options.input

	if options.radius:
		options.align = None
		options.falign = None
	
	#print "\n\n\n(e2spt_refinemulti.py) options.refpreprocess is", options.refpreprocess
	#print "\n\n\n"
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	
	from e2spt_classaverage import sptmakepath 
	
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

	
	logger = E2init(sys.argv,options.ppid)	
	
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
	
	
	'''
	Store parameters in parameters.txt file inside --path
	'''
	from e2spt_classaverage import writeParameters
	cmdwp = writeParameters(options,'e2spt_refinemulti.py', 'refinemulti')
	
	
	
	'''
	Determine how many references there are and put them into one file classAvg.hdf, or 
	classAvg_iterXX.hdf, if they come from separate files
	'''
	nrefs=0
	
	
	
	from e2spt_classaverage import sptRefGen 
	
	nptcls = EMUtil.get_image_count( options.input )
	
	reffilesrefine = {}
	
	if options.nref == 1:
		print "ERROR: --nref must be > 1. For single reference refinement, use e2spt_classaverage.py"""
		sys.exit(1)
	
	
	
	classmxFile = options.path + '/classmx_' + str( 0 ).zfill( len (str (options.iter))) + '.hdf'
	
	ncls = 0
	if options.classmx:
		classmxFile = options.classmx
		classmx = EMData.read_images( classmxFile )		# we keep the entire classification matrix in memory, since we need to update it in most cases
		#ncls = int(classmx[0]["maximum"])
		ncls = int( classmx[0]['nx'] )	
	
		if options.nref: 
			if options.nref != ncls:
				print """ERROR: It makes no sense to provide both --classmx and --nref that don't 
				match. While --nref is %d, the number of classes specified by the classmx 
				file is %d""" %(options.nref,ncls)
				sys.exit()
				
		if options.ref:
			print """ERROR: It makes no sense to provide both --classmx and --ref. 
			This program PRODUCES classmx files. If you provide a classmx file, it will be 
			used to separate the data in --input into subgroups for initial model generation
			to seed multiple model refinement. It's an alternative to --nref, available in case
			the classmx file has useful information and it "knows" how many classes there 
			should be (opposed to providing an arbitrary value through --nref)."""
			sys.exit()
	
	elif options.nref:
		ncls = options.nref			
			
		
	if options.ref:
		
		reffiles = options.ref.split(',')
		nreffiles = len(reffiles)
		
		print "there are %d external references" %( nreffiles )
		
		if len( reffiles ) == 1:
			if options.nref > 1:
				print "%d copies of the single reference %s being prepared" %( options.nref, options.ref) 
				reffiles *= options.nref
				
				#for i in range( options.nref ):
				#	reffiles.append( options.ref )		
					
			elif options.nref < 2:
				print """ERROR: you've provided a single reference file, and specified 
				the number of references to refine as less than 2 through --nref.
				For single reference iterative refinement/alignment, use
				e2spt_classaverage.py"""
				sys.exit()
				
				
		elif len( reffiles ) > 1:
			if len( reffiles ) == options.nref:
				pass #everything matches, no problems
				
			elif len( reffiles ) > options.nref:
				
				if options.nref < 2:
					print """WARNING: --nref was < 2, but you provided more than one reference
					through --ref, therefore the number of classes to seed refinement will match the number
					of files in --ref."""
					ncls = len( reffiles )
				
				elif options.nref > 1:
					print """WARNING: some files provided through --ref will NOT be used to 
					refine the data, since the number of files in --ref is %d but you've selected
					to refine only %d references""" %( len( reffiles ), options.nref )
					ncls = options.nref
						
			elif len( reffiles ) < options.nref:
				print """WARNING: the number of references provided through --ref, %d, 
				is lower than the number of references requested for refinement --nref, %d""" %( len(reffiles), options.nref )
				print """Therefore, copies will be made from the files provided through 
				--ref, until enough references are seeded to match --nref."""
				ncls = options.nref

				for j in range( options.nref ):
					reffiles *= math.ceil( float(options.nref) / float(len(reffiles) ) )
					reffiles = reffiles[:options.nrefs]
				
		
		print "there are %d references to prepare" % ( len(reffiles) )
		print "from these files", reffiles
		
		
		ptclhdr = EMData( options.input, 0, True)
		ptclnx = ptclhdr['nx']
		ptclny = ptclhdr['ny']
		ptclnz = ptclhdr['nz']
		rr=0
		for rf in reffiles:
			refhdr = EMData( rf, 0, True)
			refnx = refhdr['nx']
			refny = refhdr['ny']
			refnz = refhdr['nz']
			print "ref size and apix are", refnx, refny, refnz, refhdr['apix_x']
			print "ptcl size and apix are", ptclnx, ptclny, ptclnz, ptclhdr['apix_x']
			if ptclnx != refnx or ptclny != refny or ptclnz != refnz:
				print """ERROR: particles in --input=%s, size %d x %d x %d, are not the 
				same size as reference %s of size %d x %d x %d""" % ( options.input, ptclnx,ptclny,ptclnz, rf, refnx,refny,refnz )
				sys.exit(1)
			
			print "preparing ref number", rr
			ptclnumsdummy = {0:[]}
			options.ref = rf
			ret = sptRefGen( options, ptclnumsdummy, cmdwp, rr+1, method=options.refgenmethod ) 	#pass in the reference (or 'class') number by sending rr; #This returns a dictionary with { klass_indx:img } pairs; klass_indx is always zero in this case though, since we're retrieving one reference at a time
			refimg = ret[0]
			print "returned from sptRefGen", refimg
			
			reffile = rootpath + '/' + options.path + '/ref' + str( rr ).zfill( len( str( len( reffiles)))) + '.hdf'
			
			refimg.write_image( reffile, 0 )
			
			reffilesrefine.update( {rr: reffile} )
			
			rr+=1
		
	elif options.classmx:
		pass
		#building initial models using alignment info in classmx file
	
	else:

		ptclnumsdict = {}
		
		groupsize = nptcls / ncls
		
		for i in range( ncls):

			ptclist = [j for j in xrange(groupsize*i, groupsize*(i+1))]	
			if i == ncls - 1:
				ptclist = [j for j in xrange(groupsize*i, nptcls) ]
				#ptclnumsdict.update( { i: ptclist } )
			
			ptclnumsdict.update( { 0: ptclist  } )
		
			print "ptclnumsdict to send is", ptclnumsdict
				
			ret = sptRefGen( options, ptclnumsdict, cmdwp, 1, method=options.refgenmethod, subset4ref=options.subset4ref ) #This returns a dictionary with { klass_indx:img } pairs; klass_indx is always zero in this case though, since we're retrieving one reference at a time
			refimg = ret[0]
			
			reffile = rootpath + '/' + options.path + '/ref' + str( i ).zfill( len( str( ncls ))) + '.hdf'
			
			refimg.write_iamge( reffile, 0 )
		
			reffilesrefine.update( {i:reffile} )


	print "There are these many references", len(reffilesrefine)	
	#sys.exit()
	
	nclstest = len( reffilesrefine )
	print "nclstest %d should match ncls %d" %(nclstest,ncls)
		
	
	'''
	if classmxfiles do not exist, generate/seed them
	'''
	
	#weights = [1.0]*nptcls
	if not options.classmx:
		options.classmx = classmxFile
			
		#C: classmx images are put into a single stack of 2-D images with 9 images in it (0-8)
		#C: and written out so that the file exists when accessed later by the code
			
		classmxScores = EMData(ncls,nptcls)
		classmxWeights = EMData(ncls,nptcls)
		classmxXs = EMData(ncls,nptcls)
		classmxYs = EMData(ncls,nptcls)
		classmxZs = EMData(ncls,nptcls)
		classmxAzs = EMData(ncls,nptcls)
		classmxAlts = EMData(ncls,nptcls)
		classmxPhis = EMData(ncls,nptcls)
		classmxScales = EMData(ncls,nptcls)

		classmxScores.to_zero()
		
		groupsize = nptcls / ncls
		
		for ii in range(ncls):
			ptclist = [jj for jj in xrange(groupsize*ii, groupsize*(ii+1))]	
			if ii == ncls - 1:
				ptclist = [j for j in xrange(groupsize*ii, nptcls) ]
			
			
			#ptclnumsdict.update( { i: ptclist } )
			
			klassid = ii 
		
		
			for p in ptclist:				
				print "\n(e2spt_refinemulti)(main) - Particle %d will belong to classid %d" %( p, klassid )
				classmxScores.set_value_at( klassid, p, 1.0 )
		
		classmxWeights.to_one() 	#Particles contribute entirely and equally to the class to which they are assigned
		classmxXs.to_zero()
		classmxYs.to_zero()
		classmxZs.to_zero()
		classmxAzs.to_zero()
		classmxAlts.to_zero()
		classmxPhis.to_zero()
		classmxScales.to_one()		#One means no scaling
	
		classmxScores.write_image(classmxFile,0)
		classmxWeights.write_image(classmxFile,1)
		classmxXs.write_image(classmxFile,2)
		classmxYs.write_image(classmxFile,3)
		classmxZs.write_image(classmxFile,4)
		classmxAzs.write_image(classmxFile,5)
		classmxAlts.write_image(classmxFile,6)
		classmxPhis.write_image(classmxFile,7)	
		classmxScales.write_image(classmxFile,8)
			
		print "\n(e2spt_refinemulti)(main) - classmx files initialized."
		
	
	'''
	Initialize parallelism
	'''
	if options.parallel:

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

	
	
	'''
	Refine the data against each reference and loop over iterations
	'''
	for it in range( options.iter ):
		print "\n\nrefining all references in iteration", it
		print "\n\n"
				#The program goes into --path to execute the alignment command; therefore, --input will be one level furtherback
	
		if it > 0:
			newreffiles = {}				
			#avgsName =  'classAvgs.hdf'
			#if options.savesteps and int( options.iter ) > 1 :
			#	avgsName = 'classAvgs_iter' + str().zfill( len( str( options.iter ))) + '.hdf'
	
			newreffile = ''
			for reftag in avgs:
				if avgs[ reftag ] and avgs[reftag] != None and avgs[reftag]!='None':
					newref = avgs[ reftag ]
					
					for refindx in reffilesrefine:
						ref = reffilesrefine[refindx]
						if reftag in ref:
							if '_iter' in ref:
								newreffile = ref.split('_iter')[0] + '_iter' + str(it).zfill( len( str( options.iter ))) + '.hdf'
							else:
								newreffile = ref.replace('.hdf','_iter' + str(it).zfill( len( str( options.iter ))) + '.hdf')
						
							newref.write_image( newreffile, 0)
		
				
				newreffiles.update({ reftag: newreffile } )
			
			reffilesrefine = newreffiles
			
			if len(reffilesrefine) < 2:
				finalize = 1
				print "(e2spt_refinemulti.py, line 268) All particles preferred one average and therefore multirefine has failed/converged."
				sys.exit()
				


		reftags = []
		masterInfo = {}
		
		if options.syms and it == 0:
			options.syms = options.syms.split(',')
		
			if len(options.syms) != len(reffilesrefine):
				if len(options.syms) > len(reffilesrefine):
					options.syms = options.syms[0,len(reffilesrefine)]
				elif len(options.syms) < len(resfiles):
				
					howMany = len(reffilesrefine) - len(options.syms)
				
					for pi in range(howMany):
						options.syms.append('c1')
			
		
		
		classmxFile = options.path + '/classmx_' + str( it ).zfill( len (str (options.iter))) + '.hdf'
		
		ic = 0
		for refindx in reffilesrefine:
			
			#results = refineref ( options, reffilesrefine[refindx], nptcls, it )
			
			
			tasks = []

			results = []
			transform = None
			# loop over volumes
	
			#ref.write_image(os.path.join(options.path,"tmpref.hdf"),0)
			reffile = reffilesrefine[refindx]
			'''
			set up tasks
			'''
			for i in range( nptcls ):
				ptclnum = i
				if options.parallel:
					task=Align3DTask( ["cache", reffile , 0], ["cache",options.input,ptclnum],ptclnum, "ptcl %d in iter %d" % (ptclnum, it), options, transform, it )
					tasks.append(task)
				else:
					#print "No parallelism specified"
					result=align3Dfunc( ["cache", reffile , 0], ["cache",options.input,ptclnum],ptclnum, "Ptcl %d in iter %d" % (ptclnum, it), options, transform, it )
					results.append(result['final'])


			'''
			start alignments (execute tasks)
			'''
			if options.parallel:
				tids=etc.send_tasks(tasks)
				if options.verbose: 
					print "%d tasks queued in seedtree level %d"%(len(tids),i) 

				"""Wait for alignments to finish and get results"""
				results = get_results(etc,tids,options.verbose, nptcls ,'refinemulti')
			
			
				
				
			'''
			Add info per particle for results from all references to a master dictionary, 'masterInfo',
			and write out results to .json database and classmx files
			'''
			
			reftag = str( refindx ).zfill( len( str( len( reffilesrefine ))))
			reftags.append( reftag )
			
			itertag = str( it ).zfill( len( str( options.iter )))
			
			'''
			Define and open the .json dictionaries where alignment and score values will be stored, for each iteration,
			and for each reference if using multiple model refinement
			'''
			scoresFile = originalCompletePath +'/sptali_rm_' + itertag + '_' + reftag + '.json'
			jsScores = js_open_dict(scoresFile)			
			
			
			'''
			Iterate over alignment results to write them to classmx.hdf and .json files
			'''
			print "len results is", len(results)
			print "should match nptcls", nptcls
		
			print "results are", results
			iii = 0
		
			classScoresList = [] 
			for r in results:
				
				ptclindx = r[1]
				
				score = r[0][0]['score']
				if options.verbose > 3:
					print "for particle %d score is %.4f" %(ptclindx,score)
				
				classScoresList.append(score)
			
				t = r[0][0]['xform.align3d']
			
				if options.verbose > 3:
					print "and transform is",t
		
				ptclID = 'subtomo_' + str( ptclindx ).zfill( len( str( nptcls ) ) )			
				aliparams = [t,score]
				jsScores.setval( ptclID, aliparams )
				
		
				if options.verbose > 3:
					print "wrote info to .json file"
			
			
				#ptclID = "tomo_" + str(i).zfill( len(str( nptcls )) )
				ptclScore = float( aliparams[-1] )
				ptclAliParams = aliparams
			
				
				infolist = [ ptclScore, ptclAliParams, reftag]
				if ic==0:									#In the first iteration and for the first reference, the ptclID key does not exist; create a dummy to update below in further ic and it iterations.
					masterInfo.update({ ptclID: [] })
			
				print "\n\niteration %d, reference %d, ptclID to update is %s" %( it,ic,ptclID)
				print "infolist to append is", infolist
				print "BEFORE appending, masterInfo[ ptclID ] is", masterInfo[ ptclID ]
				print "Of type", type(masterInfo[ ptclID ])
			
				value = masterInfo[ ptclID ]				#Retrieve results from previous alignments against previous references
				
				
				value.append(infolist)						#Append results from aligning against current reference 
			
				#print "Therfore value is", value
				#print "Of type", type(value)
				#print "\n\n"
				
				masterInfo.update({ ptclID: value })		#Update masterInfo with expanded results. In the end, a particle key under ptclID should have a list with results for as many references as there are
				
				print "AFTER appending and updating, masterInfo[ ptclID ] is", masterInfo[ ptclID ]
				
	
		
			ic+=1
		
		
		
		print "reftags are", reftags
			
		'''
		Analyze all results and classify particles based on them
		'''			
		from operator import itemgetter						
	
		print "I've aligned all the particles in the data set to all the references for iter %d and will now classify them from the masterInfo dict" %(it), masterInfo
	
		classes = {}			#classes dictionary to store alignment info per class, depending on which reference each particle preferred
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
		print "(e2spt_refinemulti.py) there are these many classes with particles",len(classes)
		print "classes are", classes
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
				print "\n\nsending particles in class %d to averaging" %( klassIndx ) 
				print classes[klass]
				print "\n\n"
				ret = makeAverage( options, classes[klass], klassIndx, klassesLen, it, finalize, originalCompletePath)
					
				#ret = makeAverage(options, ic,results,it)
				ref = ret[0]
				weights = ret[1]
				
				print "returned weights are", weights
				
				if ref:
					avgsName =  originalCompletePath + '/class_avgs.hdf'
					if options.savesteps and int( options.iter ) > 1 :
						avgsName = originalCompletePath + '/class_avgs_iter' + str( it ).zfill( len( str( options.iter ))) + '.hdf'						
					ref.write_image(avgsName,-1)
			else:
				print "The klass %d was empty (no particles were assgined to it). You might have too many classes." % ( klassIndx )	
				#dummyClassAvg=EMData(boxsize,boxsize,boxsize)
				#dummyClassAvg.to_zero()
				ret = None
				ref = None
				weights = None
			
			avgs.update({ klass : ref })
		
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
		
		#print "gout to set value for classmxweights at", ic, ptclindx, weights[ptclindx]
		#classmxWeights.set_value_at(ic,ptclindx,weights[ptclindx])
			
		#print "set value for classmxweights at", ic, ptclindx, weights[ptclindx]
		
	
	E2end(logger)
	
	return





def refineref ( options, reffile, nptcls, it ):

	'''
	#Outer loop covering levels in the converging binary tree
	'''
	#for i in range( nptcls ):
		
	#print "Outfile will be", outfile

	

	return results
	
	


	
def makeAverage(options, klass, klassIndx, klassesLen, iterNum, finalize, originalCompletePath):
	"""Will take a set of alignments and an input particle stack filename and produce a new class-average.
	Particles may be excluded based on the keep and keepsig parameters. If keepsig is not set, then keep represents
	an absolute fraction of particles to keep (0-1). Otherwise it represents a sigma multiplier akin to e2classaverage.py"""
	path = options.path
	weights = {}
	'''	
	writeali = 0
	aliptcls = path + '/aliptcls' + klassid + '.hdf'
	
	weights={}
	
	try:
		if options.saveallalign:
			writeali = 1
			aliptcls = path + '/aliptcls' + klassid + '_' + str(it).zfill( len(str(options.iter)) ) + '.hdf'

		elif saveali and it == options.iter - 1:
			writeali = 1
	except: #The exception should be triggered when e2spt_hac.py is called since it doesn't have the --iter parameter.
		if options.saveali:
			writeali = 1		#This saves the aligned particles across ALL iterations for HAC -probably shouldn't be done.
	'''
	
			
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
	
	avgr = Averagers.get(options.averager[0], options.averager[1])
	included = []
	
	#print "The path to save the class average is", options.path
			
	#jsdict = path + '/tomo_xforms.json'
	#js = js_open_dict(jsdict)
			
	ptclsAdded = 0
	for k in klass:
		print "klass is", klass
		weight = 1.0
		if klass and len(klass) > 0:
			#print "\n\nk in klass is", k
			#print "\nThe index of the particle to add is",k[0]
			#print "\nAnd this its transform", k[1][0]
			ptclindx = k[0]
			ptcl = EMData(options.input,ptclindx)
		
			#print "\n\n\n(e2spt_refinemuti.py) in makeAverage, ptcl and its type are",ptcl,type(ptcl)
		
			ptclTransform =k[1][0]
			#print "And the ptcl transform is", ptclTransform
			ptcl.process_inplace("xform",{"transform" : ptclTransform})
		
			#print "I've applied the transform"
			#print "I have applied this transform before averaging", ptcl_parms[0]["xform.align3d"]			
		
			if k[-1] <= thresh: 
			
				if options.weighbytiltaxis:
					px = x = int(ptcl['ptcl_source_coord'][0])
				
					tiltaxis = int( options.weighbytiltaxis.split(',')[0] )
					minweight = float( options.weighbytiltaxis.split(',')[1] )
				
					if px > tiltaxis:
						px = -1 *( px - 2*tiltaxis )	#This puts te particle at the same distance from te tilt axis, but to the left of it.
						
					X = tiltaxis				#This models a line in 'weight space' (x, w), that passes through (0, minweight) and ( tiltaxis, maxweight ) 
					W = 1.0 - minweight
					slope = W/X
											#Having the slope of the line and its y-axis (or w-axis in this case) crossing we predict the weight of any particle depending on its dx distance to the tiltaxis
					print "Tiltaxis is", X
					print "W is", W
					print "Therefore slope is", slope
				
					dx = tiltaxis - px 
					#if px > tiltaxis:
					#	dx = px - tiltaxis
						
					taweight = slope * px + minweight 
					weight = weight * ( taweight )
					print "tiltaxis weight was %f because it's distance from the tilt axis is %d, because it's x coordinate was %d" % (taweight, dx, x)

				if options.weighbyscore:
					scoreweight = score / maxscore
					print "the score weight is %f because score was %f and the best score was %f" % (scoreweight, score, maxscore )
					weight = weight * scoreweight
			
				weights.update( {ptclindx:weight} )
				
				print "therefore the final weight for particle %d is %f" %(ptclindx, weight )
				
				ptcl.mult( weight )
			
			
			
			
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
	avg["spt_multiplicity"] = len(included)
	avg['spt_ptcl_indxs']=included
	
	#if options.averager[0] == 'mean' and variance:
	#	variance.write_image(path+"/class_varmap.hdf",it)
				
	if options.autocenter:
		print "\n\n\n\nYou have selected to autocenter!\n", options.autocenter
		
		avgac = avg.copy()
		if options.autocentermask:
			avgac.process_inplace( options.autocentermask[0],options.autocentermask[1] )
			
		if options.autocenterpreprocess:
			apix = avgc['apix_x']
			halfnyquist = apix*4
			highpassf = apix*a['nx']/2.0
			
			avgac.process_inplace( 'filter.highpass.gauss',{'cutoff_freq':highpassf,'apix':apix})
			avgac.process_inplace( 'filter.lowpass.gauss',{'cutoff_freq':halfnyquist,'apix':apix})
			avgac.process_inplace( 'math.meanshrink',{'n':2})
			
		avgac.process_inplace(options.autocenter[0],options.autocenter[1])
		
		tcenter = avgac['xform.align3d']
		print "Thus the average HAS BEEN be translated like this", tcenter
		avg.transform(tcenter)

	avg['origin_x']=0
	avg['origin_y']=0
	avg['origin_z']=0
	
	
	if ptclsAdded > 0:
		#avg.write_image(avgsName,klassIndx)
		
		return [avg,weights]
	else:
		return None,None



	
if __name__ == '__main__':
	main()

