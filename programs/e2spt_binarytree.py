#!/usr/bin/env python

#
# Author: Jesus Galaz-Montoya 01/Nov/2014
# Last modification: 16/Sept/2014
#
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
	
	parser.add_argument("--align",type=str,help="""This is the aligner use for alignments. 
		Default is rotate_translate_3d:search=8:delta=12:dphi=12""", default="rotate_translate_3d:search=8:delta=12:dphi=12")

	parser.add_argument("--aligncmp",type=str,help="""The comparator used for the --align aligner. 
		Default is the internal tomographic ccc. Do not specify unless you need to use another specific aligner.""",default="ccc.tomo")	
	
	parser.add_argument("--autocenter",type=str, default='',help="""Autocenters each averaged pair on
		all rounds. Options are --autocenter=xform.centerofmass (self descriptive), or
		--autocenter=xform.centeracf, which applies auto-convolution on the average.
		Default=None.""")
	
	parser.add_argument("--autocentermask",type=str, default='',help="""Masking processor 
		to apply before autocentering. Default=None.
		'e2help.py processors -v 10' at the command line.""")
	
	parser.add_argument("--autocenterpreprocess",action='store_true', default=False,help="""This 
		will apply a highpass filter at a frequency of half the box size times the apix, 
		shrink by 2, and apply a low pass filter at half nyquist frequency to any computed
		average for autocentering purposes if --autocenter is provided. Default=False.""")
	
	parser.add_argument("--averager",type=str,help="""The type of averager used to produce the class average. 
		Default=mean""",default="mean")
	
	parser.add_argument("--clipali",type=int,default=0,help="""Boxsize to clip particles as part of preprocessing
		to speed up alignment. For example, the boxsize of the particles might be 100 pixels, but the particles are only 50 pixels 
		in diameter. Aliasing effects are not always as deleterious for all specimens, and sometimes 2x padding isn't necessary;
		still, there are some benefits from 'oversampling' the data during averaging; so you might still want an average of size
		2x, but perhaps particles in a box of 1.5x are sufficiently good for alignment. In this case, you would supply --clipali=75""")
	
	parser.add_argument("--falign",type=str,help="""This is the second stage aligner used to refine the first alignment. 
		Default is refine_3d_grid:range=12:delta=4:search=2, specify 'None' to disable""", default="refine_3d_grid:range=12:delta=4:search=2")
	
	parser.add_argument("--faligncmp",type=str,help="""The comparator used by the second stage aligner. 
		Default is the internal tomographic ccc""",default="ccc.tomo")
		
	parser.add_argument("--highpass",type=str,help="""A processor (as in e2proc3d.py; 
		could be masking, filtering, etc.) to be applied to each volume prior to alignment. 
		NOT applied to aligned particles before averaging.""",default=None)
	
	parser.add_argument("--highpassfine",type=str,help="""A highpass filtering processor (as in e2proc3d.py) 
		to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.""", default=None)
	
	parser.add_argument("--input", type=str, help="""The name of the input volume stack. 
		MUST be HDF or BDB, since volume stack support is required.""", default='')
		
	parser.add_argument("--iter", type=int, help="""The number of iterations to perform. 
		Default is automatically determined depending on the size of the input stack.""", default=0)
	
	parser.add_argument("--lowpass",type=str,help="""A processor (as in e2proc3d.py; 
		could be masking, filtering, etc.) to be applied to each volume prior to alignment. 
		NOT applied to aligned particles before averaging.""",default=None)
		
	parser.add_argument("--lowpassfine",type=str,help="""A lowpass filtering processor (as in e2proc3d.py) 
		to be applied to each volume prior to FINE alignment. Not applied to aligned particles before averaging.""", default=None)
	
	parser.add_argument("--mask",type=str,help="""Mask processor applied to particles before alignment. 
		Default is mask.sharp:outer_radius=-2""", default="mask.sharp:outer_radius=-2")
		
	parser.add_argument("--maskfile",type=str,default=None,help="""Mask file (3D IMAGE) applied to particles 
		before alignment. Must be in HDF format. Default is None.""")

	parser.add_argument("--normproc",type=str,help="""Normalization processor applied to particles before alignment. 
		Default is to use 'normalize'. If normalize.mask is used, results of the mask option will be passed in automatically. 
		If you want to turn this option off specify \'None\'""", default="normalize")

	parser.add_argument("--npeakstorefine", type=int, help="""The number of best coarse alignments 
		to refine in search of the best final alignment. Default=4.""", default=4)
	
	parser.add_argument("--nseedlimit",type=int,default=0,help="""Maximum number of particles
		to use. For example, if you supply a stack with 150 subtomograms, the program will
		automatically select 128 as the limit to use because it's the largest power of 2 that is
		smaller than 150. But if you provide, say --nseedlimit=100, then the number of particles
		used will be 64, because it's the largest power of 2 that is still smaller than 100.""")
	
	parser.add_argument("--parallel",  help="""Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel""", default="thread:1")
	
	parser.add_argument("--path",type=str,default='',help="""Directory to store results in. 
		The default is a numbered series of directories containing the prefix 'spthac'; 
		for example, spthac_02 will be the directory by default if 'spthac_01' already exists.""")
		
	parser.add_argument("--plotccc", action='store_true', help="""Turn this option on to generate
		a plot of the ccc scores for all comparisons for the FIRST iteration of all vs all.
		Running on a cluster or via ssh remotely might not support plotting.""",default=False)
	
	parser.add_argument("--postprocess",type=str,help="""A processor to be applied to the volume 
		after averaging the raw volumes, before subsequent iterations begin.""",default=None)
			
	parser.add_argument("--ppid", type=int, help="""Set the PID of the parent process, used for cross platform PPID""",default=-1)
	
	parser.add_argument("--precision",type=float,default=1.0,help="""Precision in pixels to use
		when figuring out alignment parameters automatically using --radius. Precision 
		would be the number of pixels that the the edge of the specimen is moved (rotationally) during the 
		finest sampling, --falign. If precision is 1, then the precision of alignment will be that of 
		the sampling (apix of your images) times the --shrinkfine factor specified.""")
	
	parser.add_argument("--preprocess",type=str,help="""A processor (as in e2proc3d.py; could be masking, filtering, etc.) 
		to be applied to each volume prior to alignment. Not applied to aligned particles before averaging.""",default=None)
	
	parser.add_argument("--preprocessfine",type=str,help="""Any processor (as in e2proc3d.py) 
		to be applied to each volume prior to FINE alignment. NOT applied to aligned particles before averaging.""", default=None)
	
	parser.add_argument("--procfinelikecoarse",action='store_true',default=False,help="""
		If you supply this parameters, particles for fine alignment will NOT be preprocessed
		identically to particles for coarse alignment by default.  
		If you supply this, but want specific parameters for preprocessing particles for 
		also supply: fine alignment, nd supply fine alignment parameters, such as 
		--lowpassfine, --highpassfine, etc; to preprocess the particles for FINE alignment 
		differently than for COARSE alignment.""")
	
	parser.add_argument("--radius", type=float, help="""Hydrodynamic radius of the particle in Angstroms. 
		This will be used to automatically calculate the angular steps to use in search of the best alignment.
		Make sure the apix is correct on the particles' headers, sine the radius will be converted from Angstroms to pixels.
		Then, the fine angular step is equal to 360/(2*pi*radius), and the coarse angular step 4 times that""", default=0)
	
	parser.add_argument("--randomizewedge",action="store_true", help="""This parameter is EXPERIMENTAL. 
		It randomizes the position of the particles BEFORE alignment, to minimize missing wedge bias 
		and artifacts during symmetric alignment where only a fraction of space is scanned""", default=False)
	
	parser.add_argument("--search", type=int,default=8,help=""""During COARSE alignment
		translational search in X, Y and Z, in pixels. Default=8.
		This WILL overwrite any search: provided through --align,
		EXCEPT if you provide --search=8, which is the default. In general, just avoid
		providing search twice (through here and through the aligner, --align). If you do,
		just be careful to make them consistent to minimize misinterpretation and error.""")
	
	parser.add_argument("--searchfine", type=int,default=2,help=""""During FINE alignment
		translational search in X, Y and Z, in pixels. Default=2.
		This WILL overwrite any search: provided through --falign,
		EXCEPT if you provide --searchfine=2, which is the default. In general, just avoid
		providing search twice (through here and through the fine aligner --falign). If you do,
		just be careful to make them consistent to minimize misinterpretation and error.""")
	
	parser.add_argument("--savesteps",action="store_true", help="""If set, this will save 
		the averages after each iteration.""", default=False)
	
	#parser.add_argument("--saveali",action="store_true", help="""If set, this will save the 
	#	aligned/averaged volumes from the immediately PREVIOUS round 
	#	that went into the NEW particles in the "current" round, to round#_particles.hdf. 
	#	It will also save the latest state of alignment (for the LAST iteration only) of ALL particles provided in the input stack.
	#	Overwrites existing files.""",default=False)

	parser.add_argument("--sym", dest = "sym", default='c1', help = """Symmetry to impose - choices are: c<n>, d<n>, h<n>, tet, oct, icos""")
		
	parser.add_argument("--shrink", type=int,default=1,help="""Optionally shrink the input volumes 
		by an integer amount for coarse alignment.""")
	
	parser.add_argument("--shrinkfine", type=int,default=1,help="""Optionally shrink the input volumes 
		by an integer amount for refine alignment.""")
	
	parser.add_argument("--threshold",type=str,help="""A threshold applied to the subvolumes after normalization. 
		For example, --threshold=threshold.belowtozero:minval=0 makes all negative pixels equal 0, so that they do not contribute to the correlation score.""", default=None)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="""verbose level [0-9], higner number means higher level of verboseness""")

	parser.add_argument("--savepreprocessed",action="store_true", help="""Will save stacks 
		of preprocessed particles (one for coarse alignment and one for fine alignment if 
		preprocessing options are different).""", default=False)	

	(options, args) = parser.parse_args()
	
	'''
	If --radius of the particle is provided, we calculate the optimal alignment steps for 
	coarse and fine alignment rounds using --shrink and --shrinkfine options and apix info
	'''
	if options.radius:
		from e2spt_classaverage import calcAliStep
		options = calcAliStep(options)
		
	
	if options.align:
		if options.sym and options.sym is not 'c1' and options.sym is not 'C1' and 'sym' not in options.align and 'grid' not in options.align:
			if 'rotate_translate_3d' in options.align or 'rotate_symmetry_3d' in options.align: 
				options.align += ':sym=' + str( options.sym )
						
		if 'search' not in options.align:
			if 'rotate_translate_3d' in options.align:		
				options.align += ':search=' + str( options.search )
			
			#sys.exit()
			
		elif 'rotate_translate_3d' in options.align:
			searchA = options.align.split('search=')[-1].split(':')[0]
			searchdefault = 8
			
			#print "There was search in --align", searchA
			#sys.exit()
			if options.search != searchdefault:
						
				prefix = options.align.split('search=')[0]
				trail = options.align.split('search=')[-1].split(':')[-1]
			
				options.align =  prefix + 'search=' + str(options.search)
				if len(trail) > 2 and '=' in trail:
					options.align += ':' + trail 
			
				print """\nWARNING: --search is different from search= provided through
				--align or its default value of 8. There's no need to specify both, 
				but if you did, --search takes precedence :-) ."""
				#sys.exit()
			elif options.search == searchdefault:
				options.search = searchA
				

	
		
		if "rotate_translate_3d_grid" in options.align:
			if "alt0" and "alt1" in options.align:
				alt0 = int(options.align.split('alt0')[-1].split(':')[0].replace('=',''))	
				alt1 = int(options.align.split('alt1')[-1].split(':')[0].replace('=',''))
				
				print "alt0 and alt1 are", alt0,alt1, type(alt0), type(alt1)
				print alt1-alt0 == 0
				#sys.exit()
				
				if alt1-alt0 == 0:
					print """\nERROR: alt0 and alt1 cannot be equal for rotate_translate_3d_grid.
					If you want to inactivate searches in this angle, provide a alt0 and alt1
					such that alt1-alt0 is NOT ZERO, and provide a step size for dalt that is larger
					than this difference. For example: 
					alt0=0:alt1=1:dalt=2."""
					sys.exit()
					
			if "phi0" and "phi1" in options.align:
				phi0 = int(options.align.split('phi0')[-1].split(':')[0].replace('=',''))	
				phi1 = int(options.align.split('phi1')[-1].split(':')[0].replace('=',''))
				
				print "phi0 and phi1 are", phi0,phi1, type(phi0), type(phi1)
				print phi1-phi0 == 0
				#sys.exit()
				
				if phi1-phi0 == 0:
					print """\nERROR: phi0 and phi1 cannot be equal for rotate_translate_3d_grid.
					If you want to inactivate searches in this angle, provide a phi0 and phi1
					such that phi1-phi0 is NOT ZERO, and provide a step size for dphi that is larger
					than this difference. For example: 
					phi0=0:phi1=1:dphi=2."""
					sys.exit()
					
			if "az0" and "az1" in options.align:
				az0 = int(options.align.split('az0')[-1].split(':')[0].replace('=',''))	
				az1 = int(options.align.split('az1')[-1].split(':')[0].replace('=',''))
				
				print "az0 and az1 are", az0,az1, type(az0), type(az1)
				print az1-az0 == 0
				#sys.exit()
				
				if az1-az0 == 0:
					print """\nERROR: az0 and az1 cannot be equal for rotate_translate_3d_grid.
					If you want to inactivate searches in this angle, provide a az0 and az1
					such that az1-az0 is NOT ZERO, and provide a step size for daz that is larger
					than this difference. For example: 
					az0=0:az1=1:daz=2."""
					sys.exit()
	
	
	if options.falign and options.falign != None and options.falign != 'None' and options.falign != 'none':
		if 'search' not in options.falign and 'refine_3d_grid' in options.falign:	
			options.falign += ':search=' + str( options.searchfine )
			
		else:
			searchF = options.falign.split('search=')[-1].split(':')[0]
			searchfinedefault = 2
			
			if options.searchfine != searchfinedefault:
						
				prefix = options.falign.split('search=')[0]
				trail = options.falign.split('search=')[-1].split(':')[-1]
			
				options.falign =  prefix + 'search=' + str(options.searchfine)
				if len(trail) > 2 and '=' in trail:
					options.falign += ':' + trail 
			
				print """\nWARNING: --searchfine is different from search= provided through
				--falign or its default value of 2. There's no need to specify both, but 
				if you did, --searchfine takes precedence :-) ."""
				#sys.exit()
				
			elif options.searchfine == searchfinedefault:
				options.searchfine = searchF
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'spt_bt')
	
	'''
	Store parameters in parameters.txt file inside --path
	'''
	from e2spt_classaverage import writeParameters
	writeParameters(options,'e2spt_binarytree.py', 'bt')
	
	'''
	Parse parameters
	'''
	
	if options.autocenter:
		options.autocenter=parsemodopt(options.autocenter)
		
	if options.autocentermask:
		options.autocentermask=parsemodopt(options.autocentermask)

	if options.align:
		options.align=parsemodopt(options.align)
	
	if options.falign and options.falign != None and options.falign != 'None' and options.falign != 'none': 
		options.falign=parsemodopt(options.falign)
	
	if options.aligncmp: 
		options.aligncmp=parsemodopt(options.aligncmp)
	
	if options.faligncmp: 
		options.faligncmp=parsemodopt(options.faligncmp)
	
	if options.averager: 
		options.averager=parsemodopt(options.averager)		
		
	if options.normproc and options.normproc != 'None' and options.normproc != 'none':
		options.normproc=parsemodopt(options.normproc)	
		
	if options.mask and options.mask != 'None' and options.mask != 'none':
		#print "parsing mask", sys.exit()
		options.mask=parsemodopt(options.mask)
	
	if options.preprocess and options.preprocess != 'None' and options.preprocess != 'none': 
		options.preprocess=parsemodopt(options.preprocess)
		
	if options.threshold and options.threshold != 'None' and options.threshold != 'none': 
		options.threshold=parsemodopt(options.threshold)
		
	if options.preprocessfine and options.preprocessfine != 'None' and options.preprocessfine != 'none': 
		options.preprocessfine=parsemodopt(options.preprocessfine)
		
	if options.lowpass and options.lowpass != 'None' and options.lowpass != 'none': 
		options.lowpass=parsemodopt(options.lowpass)
		
	if options.lowpassfine and options.lowpassfine != 'None' and options.lowpassfine != 'none': 
		options.lowpassfine=parsemodopt(options.lowpassfine)
	
	if options.highpass and options.highpass != 'None' and options.highpass != 'none': 
		options.highpass=parsemodopt(options.highpass)
		
	if options.highpassfine and options.highpassfine != 'None' and options.highpassfine != 'none': 
		options.highpassfine=parsemodopt(options.highpassfine)
		
	if options.postprocess and options.postprocess != 'None' and options.postprocess != 'none': 
		options.postprocess=parsemodopt(options.postprocess)
		
	
	
	if options.shrink < options.shrinkfine:
		options.shrink = options.shrinkfine
		print "It makes no sense for shrinkfine to be larger than shrink; therefore, shrink will be made to match shrinkfine"
					
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
		
	ref = binaryTreeRef(options,nptclForRef,nseed,-1,etc)

		
	print "Will end logger"	
	E2end(logger)
	
	print "logger ended"
	sys.stdout.flush()
	
	return
	

def binaryTreeRef(options,nptclForRef,nseed,ic,etc):
	
	from e2spt_classaverage import *
	
	if nptclForRef==1: 
		print "Error: More than 1 particle required if no reference provided through --ref."
		sys.exit(1)
			
	# we need to make an initial reference. Due to the parallelism scheme we're using in 3-D and the slow speed of the
	# individual alignments we use a slightly different strategy than in 2-D. We make a binary tree from the first 2^n particles and
	# compute pairwise alignments until we get an average out. 

	
	
	#if nseed>64 : 
	#	nseed=64
	#	print "Limiting seeding to the first 64 images"

	nseediter=int(log(nseed,2))			# number of iterations we'll need
	
	if options.iter:
		nseediter=options.iter
	
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
	for i in range(nseediter):
		infile="%s/seedtree_%d_cl_%d.hdf"%(options.path,i,ic)
		if ic < 0:
			infile="%s/seedtree_%d.hdf"%(options.path,i)
			
		print "Infile will be", infile
		
		outfile="%s/seedtree_%d_cl_%d.hdf"%(options.path,i+1,ic)
		if ic < 0:
			outfile="%s/seedtree_%d.hdf"%(options.path,i+1)
	
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
			results=get_results(etc,tids,options.verbose,{},nseed,0,'binarytree')

			#results=get_results(etc,tids,options.verbose,{},len(ptclnums),0,'binarytree')

			#results=get_results(etc,tids,options.verbose,{},nptclForRef,0)



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
		
	ref=EMData(outfile,0)		# result of the last iteration
	
	return ref
	

def makeAveragePairs(options,ptcl_file,outfile,align_parms):
	"""Will take a set of alignments and an input particle stack filename and produce a new set of class-averages over pairs"""
	
	current = os.getcwd()
	print "\n(e2spt_classaverage.py) (make_average_pairs) current directory is", current
	findir = os.listdir(current)
	print "\noptions.path is", options.path
	findirpath = os.listdir(options.path)
	print "\nThe particle file where the particles ought to be read from is", ptcl_file
	print "\nLets see if ptcl_file is in path. Files in path are", findirpath
	
	print "\nalign_parms are", align_parms
	print "\nTheir len", len(align_parms)
	
	for i,ptcl_parms in enumerate(align_parms):
		
		print "\ni, ptcl_parms are", i, ptcl_parms
		
		#if ptcl_parms:
		ptcl0=EMData(ptcl_file,i*2)
		
		ptcl1=EMData(ptcl_file,i*2+1)
		ptcl1.process_inplace("xform",{"transform":ptcl_parms[0]["xform.align3d"]})
	
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
	
		avg.write_image(outfile,i)
	return


if __name__ == '__main__':
	main()

