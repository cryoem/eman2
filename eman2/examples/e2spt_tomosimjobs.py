#!/usr/bin/env python

#
# Author: Jesus Galaz-Montoya, 2012. Last update: july/24/2014.
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

from optparse import OptionParser
from EMAN2 import *
from sys import argv
import EMAN2
import heapq
import operator
import random
import numpy
import colorsys
import os
import math

from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import matplotlib
matplotlib.use('Agg',warn=False)

import matplotlib.pyplot as plt

from e2spt_classaverage import sptmakepath


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <output> [options]

	This program takes in a model in .hdf format, calls e2spt_simulation.py to generate a simulated set of subtomograms from it,
	and characterizes the ability of EMAN2's e2classaverage3d.py to align the particles under varying Single Particle Tomography
	parameters.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	#parser = OptionParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--path",type=str,default=None,help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsimjob'; for example, sptsimjob_02 will be the directory by default if 'sptsimjob_01' already exists.")
	
	parser.add_argument("--snrlowerlimit", type=float,default=0.00,help="Minimum weight for noise compared to singal.")
	parser.add_argument("--snrupperlimit", type=float,default=1.00,help="Maximum weight for noise compared to singal.")

	parser.add_argument("--snrchange", type=float,default=1.00,help="""Step to vary snr from one run to another. 
										For example, if this parameter is set to 2.0, snr will be tested from --snrlowerlimit (for example, 0.0), increasing by --snrchange, 0.0,2.0,4.0... up to --snrupperlimit.""")

	parser.add_argument("--tiltrangelowerlimit", type=int,default=60,help="""Minimum value for imaging range (at a value of 90, there's no missing wedge. 
											60 would mean the data will be simulated, as if it came from a tomogram that was reconstructed from a til series
											collected from -60deg to 60deg)""")
	parser.add_argument("--tiltrangeupperlimit", type=int,default=61,help="""Maximum value to simulate the imaging range. Simulations will be run starting at --tiltrangelowerlimit, and will increase
											by --tiltrangestep, until the last simulation is done at --tiltrangeupperlimit.""")
	parser.add_argument("--tiltrangechange", type=int,default=1,help="Amount (in degrees) to decrease the size of the missing wedge from one run to another.")	
	
	
	#parser.add_argument("--tiltsteplowerlimit", type=int,default=1,help="""Within each tiltrange simulated, you can simulate individual pictures taken with different tilt steps.
	#									For example, if you collect images from -60deg to 60deg with a 2deg tilt step, the tilt series will have 61 images.
	#									If, on the other hand, the tilt step was 4deg, the tilt series will only have 31 images.
	#									--tiltstepupperlimit is the largest step size you want to simulate.""")
	#parser.add_argument("--tiltstepupperlimit", type=int,default=2,help="""Within each tiltrange simulated, you can simulate individual pictures taken with different tilt steps.
	#									For example, if you collect images from -60deg to 60deg with a 2deg tilt step, the tilt series will have 61 images.
	#									If, on the other hand, the tilt step was 4deg, the tilt series will only have 31 images.
	#									--tiltstepupperlimit is the largest step size you want to simulate.""")
	#parser.add_argument("--tiltstepchange", type=int,default=1,help="""Increase in size of tilt step from one run to another. 
	#									Jobs will be run using --tiltstepstep as the first value, and then adding that value on subsequent runs until --tiltsteplimit is reached""")
	
	parser.add_argument("--tiltstep", type=int,default=None,help="""This is the angular distance between different tilts. 
									Should be supplied opposed to --nslices parameters, when variations in missing wedge size want to be tested, 										keeping a constant distance between tilts. --nslices should be used when the wedge size, determined by --tiltrange 										parameters, is kept constant.""")									
	
	parser.add_argument("--nsliceslowerlimit", type=int,default=None,help="Lowest number of slices to divide the tiltrange in. If on and --nslicesupperlimit is ALSO on (any number other than zero), this will turn off all tiltstep parameters.")	

	parser.add_argument("--nslicesupperlimit", type=int,default=None,help="Largest number of slices to divide the tiltrange in. If on and --nsliceslowerlimit is ALSO on (values other than zero), this will turn off all tiltstep parameters.")	

	parser.add_argument("--nsliceschange", type=int,default=None,help="Change with which the nslices parameter in e2spt_simulation.py will be varied. Will only work if --nslicesupperlimit and --nsliceslowerlimit are different than zero.")
	
	#parser.add_argument("--simstack",type=str,default=None,help="If you already have an .hdf stack of subtomograms generated with e2spt_simulation.py, you can provide it here to test for different alignment parameters.")

	"""
	Parameters to be passed on to e2spt_simulation.py
	"""
	
	#parser.add_argument("--interpolator",default='',help="""What interpolation scheme 
	#	to use for reconstruction. Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 
	#	'gauss_5_slow', 'gypergeom_5', 'experimental' """)
	
	parser.add_argument("--fillwedge",action='store_true',default=False,help="""This option will
		fill the region of the missing wedge with evenly spaced images of Gaussian noise
		with the same tiltstep used to simulate the particles (as determined through the
		parameter --nslices).""")
		
	parser.add_argument("--model", type=str, default='',help="""The name of the input volume from which 
							simulated subtomograms will be generated. 
							The output will be in HDF format, since volume stack support is required. 
							The input CAN be PDB, MRC or and HDF stack. 
							If the input file is PDB or MRC, a version of the supplied 
							model will be written out in HDF format.
							If the input file is a stack, simulatd subvolumes will be 
							generated from each model in the stack and written to different output stacks.
							For example, if the input file contains models A and B, 
							two output stacks with simulated subvolumes will be generated.""")
							
	parser.add_argument("--ptclstack", type=str, default='', help="""You can provide a stack
		of simulated particles to skip generation of the simulated stack. This can be useful
		in testing the alignment of the stack under different particle preprocessing 
		(--lowpass/--testlowpassfreq, --highpass, --mask/--testmaskradius, --shrinkalign/--testsampling) 
		and alignment (--align, --falign, --npeakstorefine) parameters. 
		""")
		
	parser.add_argument("--randstack", type=str, default='', help="""You can provide a stack
		of randomly oriented particles. These still need to be simulated though. 
		This can be useful in testing the alignment of the SAME particles under different 
		particle simulation (tilt range, tilt step, snr) and alignment parameters (shrink, mask, etc). 
		""")
				
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_argument("--nptcls", type=int,default=10,help="Number of simulated subtomograms tu generate per referece.")
	parser.add_argument("--tx", type=int,default=0,help="""Maximum number of pixels to randomly translate each subtomogram in X. The random translation will be picked between -txrange and +txrange. 
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--ty", type=int,default=0,help="""Maximum number of pixels to randomly translate each subtomogram in Y. The random translation will be picked between -tyrange and +tyrange.
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--tz", type=int,default=0,help="""Maximum number of pixels to randomly translate each subtomogram in Z. The random translation will be picked between -tzrange and +tzrange.
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--transrange", type=int,default=4,help="""Maximum number of pixels to randomly translate each subtomogram in all X, Y and Z. 
									The random translation will be picked between -transrage and +transrange; --txrange, --tyrange and --tzrange overwrite --transrange for each specified direction.""")
	
	parser.add_argument("--search", type=int,default=4,help="""Maximum number of pixels 
		to search translationally during alignment in all three X, Y and Z directions.""")
	
	
	
	parser.add_argument("--applyctf", action="store_true",default=False,help="If on, it applies ctf to the projections in the simulated tilt series based on defocus, cs, and voltage parameters.")

	parser.add_argument("--defocus", type=int,default=3,help="Intended defocus at the tilt axis (in microns) for the simulated tilt series.")
	parser.add_argument("--voltage", type=int,default=200,help="Voltage of the microscope, used to simulate the ctf added to the subtomograms.")
	parser.add_argument("--cs", type=int,default=2,help="Cs of the microscope, used to simulate the ctf added to the subtomograms.")

	parser.add_argument("--gridholesize", type=int,default=2,help="""Size of the carbon hole for the simulated grid (this will determine shifts in defocus for each particle at 
									each tilt step, depending on the position of the particle respect to the tilt axis, which is assigned randomly.""")
	parser.add_argument("--saverandstack", action="store_true",default=False,help="Save the stack of randomly oriented particles, before subtomogram simulation (before the missing wedge and noise are added).")
	
	
	parser.add_argument("--nosim", action="store_true",default=False,help="""If on
		the program will generate stacks of "perfect particles" in different random 
		orientations, but with no missing wedge, no noise, no ctf parameters, etc.""") 
		
	
	parser.add_argument("--saveprjs", action="store_true",default=False,help="Save the projections (the 'tilt series') for each simulated subtomogram.")

	parser.add_argument("--reconstructor", type=str,default="fourier",help="""The reconstructor 
		to use to reconstruct the tilt series into a tomogram. Type 'e2help.py reconstructors' 
		at the command line to see all options and parameters available.
		To specify the interpolation scheme for the fourier reconstruction, specify 'mode'.
		Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 
		'gauss_5_slow', 'gypergeom_5', 'experimental'.
		For example --reconstructor=fourier:mode=gauss_5 """)

	parser.add_argument("--pad", type=int,default=0,help="""If on, it will increase the box size of the model BEFORE generating projections and doing 3D reconstruction of simulated sutomograms.""")								
	
	parser.add_argument("--finalboxsize", type=int,default=0,help="""The final box size to clip the subtomograms to.""")								

	parser.add_argument("--snr",type=int,help="Weighing noise factor for noise added to the image.",default=0)
	#parser.add_argument("--addnoise",action="store_true",default=False,help="If on, it adds random noise to the particles")
	
	parser.add_argument("--sym",type=str,default='c1',help="If your particle is symmetrical, you should randomize it's orientation withing the asymmetric unit only. Thus, provide the symmetry.")

	parser.add_argument("--notrandomize",action="store_true",default=False,help="This will prevent the simulated particles from being rotated and translated into random orientations.")
	parser.add_argument("--simref",action="store_true",default=False,help="""This will make a simulated particle in the same orientation as the original input (or reference), 
																			and will be used as the model/reference to test alignments if --testalignment is on.
																			To have e2spt_classaverage.py make its own reference, simply do not supply this parameter.""")
	parser.add_argument("--negativecontrast",action="store_true",default=False,help="This will make the simulated particles be like real EM data before contrast reversal. Otherwise, 'white protein' (positive density values) will be used.")

	parser.add_argument("--testalignment",action="store_true",default=False,help="This will run e2spt_classaverage.py to test the alignment of the particles against the simulated reference.")

	parser.add_argument("--files",type=str,default='',help="""Text files to be plotted, with information in the correct format. To enter multiple files, separate them with commas: file1,file2,etc...
								Each file must contain lines with the following values:
								tr= ts= snr= angular_error= translational_error=""")
	parser.add_argument("--plotonly",action="store_true",default=False,help="""Skip simulation and generate plots based on .txt files of each set of parameters.
										You must also provide --files or --filesdir.""")
	parser.add_argument("--filesdir",type=str,default='',help="The directory where the files to analyze are")
	
	"""
	Parameters to be passed on to e2spt_classaverage.py
	"""
	parser.add_argument("--faligncmp",type=str,default='ccc.tomo',help="Comparator to use for missing wedge compensation during fine alignment.")
	parser.add_argument("--aligncmp",type=str,default='ccc.tomo',help="Comparator to use for missing wedge compensation during coarse alignment.")
	
	parser.add_argument("--comparators",type=str,default='',help="""Supply comparators to COMPARE, separated my a comma. For example --comparators=ccc,ccc.tomo,fsc.tomo:sigmas=0.005 .
																	If supplied, --faligncmp and --aligncmp will be ignored.""")
	
	parser.add_argument("--lowpass",type=str,help="""A lowpass filter (as in a processor 
		from e2proc3d.py) apply to the particles during alignment. Type 'e2help.py processors' 
		at the command line to find the options available from the processors list)""",default='')
	
	parser.add_argument("--highpass",type=str,help="""A highpass filter (as in a processor 
		from e2proc3d.py) to apply to the particles before alignment. Type 'e2help.py processors' 
		at the command line to find the options available from the processors list)""",default='')
	
	parser.add_argument("--mask",type=str,help="""A mask (as in a processor from e2proc3d.py)
		to apply to the particles before alignment. Type 'e2help.py processors' 
		at the command line to find the options available from the processors list)""",default='')
	
	parser.add_argument("--shrinkalign", type=int,default=0,help="Optionally shrink the input volume before the simulation if you want binned/down-sampled subtomograms.")
	parser.add_argument("--shrinksim", type=int,default=0,help="Optionally shrink the input volume before the simulation if you want binned/down-sampled subtomograms.")
	
	parser.add_argument("--npeakstorefine", type=int, help="The number of best coarse alignments to refine in search of the best final alignment. Default=1.", default=1)
	parser.add_argument("--iter",type=str,default='1',help="Number of iterations of alignment. If you want to test MULTIPLE iteration values, supply a comma separated list; e.g., --iter=1,2,3,4 or --iter=2,4,8 ...etc.")
	
	parser.add_argument("--parallel",type=str,default='thread:1',help="Parallelization to use.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--normproc",type=str,help="""Normalization processor applied to 
		particles before alignment. Default is to use None. If normalize.mask is used, 
		results of the mask option will be passed in automatically. If you want to turn 
		this option on specify any normalizing processor such as 'normalize' or '
		normalize.edgemean'. To see all processors available, type e2help.py propcessors
		at the commnad line.""", default='')


	#parser.add_argument("--fitwedgepost", action="store_true", help="Fit the missing wedge AFTER preprocessing the subvolumes, not before, IF using the fsc.tomo comparator for --aligncmp or --faligncmp.", default=False)


	parser.add_argument("--radius", type=float, help="""Hydrodynamic radius of the particle in Angstroms. 
													This will be used to automatically calculate the angular steps to use in search of the best alignment.
													Make sure the apix is correct on the particles' headers, sine the radius will be converted from Angstroms to pixels.
													Then, the fine angular step is equal to 360/(2*pi*radius), and the coarse angular step 4 times that""", default=0)

	
	parser.add_argument("--precision",type=float,default=1.0,help="""Precision in pixels to use
		when figuring out alignment parameters automatically using --radius. Precision 
		would be the number of pixels that the the edge of the specimen is moved (rotationally) during the 
		finest sampling, --falign. If precision is 1, then the precision of alignment will be that of 
		the sampling (apix of your images) times the --shrinkfine factor specified.""")
	
	#parser.add_argument("--wedgeangle",type=float,help="Missing wedge angle",default=60)
	#parser.add_argument("--wedgei",type=float,help="Missingwedge begining", default=0.10)
	#parser.add_argument("--wedgef",type=float,help="Missingwedge ending", default=0.9)
	parser.add_argument("--writewedge", action="store_true", help="Write a subvolume with the shape of the fitted missing wedge if --faligncmp or --aligncmp are fsc.tomo. Default is 'True'. To turn on supply --writewedge", default=False)

	parser.add_argument("--quicktest",action="store_true",default=False,help="This will run e2spt_classaverage.py with minimal parameters to quickly test the program.")
	
	parser.add_argument("--testsampling",type=str,default='',help="""Enter a number for 
		the maximum shrinking factor to test. For example, if your particles have an apix 
		of 2angstroms/pixel and you enter --testsampling=4, your particle will be shrunk by 
		1,2,3 and 4 in independent runs of the program (1 needs no shrinking, but still is
		the initial point to test the alignments) and thus alignment will be tested at
		samplings of 2,4,6, and 8 A/pixel (corresponding to shrink factors of 1,2,3 and 4).
		This parameter will overwrite the --shrinkalign parameter.""")
	
	parser.add_argument("--testlowpassfreq",type=str,default='',help="""
		You have to provide the lowpass processor through --lowpass for this parameter to work.
		Enter three numbers in angstroms: the lightest lowpass to start with, the harshest limit 
		for lowpassing, and the step with which to vary the lowpass.
		For example, you can test lowpass filtering between 10 and 50 angstroms, every 5
		angstroms, by providing --testlowpasssfreq=10,50,5""")

	parser.add_argument("--testboxsize",type=str,default='',help="""Enter two numbers:
		1) The largest boxsize to use, and a positive increment.
		OR
		2) The smallest boxisze to use, and a negative decrement.
		In case 1, it is assumed that the particle is in a 'tight box' with no padding,
		and increasingly large boxes will be tested. For example, if you're particles/model are/is 
		in a box of size 64, and you supply --testboxisze=128,4, successive boxes will be
		tested starting at 64 and up to 128: 64,68,72,76,... ...128.
		In case 2, it is assumed that the particles/model are/is already padded, therefore smaller
		boxes will be tested, as follows. For example, if the box is 128, and you supply:
		--testboxsize=64,-4, the boxes to test will be: 128, 124, 120... ...64.
		""")
		
	parser.add_argument("--testmaskradius",type=str,default='',help="""
		--mask has to be supplied for this parameter to work.
		This assumes the particles/model are/is already adequately padded (that is, 
		the box size is ~2x the size of the molecule it contains).
		The mask's radius is expressed in pixels from the edge, and starts at -2.
		Provide two numbers: The smallest radius to use, and the decrement.
		For example, if your particles/model are/is in a box of 128, and you provide
		--testmaskradius=-64,-4, the masks tried will be: -2,-6,-10,-14 ... ... -62.""")
	
	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)
	
	rootpath = os.getcwd()
	
	'''
	If PERFORMING alignments to compute ALIGNMENT ERROR
	'''
	if not options.plotonly:
		
		skip=0
		if options.comparators:
			if len( options.comparators.split(',') ) > 1:
				skip=1
		
		if not skip:
		
			options = sptmakepath( options, options.path )
		
		#print "This is the path RETURNED after the first path made", options.path
		
		rootpath = os.getcwd() + '/' + options.path
		originalpathbase = options.path
		
		randstack=''
		if not options.ptclstack:
			
			#if options.randstack:
			#	os.system('e2proc3d.py ' + options.randstack + ' ' + rootpath + '/model.hdf --first=0 --last=0')
			#	options.model = rootpath + '/model.hdf'
			#	
			#elif options.model:	
			#	os.system('cp ' + options.model + ' ' + rootpath + '/' + os.path.basename(options.model) )
			
			wildcard = ''
			if options.testsampling or options.testboxsize or options.testlowpassfreq or options.testmaskradius:
				wildcard = 'tmp'
			
			
			if options.randstack:
				randstack = options.randstack
			
			ret=simloop(options,rootpath,randstack,wildcard)
		
			randstack= ret[10]
			simstack= ret[-3][0]
			#print "\n\n\n\nSize of randomized randstack is", EMData(randstack,0,True)['nx']
			#print "Randstack is", randstack
			#print "!@#$%^&*()_+(*&^%$#@!@#$%^&*()(*&^%$#@!@#$%^&*()(*&^%$#@#$%^&*(*&^%$#@#$%^&*()*&^%$#@#$%^&*\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
		
			#if options.path not in randstack:
			#	randstack = options.path + '/' + randstack
			
			#if rootpath not in randstack:
			#	randstack = rootpath + '/' + randstack
				
			#if options.path not in simstack:
			#	randstack = options.path + '/' + randstack
			
			#if rootpath not in simstack:
			#	randstack = rootpath + '/' + randstack
		
			#print "\n\n(e2spt_tomosimjobs.py) in the first round,randstack fixed after returned is", randstack
			#print "simstack fixed after returned is", simstack
			
		else:
			randstack=options.ptclstack
			#print "You've provided randstack via ptclstack", options.ptclstack
			
		
		if randstack and (options.testsampling or options.testboxsize or options.testlowpassfreq or options.testmaskradius):
			
			
			#print "\n\n\n\nInitial randstack for testboxsize is", randstack
			#print "And rootpath is", rootpath
			
			#randstack = randstack.replace(options.path,options.path + '/' + options.path)
			#print "Randstack should be", randstack
			
			filesindir = os.listdir( os.getcwd() )
			if originalpathbase not in filesindir:
				os.system('mkdir ' + rootpath)
			
			
			randstack = rootpath + '/' + options.path + '/randstack.hdf'
			
			options.model = rootpath + '/simmodel.hdf'
			#print "\n\n\n\nSize of randomized randstack is", EMData(randstack,0,True)['nx']
			#print "Randstack is", randstack
			#print "!@#$%^&*()_+(*&^%$#@!@#$%^&*()(*&^%$#@!@#$%^&*()(*&^%$#@#$%^&*(*&^%$#@#$%^&*()*&^%$#@#$%^&*\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
			
			#print "Actual randstack after adding rootpath is", randstack
			#print "Because rootpath plus randstack are", rootpath, randstack, rootpath + randstack
			#print "\n\n\n\n"
			
			
			if options.ptclstack:
				randstack=options.ptclstack
				#print "You've provided randstack via ptclstack", options.ptclstack
			
			os.system( 'e2proc3d.py ' + randstack + ' ' + rootpath + '/randstack.hdf')
						
			os.system( 'e2proc3d.py ' + randstack + ' ' + options.model + ' --first=0 --last=0' )
	
			if options.ptclstack:
				randstack = rootpath + '/randstack.hdf'
			
					
			nrefs=1
			
			nslices = ''
			if options.tiltstep:
				nslices = ((2 * options.tiltrangelowerlimit) / options.tiltstep ) + 1						#You have to add 1 to account for the slice at 0deg
				#nslicesu = nslices + 1																#Add 1 because the loop will exclude the upper limit
			elif options.nsliceslowerlimit:
				nslices = options.nsliceslowerlimit		
			else:
				print "ERROR: You have to supply either tiltrange parameters including --tiltstep, or nslices parameters"
			
			if not nslices:
				print "ERROR: Couldn't figure out nslices"
				sys.exit()
			
			
			tiltrange = options.tiltrangelowerlimit
			snr = options.snrlowerlimit
			
			tiltrangetag = ("%d" %( tiltrange ) ).zfill(3)
			nslicestag = ("%d" %( nslices ) ).zfill(3)
			snrtag = ("%.2f" %( snr) ).zfill(5)
			
			
			simround = 1
			samestackformany = 1
			firstrandstack = randstack
			
			thesemodels = [ options.model ]
			thesestacks = [ randstack ]
			
			'''
			TEST PADDING BY CLIPPING
			'''
			if options.testboxsize:
				limit = int(options.testboxsize.split(',')[0])
				boxstep = int(options.testboxsize.split(',')[-1])
				
				boxes=[]
				currentbox = int(EMData(randstack,0,True)['nx'])
				if boxstep < 0:
					boxes = [ x for x in xrange( limit, currentbox + -1*boxstep, int(math.fabs( boxstep )) ) ]
				elif boxstep > 0:
					boxes = [ x for x in xrange( currentbox, limit + 1*boxstep, int(math.fabs( boxstep )) ) ]
				
				#print "\n\n\n\nThe boxes to use are", boxes
				#print "Current box is", currentbox
				#print "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBVVVVBVVVVBBBBBBBBBBBBBB\n\n\n\n"
				if boxes:
					for box in boxes:
						
						options.finalboxsize = box
						
						options.path = 'bx' + str(box).zfill(3)
						
						#os.system('mkdir ' + rootpath + '/' + options.path)
						#wildcard = options.path
						
						gencmds(options,rootpath,nrefs,tiltrangetag,tiltrange,nslicestag,nslices,snrtag,snr,simround,firstrandstack,samestackformany=0,thesestacks=[],thesemodels=[],itersMode=0,wildcard=options.path)

				if options.testalignment:
					analyzeresults(options,rootpath,nrefs,{},'testbx')
			
			
			'''
			TEST SAMPLING BY SHRINKING
			'''
			if options.testsampling:
				shrinkvals = [ x for x in range(options.testsampling) ]
				
				for shrink in shrinkvals:
					if shrink > 1:
						if options.finalboxsize:
							options.finalboxsize /= shrink
					
						options.shrinkalign = shrink
					
					options.path = 'sh' + str(shrink).zfill(2)
					#wildcard = options.path
					
					gencmds(options,rootpath,nrefs,tiltrangetag,tiltrange,nslicestag,nslices,snrtag,snr,simround,firstrandstack,samestackformany=0,thesestacks=[],thesemodels=[],itersMode=0,wildcard=options.path)
				
				
				
				if options.testalignment:
					analyzeresults(options,rootpath,nrefs,{},'testsh')
			
			
			'''
			'''
			if options.testmaskradius:
				pass
			
				if options.testalignment:
					analyzeresults(options,rootpath,nrefs,{},'testlp')
			
			'''
			'''
			if options.testlowpassfreq:
				pass
						
				if options.testalignment:
					analyzeresults(options,rootpath,nrefs,{},'testhp')

			
			
			
		
		'''
		If ONLY PLOTTING previously computed alignment error files
		'''		
	elif options.plotonly:
		print "\n\nI'm in PLOTONLY \n\n"
		
		'''
		Make the directory where to store the results. Sptakepath exists in e2spt_classaverage.py
		'''
		
		
		
		if options.filesdir:
			options.path = options.filesdir
		else:
			options = sptmakepath(options,'spt_tomosimjobs')
		
		files=[]
		if options.files:
			files=options.files.split(',')
		
		elif options.filesdir:
			findir = os.listdir(options.filesdir)
			for f in findir:
				if 'TR' in f and 'SNR' in f and 'NS' in f and '.txt' in f:
					print "Found a file to analyse", f
					files.append(f)
		else:
			print "ERROR: If --plotonly is provided, then you also must supply EITHER --files or --filesdir"
			sys.exit()		
			
		if files:
			resultsdir = rootpath + '/' + options.path
			files.sort()
			print "\n\nAnd these the files", files
			print "\n\n"
			extra=''
			resfiles_analysis(options,files,resultsdir,extra,0)
		
		
	E2end(logger)
	return()
	

def simloop(options,rootpath,randstack,wildcard=''):
	
	

	rootpathOrig = rootpath
	
	compsPATHS = {}
	comps = options.comparators.split(',')
	
	for c in comps:
		compsPATHS.update( { c:'' } )
	
	#compIDS={}
	#if comps and len(comps) >1:
	#	for comp in comps:
	#		compID = comp.split(':')[0].replace('.','P')
	#		compIDS.update({ compID : comp })
	
	iters=str(options.iter).split(',')
	iterPATHS={}
	
	
	#print "ITers are", iters
	#print "Therefore len(iters) is", len(iters)
	
	if not options.comparators or len(comps) == 1:
		if len(iters) < 2:
			if not wildcard:
			#	print "No wildcard ." + wildcard + " Will make path directly"
			#	options = sptmakepath(options,'spt_tomosimjob')
			#else:
			#	print "Wildcard. " + wildcard + " .Will make path in rootpath"
			#	options = sptmakepath(options,'spt_tomosimjob',rootpath)

			#print "\n\n\n\n\n\n\n\nPPPPPPPPPPPPP\nPath made was", options.path
		
				pass
				
	#print "\nInside simloop, options.path is", options.path

	snrl = options.snrlowerlimit
	snru = options.snrupperlimit
	snrch = options.snrchange

	tiltrangel = options.tiltrangelowerlimit
	tiltrangeu = options.tiltrangeupperlimit
	tiltrangech = options.tiltrangechange
	
	tiltstep = options.tiltstep
	
	if options.nsliceschange and options.nsliceslowerlimit and options.nslicesupperlimit:
		nslicesl = options.nsliceslowerlimit
		nslicesu = options.nslicesupperlimit
		nslicesch = options.nsliceschange
		nslices = nslicesl
		
		if options.tiltstep:
			print """ERROR: You must supply EITHER parameters for --nslices change and lower and upper limits, OR --tiltstep. 
			DO NOT supply parameters for both."""
			sys.exit()
			
	elif not options.tiltstep:
		print """ERROR: You must supply parameters EITHER for --nslices change and lower and upper limits, OR --tiltstep."""
		sys.exit()
	
	#print "\n\nORIGINAL options.model is", options.model
	
	nrefs = 1
	
	if not options.randstack and not options.ptclstack and options.model:
		nrefs = EMUtil.get_image_count( os.path.basename(options.model) )
		
	originalpath = options.path
	
	simround=0
	firstrandstack = ''

	if randstack:
		simround = 1
		firstrandstack = randstack
	
	itersMode = 0
	if not options.simref:
		itersMode = 1
	
	tiltrange = tiltrangel
	while tiltrange <tiltrangeu:
		#print "tiltrage is", tiltrange
		tiltrangetag = ("%d" %(tiltrange) ).zfill(3)
		
		if options.tiltstep:
			#print "Tilt step is", tiltstep
			#print "Tilt range is", tiltrange
			nslices = ((2 * tiltrange) / tiltstep ) + 1						#You have to add 1 to account for the slice at 0deg
			#print "Therefore, 2*tiltrange/tiltstep is nslices",nslices
			nslicesu = nslices + 1											#Add 1 because the loop will exclude the upper limit
			#nslicestag = str(int(nslices)).zfill(3)

		while nslices < nslicesu:
			#print "nslices is", nslices
			#if not options.tiltstep:
			nslicestag = ("%d" %(nslices)).zfill(3)
			snr=snrl
			
			while snr < snru:
				#print "The conditions to simulate are tiltrange=%d, nslices=%d, snr=%.2f" % (tiltrange,nslices,snr)
				#print "Snr is", snr
				
				snrtag = ("%.2f" %(snr) ).zfill(5).replace('.','p')
				
				samestackformany=0
				thestack=[]
				themodel=[]
				
				if len(iters) > 1 and options.comparators:
					print """ERROR: You can't test MULTIPLE iteration numbers, AND MULTIPLE comparators at the same time.
							To test a SINGLE comparator with MULTIPLE ITERATION NUMBERS, supply the comparator via --aligncmp and --faligncmp."""
					sys.exit()
				
				if options.comparators:
					print "\nComparators are",options.comparators
					if len(comps) > 1:
						print "len(comps) is", len(comps)
						for comp in comps:
							print "Comparator is", comp
							print "Snr is", snr
							print "Simround is", simround
							#print "Whereas root is", originalpath
							
							options.aligncmp = comp
							options.faligncmp = comp
							
							#compID = comp.split(':')[0].replace('.','P')
							#if simround == 0:
							#options.path = originalpath + compID
							
							compID = comp.split(':')[0].replace('.','P')
							
							if simround == 0 or (simround == 1 and options.randstack):										#In the first round, the original path
								#pathstem = options.path
								
								rootpath = rootpathOrig + compID
								options.path = originalpath + compID
								
								try:
									options = sptmakepath( options, options.path )
									print "\n\n\n\n\n\n\n\nSimround is 0 and therefore I made path", options.path
								except:
									pass
									
								rootpath = rootpathOrig + options.path.replace(originalpath,'')
								
								compsPATHS.update( { comp: options.path } )
							else:
								options.path = compsPATHS[comp]
								rootpath = rootpathOrig + options.path.replace(originalpath,'')
								
							if options.randstack:
								os.system('e2proc3d.py ' + options.randstack + ' ' + rootpath + '/model.hdf --first=0 --last=0')
								options.model = rootpath + '/model.hdf'			
							
							if not options.randstack and options.model:
								os.system('cp ' + options.model + ' ' + rootpath)
								
								options.model = rootpath + "/" + options.model.split('/')[-1]
						
							nrefs = EMUtil.get_image_count( options.model )
							
							#compPATHS.update({ compID : options.path })
							
							print "THE NEW PATH IS",  options.path
							print "because compID is", compID
							print "because comp is", comp
							
							
							#elif simround > 0:
							#	print "compPATHS are", compPATHS
							#	options.path = compPATHS[ compID ]
							
							#print "\n\n$$$$$$$$$$\nThe stack thestack BEFORE sending is", thestack

							itersMode=0
							
							#ret=['    ','    ','    ']
							ret=gencmds(options,rootpath,nrefs,tiltrangetag,tiltrange,nslicestag,nslices,snrtag,snr,simround,firstrandstack,samestackformany,thestack,themodel,itersMode,wildcard)
						
							if samestackformany == 0:
								thestack = ret[0]
								themodel = ret[1]
							
							if simround == 0:
								firstrandstack = ret[2]
						
							print "Samestackformany is", samestackformany
							
							samestackformany+=1
					
					elif len(iters) == 1:
						options.faligncmp = comps[0]
						options.aligncmp = comps[0]
						
						itersMode = 0
						
						
						if options.randstack:
							os.system('e2proc3d.py ' + options.randstack + ' ' + rootpath + '/model.hdf --first=0 --last=0')
							options.model = rootpath + '/model.hdf'			
					
						if not options.randstack and options.model:
							os.system('cp ' + options.model + ' ' + rootpath)
					
							#if "/" not in options.model:								
							
							options.model = rootpath + "/" + options.model.split('/')[-1]
				
				
						nrefs = EMUtil.get_image_count( options.model )

						
						
						#print "\n\n\n\nThere's comps but length was less than one!"
						ret=gencmds(options,rootpath,nrefs,tiltrangetag,tiltrange,nslicestag,nslices,snrtag,snr,simround,firstrandstack,samestackformany,thestack,themodel,itersMode,wildcard)
						
						#if samestackformany == 0:
						#	thestack = ret[0]
						#	themodel = ret[1]
							
						if simround == 0:
							firstrandstack = ret[2]
						
						
						#if simround == 0:
						#	firstrandstack = ret[2]
				
					'''
					elif len(iters) > 1:
					
						itersMode = 1
					
						itersint =[ int(it) for it in iters]
						itmax = max(itersint)			
						itersint.sort()
					
						#iters = [ str(it) for it in itersint]
					
						itround = 0
						for it in itersint:
							#print "ITERATIONS TO DO ARE", it
						
							if itround == 0:
								options.iter = str(it)
							else:
								options.iter = str(itersint[itround] - itersint[itround-1])
						
							print "\nThe literal iterations to do were supposed to be", it
							print "And because itround is", itround
							print "The actual iterations necessary are", options.iter	
						
							iterID = 'ITERS' + str(it).zfill( len( str(itmax)))
						
							if simround == 0:
								options.path = originalpath + iterID
								pathstem = options.path
							
								options = sptmakepath( options, pathstem )
								options.path = rootpath + '/' + options.path
							
								iterPATHS.update({ iterID : options.path })
							
								print "THE NEW PATH IS",  options.path
								print "because compID is", iterID
								print "because iter is", it
						
							elif simround > 0:
								print "iterPATHS are", iterPATHS
								options.path = iterPATHS[ iterID ]
					
						
							#options.path = originalpath + '_' + itID
						
							print "The path to use is", options.path
		
							ret=gencmds(options,rootpath,nrefs,tiltrangetag,tiltrange,nslicestag,nslices,snrtag,snr,simround,firstrandstack,samestackformany,thestack,themodel,itersMode,wildcard)
						
							if samestackformany == 0:
								thestack = ret[0]
						
							themodel = ret[1]
						
							if simround == 0:
								firstrandstack = ret[2]
				
							print "Samestackformany is", samestackformany
						
							samestackformany+=1
						
							itround += 1
					'''	
				else:
					itersMode = 0
					
					#if samestackformany == 0:
					#	thestack = ret[0]
					#	themodel = ret[1]
							
					
					
					if options.randstack:
						os.system('e2proc3d.py ' + options.randstack + ' ' + rootpath + '/model.hdf --first=0 --last=0')
						options.model = rootpath + '/model.hdf'			
					
					if not options.randstack and options.model:
						os.system('cp ' + options.model + ' ' + rootpath)
					
						#if "/" not in options.model:								
							
						options.model = rootpath + "/" + options.model.split('/')[-1]
				
				
					nrefs = EMUtil.get_image_count( options.model )

					ret=gencmds(options,rootpath,nrefs,tiltrangetag,tiltrange,nslicestag,nslices,snrtag,snr,simround,firstrandstack,samestackformany,thestack,themodel,itersMode,wildcard)

					if simround == 0:
						firstrandstack = ret[2]
					
					
					
					#if simround == 0:
					#	firstrandstack = ret[2]
					#	thestack = ret[0]
				
				simround+= 1
				
				snr += snrch
				print "\n\nThe snr has increased by", snrch
				print "And thus will be for the next round", snr

			if options.tiltstep:
				nslices += tiltstep
			else:
				nslices += nslicesch
		
		tiltrange += tiltrangech
	
	if options.testalignment and wildcard != 'tmp':
	
		print "\n\n\n\n========================================================================================\n\n\n\n"
		print "All ALIGNMENTS ARE DONE! It is now a matter of analyzing the results"
		print "\n\n\n\n========================================================================================\n\n\n\n"

		analyzeresults(options,rootpathOrig,nrefs,iterPATHS)
		print "len(iterPATHS) is", len(iterPATHS)
	
	
	
	print "\n\nIn simloop, these are stack and model", thestack, themodel
	
	return (options,rootpathOrig,nrefs,tiltrangetag,tiltrange,nslicestag,nslices,snrtag,snr,simround,firstrandstack,samestackformany,thestack,themodel,itersMode)






def gencmds(options,rootpath,nrefs,tiltrangetag,tiltrange,nslicestag,nslices,snrtag,snr,simround,firstrandstack,samestackformany=0,thesestacks=[],thesemodels=[],itersMode=0, wildcard=''):
	thisstack=''
	thissimodel=''
	
	#finalAvgPath=''
	
	
	finalAvgs=[]
	
	#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n(e2spt_tomosimjobs.py)(gencmds)\nSSSSSSSSSSSSSSSSS\nsamestackformany is", samestackformany
	#print "SSSSSSSSSSSSSSSSS\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
	
	#modeldir = options.path
	#if rootpath not in modeldir:
	
	#modeldir = rootpath + '/' + modeldir
	modeldir = rootpath
	if options.model and "/" in options.model:
		modelbase = os.path.basename(options.model)
		modeldir = options.model.replace( modelbase , '')
			
	#resultsdir = options.path + '/results'
	#resultsdir = 'results'
	
	if wildcard:
		modeldir = rootpath
		#resultsdir = rootpath + '/results'
	
	#if rootpath not in resultsdir:
	
	resultsdir = rootpath + '/results'

	print "\n\n! ! ! ! ! ! ! ! ! ! ! ! ! ! \n The results dir to make should be", resultsdir
	print "(e2spt_tomosimjobs) (function gencmds) Options path is", options.path
	#print"\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
	
	inputdata = options.model
	
	
	for d in range(nrefs):
		modeltag = ''
		#subpath = rootpath + '/' + options.path + '/' +'TR' + str(tiltrange).zfill(5) + '_TS' + tiltsteptag + '_SNR' + str(snr).zfill(5)
		
		#subpath = options.path + '/' +'TR' + tiltrangetag + '_NS' + nslicestag + '_SNR' + snrtag
		
		subpath = 'TR' + tiltrangetag + '_NS' + nslicestag + '_SNR' + snrtag
		
		
		subdirs = ''

		if rootpath not in subpath:
			#subpath = rootpath + '/' + options.path + '/' +'TR' + tiltrangetag + '_NS' + nslicestag + '_SNR' + snrtag
			subpath = rootpath + '/' + subpath
		
		
		if wildcard:
			#	#if wildcard != 'tmp'
			#	subpath += '_' + wildcard
			#	
			 
			#if wildcard == 'ref':
			#	subpath = wildcard
			
			#else:
			subpath = options.path
		
		print "\n\n\n\nSUBPATH IS!!!!!!!!!!!!!!!!!!!!!!", subpath
		#print "\n\n\n"
		#>>subdirs = os.listdir(rootpath + '/' + options.path)
		
		if wildcard:
			subdirs = os.listdir( rootpath )
		else:
			#os.system('mkdir ' + options.path)
			subdirs = os.listdir( options.path )
			
			
		resultsbase = 'results'
		
		if nrefs > 1:
			modeltag = 'model' + str(d).zfill(len(str(nrefs)))
			#subpath += '_' + modeltag
			
			
			#print "Model tag is", modeltag
			
			modeldir = modeltag
			
			if rootpath not in modeldir:
				modeldir = rootpath + '/' + modeldir
			
			#print "(rootpath) modeldir is", modeldir		
			
			if options.path not in modeldir:
				modeldir = options.path + '/' + modeldir
			
			#print "(path) Modeldir is", modeldir
			
			
			if modeltag not in subdirs:
				os.system('mkdir ' + modeldir)
			
			#print "Subdirs are", subdirs
			
			resultsbase = 'results_' + modeltag

			resultsdir = modeldir + '/' + resultsbase
			
			#print "Results dir is", resultsdir
			#sys.exit()
			
			subdirs = os.listdir(modeldir)
			
			subpath = modeldir + '/' + 'TR' + tiltrangetag + '_NS' + nslicestag + '_SNR' + snrtag #CHANGED
			
			#os.system('mkdir ' + subpath)
			
			#print "I've just made this subpath for this specific model and conditions", subpath
			
			model = EMData(options.model,d)

			newname = modeldir + '/' + inputdata.split('/')[-1].replace('.hdf','_' + modeltag + '.hdf') #CHANGED
			
			model.write_image(newname,0)
			inputdata = newname
			print "\nModel is here", newname
			print '\n'
			
		print "\n\nHowever, for this specific model, the resultsdir to make should be", resultsdir
		#print "And it should be amongst subdirs", subdirs
		
		if resultsbase not in subdirs:
			os.system('mkdir ' + resultsdir)
		
		#>>subtomos =  subpath.split('/')[-1] + '.hdf'
		
		subtomos = 'simptcls.hdf'
		
		cmd1a = cmd1b = ''
		
		if samestackformany < 1:
			thisstack = subtomos
			thissimmodel = inputdata
			
			#print "\n\n\n@@@@@@@@@@@@@@@@@@\nstackformany is", samestackformany
			#print "Therefore I WILL run a job for e2spt_simulation"
			#print "@@@@@@@@@@@@\n\n"
			
			#print "\n(e2spt_tomosimjobs.py) (gencmds function) Simulated particles are here", subtomos
			#print '\n'
			
			#jobcmd = 'e2spt_simulation.py --input=' + inputdata + ' --output=' + subtomos + ' --snr=' + str(snr) + ' --nptcls=' + str(options.nptcls) + ' --nslices=' + str(nslices) + ' --tiltrange=' + str(tiltrange) + ' --transrange=' + str(options.transrange) + ' --pad=' + str(options.pad) + ' --finalboxsize=' + str(options.finalboxsize) + ' --verbose=' + str(options.verbose) + ' --parallel=' + str(options.parallel) + ' --path=' + subpath.split('/')[-1]
			
			
			
			
			jobcmd = 'e2spt_simulation.py --input=' + inputdata +  ' --snr=' + str(snr) + ' --nptcls=' + str(options.nptcls) + ' --nslices=' + str(nslices) + ' --tiltrange=' + str(tiltrange) + ' --transrange=' + str(options.transrange) + ' --pad=' + str(options.pad) + ' --finalboxsize=' + str(options.finalboxsize) + ' --verbose=' + str(options.verbose) + ' --parallel=' + str(options.parallel) 
			
			jobcmd += ' --path=' + subpath.split('/')[-1]
			
			if float(snr) > 0.0:
				jobcmd += ' --savenoise'
		
			if simround < 1:
				jobcmd += ' --saverandstack' + ' --shrink=' + str(options.shrinksim)
				
				#nrefs = EMUtil.get_image_count(inputdata)
				tag=''
				if nrefs>1:
					tag = str(d).zfill(len(str(nrefs)))
					
				#>>firstrandstack =  subpath + '/' + subtomos.replace('.hdf','_randst' + tag + '_n' + str(options.nptcls).zfill(len(str(options.nptcls))) + '.hdf').split('/')[-1]
				#firstrandstack =  subpath.split('/')[-1] + '/' + inputdata.replace('.hdf','_randst' + tag + '_n' + str(options.nptcls).zfill(len(str(options.nptcls))) + '.hdf').split('/')[-1]
				
				firstrandstack =  subpath + '/randstack.hdf'
				
				#print "\n\nThe firstrandstack name in e2spt_tomosimjobs when simround is 0, is", firstrandstack
			
			if simround > 0 and firstrandstack:
				#print "\n\nRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR\nThe randstack to PROVIDE because simround > 1 is", firstrandstack
				#print "\nRRRRRRRRRRRRRRRRRRRR\n\n"
				randstackcmd = ' --randstack=' + firstrandstack
				
				if options.randstack:
					randstackcmd = ' --randstack=../' + firstrandstack
					
					if simround == 1:
						jobcmd += ' --shrink=' + str(options.shrinksim)
					
				jobcmd += randstackcmd
			
			if options.simref:
				jobcmd += ' --simref'
			#if snr:
			#	jobcmd += ' --addnoise'
			if options.saveprjs:
				jobcmd += ' --saveprjs'
			if options.negativecontrast:
				jobcmd += ' --negativecontrast'
			
			if options.notrandomize:
				jobcmd += ' --notrandomize'
				
			if options.nosim:
				jobcmd += ' --nosim'			
			
			if options.reconstructor:
				jobcmd += ' --reconstructor=' + options.reconstructor
			
			if options.fillwedge:
				jobcmd += ' --fillwedge'

			
			cmd1a = 'cd ' + modeldir + ' && ' + jobcmd
			
			print "\n\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! \n "
			print "e2spt_simulation job, or cmd1a is",cmd1a
			print "\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !\n\n"
			
			runcmd( cmd1a )
	
		elif samestackformany > 0:
			#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nMMMMMMMMMMM\nThesemodels when stackformany > 0, and its type, are", thesemodels, type(thesemodels)
			#print "Therefore for this model", d
			#print "THEMODEL, thesemodels[d], would be", thesemodels[d]
			#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nMMMMMMMMMMM\n"
			
			try:
				os.system('mkdir ' + subpath)
			except:
				pass
			#>>cmd1b = 'cp ' + str(thesestacks[d]) + ' ' + subpath+ '/' + subtomos.replace('.hdf','_ptcls.hdf')
			
			cmd1b = 'cp ' + str(thesestacks[d]) + ' ' + subpath+ '/simptcls.hdf'
			
			#>>simmodel = subpath + '/' + inputdata.split('/')[-1].replace('.hdf','_sptsimMODEL_SIM.hdf') 
			
			simmodel = subpath + '/simmodel.hdf'
			
			if itersMode:
				simmodel = subpath + '/' + thesemodels[d].split('/')[-1].replace('.hdf','_previousITER.hdf')
			
			cmd1b = cmd1b + ' && cp ' + str(thesemodels[d]) + ' ' + simmodel
			
			print "\n\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! \n "
			print "e2spt_simulation job, or cmd1b is",cmd1b
			print "\n ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !\n\n"
			
			runcmd( cmd1b )
			
			#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nCCCCCCCCCC\nstackformany is", samestackformany
			#print "Therefore I WILL NOT run a job for e2spt_simulation, and will rather copy the previous MODEL",  str(thesemodels[d])
			#print "TO", simmodel 
			#print "CCCCCCCCCCCCC\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"

		#resultsfiles=[]
		
		finalAvgPath = ''
		print "finalAvgPath is", finalAvgPath
		cmd2 = extractcmd = solutioncmd = rfilecmd = ''
		extractcmd1 = extractcmd2 = refprepcmd = extractcmd0 = finalAvgPath = ''
		
		if options.testalignment and wildcard!='tmp': 
			#print "\n\n$$$$$$$$$\nI will test alignment and for that will cd into SUBPATH", subpath
			#print "This is SIMROUND number", simround 
			#print "$$$$$$$$\n\n"

			#cmd = cmd + ' && cd ' + subpath
			if rootpath not in subpath:
				subpath = rootpath + '/' + subpath
				
			cmd2 = 'cd ' + subpath

			ref=''
			if options.simref:
				#ref = inputdata.split('/')[-1].replace('.hdf','_sptsimMODEL_SIM.hdf')
				ref = 'simmodel.hdf'
				
			else:
				prepRefFromOrigModel = subpath + '/' + options.model.split('/')[-1].replace('.hdf','_REFPREP.hdf')
				#refprepcmd = 'e2spt_refprep.py --refstack=' + subpath + '/' + subtomos.replace('.hdf','_ptcls.hdf') + ' --output=' + prepRefFromOrigModel + ' --stack2process=' + options.model
				
				refprepcmd = 'e2spt_refprep.py --stack2match=' + subpath + '/simptcls.hdf --output=' + prepRefFromOrigModel + ' --stack2process=' + options.model
				
				runcmd( refprepcmd )
				ref = prepRefFromOrigModel
				
				#print "\n\n\nRRRRRRRRRRRRRRRRRRRRRRRR\n\n\n REFPREP command was", refprepcmd
				#print "\nRRRRRRRRRRRRRRRRRRRRRRRRRRRR\n\n"
			
			#print "\n\n\n\n\n\n\n\n\n\n\n rrrrrrrrrrrrrrrrrrrrrr THERE WAS SIMMREF, therefore, ref is", ref			
			if samestackformany > 0:
				ref = thesemodels[d].split('/')[-1]
	
				if itersMode:
					#ref = subpath + '/' + thesemodels[d].split('/')[-1].replace('.hdf','_previousITER.hdf')
					ref = thesemodels[d].split('/')[-1].replace('.hdf','_previousITER.hdf')

				#print "BUT THERE is ALSO SAMESTACK FOR MANY!, so ref is", ref
				#print "rrrrrrrrrrrrrrrr\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
				
			output = subtomos.replace('.hdf', '_avg.hdf')
			#print "\n\n$$$$$$$$$$$$$$$$$$$$$$\nRef name is\n$$$$$$$$$$$$$$$$$$$\n", ref

			#print "\n\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%noutput name is\n", output

			alipath1 = output.split('/')[-1].replace('_avg.hdf','_ali')
			#alipath2= subpath + '/' + output.replace('_avg.hdf','_ali')
			alipath2 = subpath + '/' + alipath1
			#print "\n##################################\nAlipath1 for results will be\n", alipath1
			#print "\n##################################\nAlipath2 for results will be\n", alipath2
			
			subtomos = subpath + '/' + subtomos
			if rootpath not in subtomos:
				subtomos = rootpath + '/' + subtomos
			
			#>>alicmd = " && e2spt_classaverage.py --path=" + alipath1 + " --input=" + subtomos.replace('.hdf','_ptcls.hdf') + " --output=" + output + " --ref=" + ref + " --radius=" + str(options.radius) + " --mask=mask.sharp:outer_radius=-2 --lowpass=" + options.lowpass + " --highpass=" + options.highpass + " --align=rotate_translate_3d:search=" + str(options.transrange) + ":delta=12:dphi=12:verbose=0 --parallel=" + options.parallel + " --falign=refine_3d_grid:delta=3:range=12:search=2 --averager=mean.tomo --aligncmp=" + options.aligncmp + " --faligncmp=" + options.faligncmp + " --shrink=" + str(options.shrinkalign) + " --shrinkfine=" + str(options.shrinkalign) +" --savesteps --saveali --normproc=" + str(options.normproc)  + ' --iter=' + str(options.iter) + ' --npeakstorefine=' + str(options.npeakstorefine) + ' --verbose=' + str(options.verbose)
			
			alicmd = " && e2spt_classaverage.py --path=" + alipath1 + " --input=" + subtomos + " --output=" + output + " --ref=" + ref + " --radius=" + str(options.radius) + " --precision=" + str(options.precision) + " --mask=mask.sharp:outer_radius=-2 --lowpass=" + options.lowpass + " --highpass=" + options.highpass + " --align=rotate_translate_3d:search=" + str(options.search) + ":delta=12:dphi=12:verbose=0 --parallel=" + options.parallel + " --falign=refine_3d_grid:delta=3:range=12:search=2 --averager=mean.tomo --aligncmp=" + options.aligncmp + " --faligncmp=" + options.faligncmp + " --shrink=" + str(options.shrinkalign) + " --shrinkfine=" + str(options.shrinkalign) +" --savesteps --saveali --normproc=" + str(options.normproc)  + ' --iter=' + str(options.iter) + ' --npeakstorefine=' + str(options.npeakstorefine) + ' --verbose=' + str(options.verbose)
			#alicmd += " --nocenterofmass"
			
			
			if options.shrinkalign or options.lowpass or options.highpass:
				alicmd += ' --refpreprocess'
			
			if options.quicktest:
				#>>alicmd = " && e2spt_classaverage.py --path=" + alipath1 + " --input=" + subtomos.replace('.hdf','_ptcls.hdf') + " --output=" + output + " --ref=" + ref + " --radius=" + str(options.radius) + " -v 0 --mask=mask.sharp:outer_radius=-2 --lowpass=" + options.lowpass + " --highpass=" + options.highpass + " --align=rotate_symmetry_3d:sym=c1:verbose=0 --parallel=" + options.parallel + " --falign=None --averager=mean.tomo --aligncmp=" + options.aligncmp + " --faligncmp=" + options.faligncmp + " --shrink=3 --normproc=" + str(options.normproc) + ' --iter=' + str(options.iter) + ' --npeakstorefine=1 --savesteps --saveali --refpreprocess'
				alicmd = " && e2spt_classaverage.py --path=" + alipath1 + " --input=" + subtomos + " --output=" + output + " --ref=" + ref + " --radius=" + str(options.radius) + " -v 0 --mask=mask.sharp:outer_radius=-2 --lowpass=" + options.lowpass + " --highpass=" + options.highpass + " --align=rotate_symmetry_3d:sym=c1:verbose=0 --parallel=" + options.parallel + " --falign=None --averager=mean.tomo --aligncmp=" + options.aligncmp + " --faligncmp=" + options.faligncmp + " --shrink=3 --normproc=" + str(options.normproc) + ' --iter=' + str(options.iter) + ' --npeakstorefine=1 --savesteps --saveali --refpreprocess'
				#alicmd += " --nocenterofmass"
			
			#You want --wedgeangle (or data collection range) to be a bit bigger than it actually is
			#so that the missing wedge is smaller than it actually is, to ensure there will be no overlap with any data.
			#This is effectively accomplished by adding 1 (or more) to 'tiltrange'
			
			#if 'fsc.tomo' in options.aligncmp or 'fsc.tomo' in options.faligncmp:
			#	print "\n\n\n\n\n\n\n\n$$$$$$$$$$$$$$$$\nYOU are selecting FSC.TOMO, therefore, wedgeangle needs to be specified", tiltrange+1
			#	alicmd += ' --wedgeangle=' + str(tiltrange+1)  + ' --wedgei=' + str(options.wedgei) + ' --wedgef=' + str(options.wedgef)
				
			#if options.fitwedgepost:
			#	print "\n\nYou have selected --fitwedgwepost"
			#	alicmd += ' --fitwedgepost'
			
			if options.writewedge:
				alicmd += ' --writewedge'
			
			
			cmd2 = cmd2 + alicmd
			runcmd( cmd2 )	
			
			print "\nCmd2+alicmd is", cmd2
			
			aliptcls = output.replace('_avg.hdf','_ali.hdf')
			
			finalAvgPath = alipath2 + '/' + output
			
			
			#print "\n\aliptcls name is\n", aliptcls

			extractcmd0 = "cd " + alipath2 + " && e2proc3d.py class_0_ptcl.hdf " + aliptcls
			runcmd( extractcmd0 )
			#finalRefVsOriginalRefPath = ''
			
			#originalref = ref
			#refprepcmd=''
					
			if int(options.iter) > 1:
				finalrefindx = int(options.iter) -2
				extractcmd1 = 'cd ' + alipath2 +' && e2proc3d.py class_0.hdf finalref.hdf --first=' + str(finalrefindx) + ' --last=' + str(finalrefindx)
				
				runcmd( extractcmd1 )
				print "\n\n\n EXTRACT command 1 was", extractcmd1
				print "\n\n\n"
				
				finalrefname = alipath2 + '/finalref.hdf'
				
				
				prepRefFromOrigModel = ''
				refToCompensate = ''
				
				'''
				#This makes NO SENSE AT ALL. With iter > 1, the reference should be the previous average. PERIOD.
				
				if not options.simref:
					prepRefFromOrigModel = subpath + '/' + options.model.split('/')[-1].replace('.hdf','_REFPREP.hdf')
					refprepcmd = 'e2spt_refprep.py --refstack=' + finalrefname + ' --output=' + prepRefFromOrigModel + ' --stack2process=' + options.model
					
					refToCompensate = prepRefFromOrigModel
					runcmd( refprepcmd )
					
					print "\n\n\n REFPREP command was", refprepcmd
					print "\n\n\n"
				
					
				else:
				'''
				
				refToCompensate = subpath + '/' + ref
					
					
				
				if not options.quicktest:
					extractcmd2 = "cd " + alipath2 + " && e2spt_classaverage.py --path=finalRef_vs_OriginalRef --input=" + finalrefname + " --output=finalrefali.hdf --ref=" + refToCompensate + " --mask=mask.sharp:outer_radius=-2 --lowpass=" + options.lowpass + " --highpass=" + options.highpass + " --align=rotate_translate_3d:search=" + str(options.search) + ":delta=12:dphi=12:verbose=0 --parallel=" + options.parallel + " --falign=refine_3d_grid:delta=3:range=12:search=2 --averager=mean.tomo --aligncmp=" + options.aligncmp + " --faligncmp=" + options.faligncmp + " --shrink=" + str(options.shrinkalign) + " --shrinkfine=" + str(options.shrinkalign) +" --savesteps --saveali --normproc=" + str(options.normproc)  + ' --iter=1 --npeakstorefine=' + str(options.npeakstorefine) + ' --verbose=' + str(options.verbose)
				
					if options.shrinkalign or options.lowpass or options.highpass:
						extractcmd2 += ' --refpreprocess'
				
				elif options.quicktest:
					extractcmd2 = "cd " + alipath2 + " && e2spt_classaverage.py --clipali=32 --path=finalRef_vs_OriginalRef --input=" + finalrefname + " --output=finalrefali.hdf --ref=" + refToCompensate + " -v 0 --mask=mask.sharp:outer_radius=-2 --lowpass=" + options.lowpass + " --highpass=" + options.highpass + " --align=rotate_symmetry_3d:sym=c1:verbose=0 --parallel=" + options.parallel + " --falign=None --averager=mean.tomo --aligncmp=" + options.aligncmp + " --faligncmp=" + options.faligncmp + " --shrink=3 --normproc=" + str(options.normproc) + " --iter=1 --npeakstorefine=1 --savesteps --saveali --refpreprocess"
				
				runcmd( extractcmd2 )
			
				print "\n\n\n EXTRACT command 2 was", extractcmd2
				print "\n\n\n"
				
			
			#>>resultsfile = aliptcls.replace('_ali.hdf','_ali_error.txt')
			
			resultsfile = subpath.split('/')[-1] + '_ali_error.txt'
			
			#print "\n@@@@@@@\n@@@@@@@\n@@@@@@@@@\n@@@@@@@ Results file is %s \n@@@@@@@\n@@@@@@@\n@@@@@@@@@\n@@@@@@@" %(resultsfile)


			#resultsdir = ''
			#if resultsdir.split('/')[-1] not in subdirs and 'results' not in subdirs:
			#	os.system('mkdir ' + resultsdir)
			
			solutioncmd = "e2spt_transformdistance.py --nolog --input=" + alipath2 + '/' + aliptcls + ' --output=' + resultsdir + '/' + resultsfile
			
			#if options.shrinkalign:
			#	solutioncmd += ' --shrink=' + str( options.shrinkalign )
			
			
			if int(options.iter) > 1:
				finalRefVsOriginalRefPath = alipath2 + '/finalRef_vs_OriginalRef'
				#finalrefaliptcls = finalRefVsOriginalRefPath + '/finalref_ali.hdf'
				
				solutioncmd += ' --transform=' + finalRefVsOriginalRefPath + '/tomo_xforms.json'
			
			runcmd( solutioncmd )
			
			print "\n\n\n SOLUTION command was", solutioncmd
			print "\n\n\n"
				
			
			#rfilecmd =  'mv ' + resultsfile + ' ' + resultsdir
			#if rootpath not in rfilecmd:
			#	rfilecmd =  'mv ' + resultsfile + ' ' +  rootpath + '/' + resultsdir

			#runcmd( rfilecmd )
			
			#if 'mpi' in options.parallel:
			#	genpbs(cmd,kk)
				#a=open('temp.pbs',r)
				#a.write(cmd)
				#a.close()
			#kk+=1
		
		#cmds = cmd2 + extractcmd + solutioncmd + rfilecmd
		
		if options.verbose:
			if cmd1a:
				print "\n\ncmd1a (jobcmd) was:\n", cmd1a
			elif cmd1b:
				print "\n\ncmd1b (to copy randstack before simulation) was:\n", cmd1b
			
			print "\n\ncmd2 (alicmd) was:\n", cmd2
			
			if not options.simref and int(options.iter) >1:
				print "\n\nrefprepcmd was:\n",refprepcmd
			
			print "\n\nextractcmd0 was:\n", extractcmd0
			
			if extractcmd1:
				print "\n\nextractcmd1 was:\n",extractcmd1
			if extractcmd2:
				print "\n\nextractcmd2 was:\n",extractcmd2
			
			print "\n\nsolutioncmd was:\n", solutioncmd

			#print "rfilecmd was:", rfilecmd
		#print "\n\n\n*********************The command to execute is \n %s \n*********************\n" %(cmd)

		#if 'mpi' in options.parallel:
			
		#	os.system('qsub temp'+str(kk)+'.pbs')
		#else:
		
		#if options.verbose:
		#	print "The command is", cmd
		#os.system(cmd)

		if samestackformany == 0:
			#>>thisstack = subpath + '/' + thisstack.split('/')[-1].replace('.hdf','_ptcls.hdf')
			
			print "\nthe stack before messing up is", thisstack
			thisstack = subpath + '/' + thisstack.split('/')[-1]
			
			print "\nThis is where I'm messing up the stack see", thisstack
			
			thesestacks.append(thisstack)			
			
			if not itersMode:
				if options.simref:			
					#>>thissimmodel = subpath + '/' + thissimmodel.split('/')[-1].replace('.hdf','_sptsimMODEL_SIM.hdf')
					
					print "\n and the model before", thissimmodel
					thissimmodel = subpath + '/' + thissimmodel.split('/')[-1]
					
					thissimmodel = subpath + '/simmodel.hdf'
					
					print "\n and the model after", thissimmodel
					thesemodels.append(thissimmodel)
				
				
				else:
					thissmodel = prepRefFromOrigModel = subpath + '/' + options.model.split('/')[-1].replace('.hdf','_REFPREP.hdf')
					thesemodels.append(thissmodel)
					
	
		if itersMode and finalAvgPath:
			finalAvgs.append(finalAvgPath)
			
	if itersMode:
		thesemodels = finalAvgs
	
	#if samestackformany == 0:
	return [thesestacks,thesemodels,firstrandstack]
	#else:
	#	return ['',thesemodels,'']








			
"""
After running simulation and alignment commands (if --testalignments is provided, 
the results need to be compiled, parsed/analyzed and plotted.
"""
	
def analyzeresults(options,rootpath,nrefs,iterPATHS,extra=''):
	
	
	resultsdir=''	
	
	for i in range(nrefs):
		#print "\n\n\n\n I will analyze alignments results for ref number", i
	
		comps = options.comparators.split(',')
		#iters = str(options.iter).split(',')
		
		print "\n\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCn\nIn analyze results, comaprators are", comps

		if options.comparators and len(comps) >1 :
			originalpath = rootpath
			compFilts = []
			for comp in comps:
			
				compID = comp.split(':')[0].replace('.','P') 
				if len(options.path.split('_')) > 1:
					compID += '_' + options.path.split('_')[-1]
			
				#if 'fsc.tomo' in comp:
				#	compID = comp.split(':')[0].replace('.','P') +  comp.split(':')[1].replace('sigmas=','Sigmas').replace('.','P')
				
				compFilt = originalpath + compID
				compFilts.append( compFilt )
		
			compFilts = set(compFilts)
		
		
			for filt in compFilts:
				rootpath = filt
				print "\nComparator filt and refnum are", filt,i
				print "Whereas originalpath, refnum are", originalpath,i
				print "But CURRENT path is", options.path
				print "And ROOTPATH is", rootpath
			
				otherFilts = compFilts - set(filt)
				#options.aligncmp = comp
				#options.faligncmp = comp
			

				#print "\nTherefore compID is", compID
				#print "\n"
			
				potentialdirs=[]
				findir = os.listdir(rootpath)

				#findir = os.listdir(filt)
				#filterID = originalpath + compID
			
				for f in findir:
					if filt == f.split('_')[0]:
						potentialdirs.append(f)
			
				thisone = filt
				if len(potentialdirs) > 1:
					potentialdirs.sort()
					thisone = potentialdirs[-1]	
			
			
				#print "In simloop after alignments, to ANALYZE results, the RESULTS dir WITHOUT hyphen should be based on filt",filt
				#print "based on the originalpath", originalpath
				#print "or now rather on the thisone path", thisone
			
				options.path = thisone
			
				if nrefs > 1:
					modname = 'model' + str(i).zfill(len(str(nrefs)))
		
					resultsdir = modname + '/results_' + modname
					if rootpath not in resultsdir:
						resultsdir = rootpath + '/' + resultsdir
	
				else:
					#resultsdir = options.path + '/results' 
					#if rootpath not in resultsdir:
					#	resultsdir = rootpath + '/' + resultsdir 
					resultsdir = rootpath + '/results'
	
				resfiles = []
			
				try:
					os.system('mkdir ' + resultsdir)
				except:
					pass
								
				findir = os.listdir(resultsdir)
				for f in findir:
					if 'error.txt' in f:
						resfiles.append(f)
	
				resfiles.sort()
	
				resfiles_analysis(options,resfiles,resultsdir,extra,modelnum=i)
	
		elif len(iterPATHS) > 1:
			#print "\n\n\nThere are these many iterPAHTS", len(iterPATHS)
			#print "Which are", iterPATHS
		
			for ele in iterPATHS:
				#print "Analyzing THIS element of iterPATHS", ele
				options.path = iterPATHS[ele]
			
				if nrefs > 1:
					modname = 'model' + str(i).zfill(len(str(nrefs)))
		
					resultsdir = modname + '/results_' + modname
					if rootpath not in resultsdir:
						resultsdir = rootpath + '/' + resultsdir
	
				else:
					resultsdir = options.path + '/results' 
					if rootpath not in resultsdir:
						resultsdir = rootpath + '/' + resultsdir 
	
				resfiles = []
			
				findir = os.listdir(resultsdir)
				for f in findir:
					if 'error.txt' in f:
						resfiles.append(f)
	
				resfiles.sort()
	
				resfiles_analysis(options,resfiles,resultsdir,extra,modelnum=i)			
	
		else:
			resultsdir='results'
			if nrefs > 1:
				modname = 'model' +  str(i).zfill(len(str(nrefs)))
			
				resultsdir = modname + '/results_' + modname
				if rootpath not in resultsdir:
					resultsdir = rootpath + '/' + resultsdir
	
			else:
				#print "(line 899) Results dir is", resultsdir
				#print "options path is", options.path 
				#if options.path not in resultsdir:
				#	resultsdir = options.path + '/results'
				#	print "options.path was not in results dir, therefore latter now is", resultsdir
				#print "rootpath is", rootpath
	
				if rootpath not in resultsdir:
					resultsdir = rootpath + '/' + resultsdir
					print "Rootpath was not in resultsdir, so latter now is", resultsdir
	
				if 'test' in extra:
					resultsdir = rootpath + '/results'
				
				
			resfiles = []
			findir = os.listdir(resultsdir)
			for f in findir:
				if 'error.txt' in f:
					resfiles.append(f)
	
			if len(resfiles) < 2:
				print "Some of your jobs failed. There seems to be only ONE results file, and thus no variability of either SNR, TR, or TS"
				sys.exit()
			else:
				resfiles.sort()		
				resfiles_analysis(options,resfiles,resultsdir,extra,modelnum=i)

	return

	
'''
FUNCTION TO PROPERLY RUN A COMMAND AT THE COMMAND LINE FROM A SCRIPT
'''
def runcmd(cmd):
	
	if cmd:
		p=subprocess.Popen( cmd , shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
		#for line in iter(p.stdout.readline, ''):
		#	print line.replace('\n','')
		
		while True:
			out = p.stdout.read(1)
			if out == '' and p.poll() != None:
				break
			if out != '':
				sys.stdout.write(out)
				sys.stdout.flush()
		
		#p.communicate()	
		p.stdout.close()
	else:
		print "Command to run was empty!"
	return



def resfiles_analysis(options,resfiles,resultsdir,extra,modelnum=0):
	print "Inside resfiles_analysis"
	
	ang_errors = []
	trans_errors = []
	var=''
	
	print "\nExtra received is", extra
	if 'test' in extra:
		var = extra.replace('test','').upper()
		
		print "Therefore var is",var
	
		vals=[]
		resfiles.sort()
		
		for ff in resfiles:
			val=0.0
			#var=''
			xaxis=''
			print "\nParsing this resfile", ff
			ffbase = ff.split('/')[-1].split('_')[0]
			if 'bx' in ffbase:
				val = float(ffbase.replace('bx',''))
				#var = 'BX'
				xaxis = 'Box size (pixels)'
			elif 'sh' in ffbase:
				val = float(ffbase.replace('sh',''))
				#var = 'SH'
				xaxis = 'Shrink factor'
			elif 'msk' in ffbase:
				val = float(ffbase.replace('msk',''))
				#var = 'MSK'
				xaxis = "Masking radius (pixels)"
			elif 'lp' in ffbase:
				val = float(ffbase.replace('lp',''))
				#var = 'LP'
				xaxis = "Low pass resolution (angstroms)"
			elif 'hp' in ffbase:
				val = float(ffbase.replace('hp',''))
				#var = 'HP'
				xaxis = "High pass resolution (angstroms)"

			if val:
				vals.append(val)
			
			
			if xaxis and var:
				resultsfilelocation = resultsdir + '/' + ff
				#print "\n\n%%%%%%%%%%%%\nThe results file is here\n%%%%%%%%%%%%%%%%%\n\n", resultsfilelocation
		
				fff = open(resultsfilelocation,'r')
				lines = fff.readlines()
	
				ang = float( lines[0].split('=')[-1].replace('\n','') )
				ang_errors.append(ang)
	
				trans = float( lines[1].split('=')[-1].replace('\n','') )
				trans_errors.append(trans)
				
			else:
				print "ERROR: There seems to be something wrong with your resfiles"
				print "xaxis is", xaxis
				print "var is", var
		
		if vals and var and ang_errors and trans_errors:
			angfilename = resultsdir+'/' + extra + '_angular_error_var' + var + '.txt'
			oneD_plot( vals , ang_errors, angfilename.replace('.txt','.png'),xaxis,'angular error')
			writeresultsfile( vals, ang_errors, angfilename)
		
			transfilename = angfilename.replace('angular','translational')
			oneD_plot(vals,trans_errors,transfilename.replace('.txt','.png'),'tilt step','translational error')
			writeresultsfile(vals,trans_errors,transfilename)
			
	
	
	
	
	elif 'test' not in extra:
		
		
		snrs = []
		trs = []
		tss = []
		
		twoD_tr_ts_points = []
		twoD_snr_tr_points = []	
		twoD_snr_ts_points = []
			
		threeD_points = []

		resfiles.sort()
		print "\n\n\nIN resfiles_analysis, sorted files are", resfiles
	
		for ff in resfiles:
			print "\nParsing this one", ff
			tr = float(ff.split('TR')[-1].split('_')[0])
			trs.append(tr)

			snr = float( ff.split('SNR')[-1].split('_')[0].replace('p','.') )
			snrs.append(snr)
		
			ts = float(ff.split('NS')[-1].split('_')[0])
			if options.tiltstep:
				ts = float( options.tiltstep )
				#tss.append(options.tiltstep)
			tss.append(ts)
		
			#3dpoints.append([tr,snr,ts])
	
		
			resultsfilelocation = resultsdir + '/' + ff
			#print "\n\n%%%%%%%%%%%%\nThe results file is here\n%%%%%%%%%%%%%%%%%\n\n", resultsfilelocation
		
			fff = open(resultsfilelocation,'r')
			lines = fff.readlines()
		
			ang = float( lines[0].split('=')[-1].replace('\n','') )
			ang_errors.append(ang)
		
			trans = float( lines[1].split('=')[-1].replace('\n','') )
			trans_errors.append(trans)
		
			threeD_points.append({'tilt range':tr,'tilt step':ts,'noise level':snr,'angular_error':ang,'translational_error':trans})
			twoD_tr_ts_points.append({'tilt range':tr,'tilt step':ts,'angular_error':ang,'translational_error':trans})
			twoD_snr_tr_points.append({'tilt range':tr,'noise level':snr,'angular_error':ang,'translational_error':trans})
			twoD_snr_ts_points.append({'tilt step':ts,'noise level':snr,'angular_error':ang,'translational_error':trans})
	
		print "The len(set(snrs)) is", len(set(snrs))
	
		print "The len(set(trs)) is", len(set(trs))
			
		print "The len(set(tss)) is", len(set(tss))
	
		if len(set(snrs)) == 1 and len(set(tss)) == 1 and len(set(trs)) == 1:
			print "Some of your jobs failed. There seems to be only ONE results file, and thus no variability of either SNR, TR, or TS"
			sys.exit()
		
	
		if len(set(snrs)) == 1: 
			if len(set(trs)) == 1 and len(set(tss)) > 1:
				angfilename = resultsdir+'/' + resultsdir.split('/')[-2] + '_angular_error_varNS.txt'
				oneD_plot(tss,ang_errors,angfilename.replace('.txt','.png'),'tilt step','angular error')
				writeresultsfile(tss,ang_errors,angfilename)
			
				transfilename = angfilename.replace('angular','translational')
				oneD_plot(tss,trans_errors,transfilename.replace('.txt','.png'),'tilt step','translational error')
				writeresultsfile(tss,trans_errors,transfilename)
			
			if len(set(tss)) == 1 and len(set(trs)) > 1:
				angfilename = resultsdir+'/' + resultsdir.split('/')[-2] + '_angular_error_varTR.txt'
				oneD_plot(trs,ang_errors,angfilename.replace('.txt','.png'),'tilt range','angular error')
				writeresultsfile(trs,ang_errors,angfilename)
			
				transfilename = angfilename.replace('angular','translational')
				oneD_plot(trs,trans_errors,transfilename.replace('.txt','.png'),'tilt range','translational error')
				writeresultsfile(trs,trans_errors,transfilename)
	
		if len(set(trs)) == 1: 
			if len(set(tss)) == 1 and len(set(snrs)) > 1:
				angfilename = resultsdir+'/' + resultsdir.split('/')[-2] + '_angular_error_varSNR.txt'
				oneD_plot(snrs,ang_errors,angfilename.replace('.txt','.png'),'noise level','angular error')
				writeresultsfile(snrs,ang_errors,angfilename)
			
				transfilename = angfilename.replace('angular','translational')
				oneD_plot(snrs,trans_errors,transfilename.replace('.txt','.png'),'noise level','translational error')
				writeresultsfile(snrs,trans_errors,transfilename)
			
		if len(set(snrs)) == 1 and len(set(trs)) > 1 and len(set(tss)) > 1:
			twoD_plot(twoD_tr_ts_points,val1='tilt range',val2='tilt step',location=resultsdir +'/')
	
		if len(set(trs)) == 1 and len(set(snrs)) > 1 and len(set(tss)) > 1:
			twoD_plot(twoD_snr_ts_points,val1='noise level',val2='tilt step',location=resultsdir +'/')
	
		if len(set(tss)) == 1 and len(set(trs)) > 1 and len(set(snrs)) > 1:
			twoD_plot(twoD_snr_tr_points,val1='noise level',val2='tilt range',location=resultsdir +'/')
	
		if len(set(tss)) > 1 and len(set(trs)) > 1 and len(set(snrs)) > 1:
			threeD_plot(threeD_points,resultsdir +'/')
	
	print "Done with resfiles_analysis"
	return()




#def genpbs(cmd,kk):
#	lines = """#!/bin/bash
#
# This is an example PBS/Torque script for testing your EMAN2 MPI installation. Clearly you should
# modify the number of nodes and ppn (processors per node) before running it on a cluster.
#
#PBS -l nodes=21:ppn=8
#PBS -l walltime=47:59:00

#cd $PBS_O_WORKDIR

#cat $PBS_NODEFILE 

#mpirun /home2/jgalaz/EMAN2/Python/bin/python mpi_test_basic.py > pbs.test.results
#mpirun /home2/jgalaz/EMAN2/Python/bin/python mpi_test.py > pbs.test.results

# NOTE : if you want to run an EMAN2 refinement or somesuch using PBS, you do NOT call mpirun directly
# as in the above statment. Instead, you use the native --parallel option, and EMAN2 will call
# mpirun itself. The number of processors you specify to the --parallel option MUST match the
# number of processors (nodes * ppn) you specify above in the PBS lines. Note the --parallel option
# in the e2refine.py command below. This is otherwise a normal EMAN2 command.

# This is just a commented out example. NOTE: the 'e2bdb.py -c' command at the end of the script. THIS IS CRITICAL !
# You MUST also run 'e2bdb.py -c' manually on the head-node before issuing the 'qsub' command, and should not work
# on any files in the project directory from the head node while the job is running. See the Wiki for more."""
#	lines+= '\n\ne2bdb.py - c &&' + cmd + '\n'
#	
#	a=open('temp'+str(kk)+'.pbs','w')
#	a.write(lines)
#	a.close()
#
#	return



def color(value):
	color =	colorsys.hsv_to_rgb( float(value) / 180.0 / (1.1), 1, 1)
	return(color)


def writeresultsfile(x,y,filename):
	print "\n\n\n"
	print "In writeresultsfile, I am going to write",filename

	aa = open(filename,'w')
	lines=[]
	for i in range(len(x)):
		line=str(x[i]) + ' ' + str(y[i]) + '\n'
		lines.append(line)
	aa.writelines(lines)
	aa.close()
	print "DONE.\n\n\n"
		
	return()


def oneD_plot(points,errors,name,xlabel,ylabel):
	print "INSIDE oneD plotter"
	titlee=name.split('/')[-1].replace('.png','').replace('_',' ')
	#plt.title(title)
	
	fig = plt.figure()
	fig.suptitle(titlee)
	
	plt.xlabel(xlabel,fontweight='bold')
	plt.ylabel(ylabel,fontweight='bold')
	#plt.xlim([min(points)-min(points)*0.1,max(points)+max(points)*0.1])
	
	plt.xlim( [0,max(points)+max(points)*0.1] )
	plt.ylim( [0,max(errors)+max(errors)*0.1] )
	plt.tick_params(axis='both', which='major', labelsize=16)
	plt.tick_params(axis='both', which='minor', labelsize=12)
	
	plt.plot(points,errors,color='k',marker='x',linewidth=2)
	plt.savefig(name,bbox_inches=0)
	print "Will save plot here", name
	plt.clf()
	
	return()
	

def twoD_plot(points,val1,val2,location):
	
	finalpoints = []
	ang_errors = []
	trans_errors = [] 
	x=[]
	y=[]
	for p in points:
		finalpoints.append([ p[val1],p[val2] ])
		x.append(p[val1])
		y.append(p[val2])
		ang_errors.append(p['angular_error'])
		trans_errors.append(p['translational_error'])
	
	plotname1 = location + 'angular_errors_2d_' + '_'.join(val1.split(' ')) + '.png'
	print "\n\n########\nI will save the plot inside 2d_plot to\n########\n\n", plotname1

	plt.title("Angular error")
	plt.xlabel(val1)
	plt.ylabel(val2)
	plt.xlim([0,max(x)+max(x)*0.1])
	plt.ylim([min(y)-min(y)*0.1,max(y)+max(y)*0.1])
	for i in range(len(finalpoints)):
		plt.plot(*zip(*[finalpoints[i]]),marker='o',markersize=4,color=color(ang_errors[i]))	
	plt.savefig(plotname1,bbox_inches=0)
	#plt.clf()
	
	'''
	plotname2 = 'translational_errors_2d_' + '_'.join(val1.split(' ')) + '.png'
	plt.title("Translational error")
	plt.xlabel(val1)
	plt.ylabel(val2)
	for i in range(len(finalpoints)):
		plt.plot(*zip(*[finalpoints[i]]),marker='o',markersize=2,color=color(trans_errors[i]))
	plt.savefig(plotname2)
	#plt.clf()
	'''
	return()
	
	
def threeD_plot(points,location):
	
	finalpoints = []
	ang_errors = []
	trans_errors = []
	x=[]
	y=[]
	z=[] 
	for p in points:
		x.append(p['tilt range'])			
		y.append(p['noise level'])
		z.append(p['tilt step'])
	
		finalpoints.append([ p['tilt range'],p['noise level'],p['tilt step'] ])
		ang_errors.append(p['angular_error'])
		trans_errors.append(p['translational_error'])
	
	plotname1 = location + 'angular_errors_3d.png'
	
	print "\n\n#########\nI will save the plot inside 3d_plot to\n###########\n\n", plotname1
	
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	
	plt.title("Angular error")
	xLabel = ax.set_xlabel('tilt range')
	yLabel = ax.set_ylabel('noise level')
	zLabel = ax.set_zlabel('tilt step')
	
	ax.set_xlim3d(min(x)-min(x)*0.1,max(x)+max(x)*0.1) 
	ax.set_ylim3d(min(y)-min(y)*0.1,max(x)+max(y)*0.1)
	ax.set_zlim3d(min(z)-min(z)*0.1,max(z)+max(z)*0.1)                              
	
	ax.xaxis.labelpad = 30
	ax.yaxis.labelpad = 30

	#ax.dist = 15
	for i in range(len(finalpoints)):
		ax.plot(*zip(*[finalpoints[i]]),marker='o',markersize=4, color=color(ang_errors[i]) )
	plt.savefig(plotname1,bbox_inches=0)
	print "\n\n**********\nI HAVE saved the plot inside 3d_plot to\n********\n\n", plotname1

	#plt.clf()
	
	
	'''
	plotname2 = 'translational_errors_3d.png'
	plt.title("Angular error")
	
	xLabel = ax.set_xlabel('tilt range')
	yLabel = ax.set_ylabel('noise level')
	zLabel = ax.set_zlabel('tilt step')
	
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	for i in range(len(finalpoints)):
		ax.plot(*zip(*[finalpoints[i]]),marker='o',color=color(trans_errors[i]) )
	plt.savefig(plotname2)
	#plt.clf()
	'''
	return()	
	
	
if __name__ == "__main__":
    main()
