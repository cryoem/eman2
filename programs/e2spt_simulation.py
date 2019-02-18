#!/usr/bin/env python
'''
====================
Author: Jesus Galaz - 2011, Last update: 07/Nov/2017
====================

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
'''
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
from EMAN2 import *
from EMAN2_utils import *
import shutil
from EMAN2jsondb import JSTask,jsonclasses

import sys
import numpy
import math
import random

def main():

	usage = """e2spt_simulation.py <options>.
	 This program DEPENDS on e2spt_classaverage.py from which it imports the sptmakepath function.
	 The options should be supplied in "--option=value", replacing "option" for a valid option name, and "value" for an acceptable value for that option. This program produces simulated sub volumes in random orientations from a given PDB or EM file. The output is ALWAYS in HDF format, since it's the only format supported by E2SPT programs.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--clip", type=int,default=0,help="""The final box size to clip the output subtomograms to.""")								
	parser.add_argument("--gridholesize", type=float,default=1.0,help="""Default=1.0. Size of the carbon hole in micrometers for the simulated grid (this will determine the shifts in defocus for each particle at each tilt angle, depending on the position of the particle respect to the tilt axis; the tilt axis by convention goes parallel to Y through the middle of the tomogram.""")
	parser.add_argument("--icethickness", type=float,default=0.4,help="""Thickness of the specimen to simulate, in microns. Default=0.4; --icethickness will be used to calculate the size of the tomogram in Z in PIXELS for the simulated tomogram. This parameter will also be used to assign a random coordinate in Z to each subtomogram.""")
	parser.add_argument("--input", type=str, default='', help="""The name of the input volume from which simulated subtomograms will be generated. The output will be in HDF format, since volume stack support is required. The input CAN be PDB, MRC or and HDF stack. If the input file is PDB or MRC, a version of the supplied model will be written out in HDF format. If the input file is a stack, simulatd subvolumes will be generated from each model in the stack and written to different output stacks. For example, if the input file contains models A and B, two output stacks with simulated subvolumes will be generated.""")
	parser.add_argument("--invert",action="store_true",default=False,help=""""This will multiply the pixel values by -1. This is intended to make the simulated particles be like real EM data before contrast reversal (black, negative contrast), assuming that they're being generated from a model/image where the protein has positive values. It not supplied, 'white protein' (positive density values) will be used by default (or whatever the original contrast is of the image supplied as a model).""")
	
	parser.add_argument("--nosim", action="store_true",default=False,help="""If on, the program will generate stacks of "perfect particles" in different random orientations, but with no missing wedge, no noise, no ctf parameters, etc. The output randstack.hdf will be identical to simptcls.hdf""")	
	parser.add_argument("--notrandomize",action="store_true",default=False,help="This will prevent the simulated particles from being rotated and translated into random orientations.")
	parser.add_argument("--nptcls", type=int,default=10,help="""Number of simulated subtomograms to generate per reference model supplied.""")	
	parser.add_argument("--nslices", type=int,default=61,help="""This will determine the tilt step between slices, depending on tiltrange. For example, to simulate a 2 deg tilt step supply --nslices=61 --tiltrange=60. Recall that --tiltrange goes from - to + the supplied value, and that there is a central slice or projection at 0 deg, for symmetrical tilt series.""")	
	
	parser.add_argument("--path",type=str,default='sptsim',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsim'; for example, sptsim_02 will be the directory by default if 'sptsim_01' already exists.""")
	parser.add_argument("--pad3d", type=float,default=0.0,help="""Factor to calculate the boxsize to use for 3D reconstruction. For example, if the model in --input has an original boxsize with its largest dimension of 64 and you enter --pad3d=1.5x, then the volume used for 3D reconstruction will be 1.5*64, that is, 96. If you provide --shrink, for example, --shrink=2, then the reconstruction box will be 64/2 * 1.5 = 48. Make sure to supply --clip to clip the simulated subtomograms to the final desired box size; otherwise they will be clipped to the current largest dimension of the supplied model/volume.""")								
	parser.add_argument("--pad2d", type=float,default=0.0,help="""Factor to pad projections in the tilt series by before reconstruction.""")								
	parser.add_argument("--parallel",type=str,default='thread:1',help="""Default=thread:1. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel""")
	parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross platform PPID")
	parser.add_argument("--preferredside",action='store_true',default=False,help='''Default=False. If supplied, this option will cause the distribution of orientations to be biased towards alt=90. Works in conjuction with --preferredtop, in which case half of the particles will be biased towards 'top' view orientations and half towards 'side' view orientations.''')
	parser.add_argument("--preferredtop",action='store_true',default=False,help='''Default=False. If supplied, this option will cause the distribution of orientations to be biased towards alt=180 and alt=0. Works in conjuction with --preferredside, in which case half of the particles will be biased towards 'top' view orientations and half towards 'side' view orientations.''')

	parser.add_argument("--randstack",type=str,default='',help="If you already have a stack of particles (presumably in random orientations) you can supply it here.")
	parser.add_argument("--reconstructor", type=str,default="fourier",help="""The reconstructor to use to reconstruct the tilt series into a tomogram. Type 'e2help.py reconstructors' at the command line to see all options and parameters available. To specify the interpolation scheme for the fourier reconstruction, specify 'mode'. Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 'gauss_5_slow', 'gypergeom_5', 'experimental'. For example --reconstructor=fourier:mode=gauss_5 """)																				
	
	parser.add_argument("--savenoise", action="store_true",default=False,help="""If on, it saves the noise stack for each particle. This can be useful for testing alignment under varying SNR, so that the same noise (just at a different ratio/level) is tested.""")
	parser.add_argument("--saveorthostack", action="store_true",default=False,help="If on, --nptcls is ignored and you get 3 subtomograms (simulated from the model supplied) which are orthogonal to each other.")
	parser.add_argument("--saverandstack", action="store_true",default=True,help="""DEPREPCATED. [This option is on by default and there's no way to turn it off. The stack of randomly oriented particles before simulating the missing wedge WILL be saved]. Save the stack of randomly oriented particles, before subtomogram simulation (before the missing wedge and noise are added).""")
	parser.add_argument("--saveprjs", action="store_true",default=False,help="""Save the projections (the 'tilt series') for each simulated subtomogram.""")	
	parser.add_argument("--savetlt",action="store_true",default=False,help="""Save a text file with .tlt extension (as in IMOD) containing the tilt angles for the simulated tomogram and/or subtomograms.""")
	parser.add_argument("--snr",type=float,default=0,help="Weighing noise factor for noise added to the image.")
	parser.add_argument("--sym",type=str,default='c1',help="If your particle is symmetrical, it is only necessary to randomize orientations within the asymmetric unit only.")
	parser.add_argument("--simref",action="store_true",default=False,help="This will make a simulated particle in the same orientation as the original input (or reference).")
	parser.add_argument("--set2tiltaxis",action='store_true',default=False,help="""Default=False. Simulate particles along the tilt axis only.""")
	
	parser.add_argument("--tiltaxis",type=str,default='y',help="""Axis to produce projections about. Default is 'y'; the only other valid option is 'x'.""")
	parser.add_argument('--tiltangles',type=str,default='',help="""File in .tlt or .txt format containing the tilt angle of each tilt image in the tiltseries.""")
	
	parser.add_argument("--txrange", type=int,default=0,help="""Maximum number of pixels to randomly translate each subtomogram in X. The random translation will be picked between -txrange and +txrange. Default value is set by --trange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--txerror", type=int,default=0,help="""Range of random translation error in pixels to perturb individual 2-D images in each subtiltseries by along x. The random translation perturbation will be picked between -txerror and +txerror. Default value is set by --terror, but --txerror will overwrite it if specified.""")
	parser.add_argument("--tyrange", type=int,default=0,help="""Maximum number of pixels to randomly translate each subtomogram in Y. The random translation will be picked between -tyrange and +tyrange. Default value is set by --trange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--tyerror", type=int,default=0,help="""Range of random translation error in pixels to perturb individual 2-D images in each subtiltseries by along y. The random translation perturbation will be picked between -tyerror and +tyerror. Default value is set by --terror, but --tyerror will overwrite it if specified.""")
	parser.add_argument("--tzrange", type=int,default=0,help="""Maximum number of pixels to randomly translate each subtomogram in Z. The random translation will be picked between -tzrange and +tzrange. Default value is set by --trange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--trange", type=int,default=0,help="""Maximum number of pixels to randomly translate each subtomogram in all X, Y and Z. The random translation will be picked between -transrage and +trange; --txrange, --tyrange and --tzrange overwrite --trange for each specified direction.""")
	parser.add_argument("--terror", type=int,default=0,help="""Range of random translation error in pixels to perturb individual 2-D images in each subtiltseries by along x, y and z. The random translation perturbation will be picked between -terror and +terror. If set, this will overwrite --txerror, --tyerror and --tzerror.""")
	parser.add_argument("--tiltrange", type=float,default=60,help="""Maximum angular value at which the highest tilt picture will be simulated. Projections will be simulated from -tiltrange to +titlrange. For example, if simulating a tilt series collected from -60 to 60 degrees, enter a --tiltrange value of 60. Note that this parameter will determine the size of the missing wedge.""")
	
	parser.add_argument("--verbose", "-v", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness", dest="verbose", action="store", metavar="n")

	'''
	CTF PARAMETERS
	'''
	parser.add_argument("--applyctf", action="store_true",default=False,help="If on, it applies ctf to the projections in the simulated tilt series based on defocus, cs, and voltage parameters.")
	parser.add_argument("--defocus", type=float,default=3.0,help="""Target defocus at the tilt axis (in microns) for the simulated tilt series. Default is 3.0. Notice that DEFOCUS (underfocus) values are POSITIVE, by convention.""")
	parser.add_argument("--voltage", type=int,default=200,help="""Voltage of the microscope, used to simulate the ctf added to the subtomograms. Default is 200 KV.""")
	parser.add_argument("--cs", type=float,default=2.1,help="""Cs of the microscope, used to simulate the ctf added to the subtomograms. Default is 2.1.""")
	parser.add_argument("--apix",type=float,default=0.0,help="""Provide accurate apix in case the header has the wrong apix info.""")	
	parser.add_argument("--bfactor",type=int,default=400,help="""Bfactor to use for CTF correction phase flipping. Default is 400.""")
	parser.add_argument("--ampcont",type=float,default=0.05,help="""Amplitude contrast to use for CTF correction phase flipping. Default is 0.05.""")	
	
	#parser.add_argument("--interpolator",default='',help="""What interpolation scheme to use for reconstruction. Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 'gauss_5_slow', 'gypergeom_5', experimental""")	
	#parser.add_argument("--fillwedge",action='store_true',default=False,help="""This option will fill the region of the missing wedge with evenly spaced images of Gaussian noise with the same tiltstep used to simulate the particles (as determined through the parameter --nslices).""")
	
	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)

	#print "e2spt_simulation received these many slices", options.nslices

	'''
	Parse the options
	'''
	#if options.filter:
	#	options.filter = parsemodopt(options.filter)

	#if options.reconstructor:
	#	options.reconstructor = parsemodopt(options.reconstructor)
	
	
	#function imported from EMAN2_utils
	options = sptOptionsParser( options )

	if options.input and options.randstack:
		print("""\n(e2spt_simulation)WARNING: No point in supplying --input and --randstack simultaneously.
		They are mutually exclusive. If --randstack is provided, --input is ignored,
		because --randstack becomes --input.\n
		--input was""",options.input)
		print('--randstack was', options.randstack)
	
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	
	#functions imported from EMAN2_utils
	options = makepath(options,'sptsim')
	writeParameters( options, 'e2spt_simulation.py', 'sptsim' )
	
	rootpath = os.getcwd()
	
	#if rootpath not in options.path:
	#	options.path = rootpath + '/' + options.path
		
	randptcls = {}
	
	
	modelhdr = EMData(options.input,0,True)
	dimension = 3
	
	if int(modelhdr['nz']) < 2:
		dimension = 2
	
	if options.randstack:
		if options.verbose > 3:
			print("\n\nI will not generate a randstack but will read it from", options.randstack)

		randstackbase = os.path.basename( options.randstack )
		randstackcopy = options.path + '/' + randstackbase

		copycmd = 'e2proc3d.py ' + options.randstack + ' ' + randstackcopy
		
		runcmd(options,copycmd)
		
		nr=EMUtil.get_image_count(randstackcopy)
		if options.verbose > 3:
			print("There are these many particles in the randstack", nr)							#ATTENTION: Randstack might still need box size changes and padding...
		
		options.input = randstackcopy
		
		#if options.preprocess or float(options.pad) > 1.0 or int(options.shrink) > 1 or options.lowpass or options.highpass or options.preprocess:
		#	ret = preprocess(options,randstackcopy,dimension) 
		#	randstackcopy = ret[-1]
				
		simptclsname = options.path + '/simptcls.hdf'
		if options.nosim:
			
			if options.clip:
				clip2cmd = 'e2proc3d.py ' + randstackcopy + ' ' + randstackcopy + ' --clip=' + str(options.clip)
				runcmd(options,clip2cmd)
			
			shutil.copy(options.input,simptclsname)
			#os.system('cp ' + options.input + ' ' + simptclsname)
		
			simrefname = options.path + '/simmodel.hdf'
			
			simrefcmd = 'e2proc3d.py ' + randstackcopy + ' ' + simrefname + ' --first=0 --last=0'
	
			runcmd(options,simrefcmd)
			
			
		else:
			for np in range(nr):
				a=EMData(randstackcopy,np)
				#>>randptcls.append(a)
				randptcls.update({np:a})
			
			subtomosim(options, randptcls, simptclsname, dimension )
		
	elif options.input:
	
		fyles = options.input.split(',')
	
		for fyle in fyles:
		
			options.input = fyle
		
			'''
			Convert pdb or mrc files to hdf
			'''
			options = sptfixformat( options )
					
			nrefs = EMUtil.get_image_count( options.input )		
			
			'''
			The program can process several FILES, each of which could be a STACK with several
			images each.
			'''
			#if len(fyles) > 1:
			#	options.output = options.input.replace('.hdf','_simSubtomos.hdf')
		
			tag = ''
	
			if options.verbose:
				print("These many particles will be simulated, for each of the supplied references/models", options.nptcls)
				print("There are these many references/models", nrefs)
		
			originalpath = options.path
			kkk=0
	
			originalinput = options.input
		
			for i in range(nrefs):
				if options.verbose:
					print("\n\nGenerating simulated subtomograms for reference number", kkk)
					
				if nrefs>1:
					modelfilename = originalinput.split('/')[-1].replace('.hdf','_model' + str(i).zfill(2) + '.hdf')

					options.path = originalpath + '/model' + str(i).zfill(2)

					#os.system('mkdir ' + options.path)
					os.mkdir(options.path)
					#cmd = 'e2proc3d.py '  + options.input + ' ' + options.path + '/' + modelfilename + ' --first=' + str(i) + ' --last=' + str(i) + ' --append'
	
					#os.system('e2proc3d.py '  + originalinput + ' ' + options.path + '/' + modelfilename + ' --first=' + str(i) + ' --last=' + str(i) + ' --append')
					
					#print("This is the command to create the model")
					cmdmakemodel = 'e2proc3d.py '  + originalinput + ' ' + options.path + '/' + modelfilename + ' --first=' + str(i) + ' --last=' + str(i) + ' --append'
					runcmd(options,cmdmakemodel)

					options.input = options.path + '/' + modelfilename
					tag = str(i).zfill(len(str(nrefs)))
	
				'''
				Make model's box cubical if it isn't
				'''
				model = EMData( options.input, 0 )
				if dimension == 3:					
					if model['nx'] != model['ny'] or model['nx'] != model['nz'] or model['ny'] != model['nz']:	
						#if options.clip:
						model = clip3d( model, max( model['nx'], model['ny'], model['nz'] ) )
						#modelhdr = EMData(options.input,0,True)
						
				elif model['nx'] != model['ny'] or model['nx'] != model['nz'] or model['ny'] != model['nz']:
					print("\nThe image is 2D")
					model = clip2d( model, max( model['nx'], model['ny'] ) )	
							
				retrand = randomizer(options, model, tag)
				randptcls = retrand[0]
				randstackname = retrand[-1]
				
				simptclsname = options.path + '/simptcls.hdf'
				
				print("\n\n\n\n\n\(e2spt_simulation) before subtomosim, simptclsname is", simptclsname)
				print("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS\n\n\n\n\n\n\n\n")
				if not options.nosim:
					subtomosim(options,randptcls,simptclsname,dimension)		
				else:
					#os.system('cp ' + randstackname + ' ' + simptclsname)
					
					shutil.copy(randstackname,simptclsname)

					simrefname = options.path + '/simmodel.hdf'
			
					simrefcmd = 'e2proc3d.py ' + options.input + ' ' + simrefname + ' --first=0 --last=0'
					
					runcmd(options,simrefcmd)
					
				if options.saveorthostack:
					orthoptcls = orthostack( options, model )
					
					if not options.nosim:
						simorthoptclsname = options.path + '/orthoptcls.hdf'
						subtomosim(options,orthoptcls,simorthoptclsname,dimension)		

				kkk+=1
	
	E2end(logger)		
	return
					

def sptfixformat( options ):
	
	if '.pdb' in options.input[-4:]:
		pdbmodel = options.input
		mrcmodel = options.path + '/' + os.path.basename(options.input).replace('.pdb','.mrc')
		
		cmdpdb2mrc = 'e2pdb2mrc.py ' + pdbmodel + ' ' + mrcmodel
		runcmd(options,cmdpdb2mrc)
		
		hdfmodel = mrcmodel.replace('.mrc','_MODEL.hdf')
		cmdmrc2hdf = 'e2proc3d.py ' + mrcmodel + ' ' + hdfmodel
		
		runcmd(options,cmdmrc2hdf)

		os.remove(mrcmodel)

		options.input = hdfmodel

	elif '.mrc' in options.input[-4:]:
		hdfmodel = options.path + '/' + os.path.basename(options.input).replace('.mrc','_MODEL.hdf')
		cmdmrc2hdf = 'e2proc3d.py ' + mrcmodel + ' ' + hdfmodel
		
		runcmd(options,cmdmrc2hdf)
		options.input = hdfmodel
	
	elif '.mrcs' in options.input[-5:]:
		hdfmodel = options.path + '/' + os.path.basename(options.input).replace('.mrcs','_MODEL.hdf')
		cmdmrc2hdf = 'e2proc3d.py ' + mrcmodel + ' ' + hdfmodel
		
		runcmd(options,cmdmrc2hdf)
		options.input = hdfmodel 		
	
	elif '.hdf' in options.input[-4:]:
		hdfmodel = options.path + '/' + os.path.basename(options.input).replace('.hdf','_MODEL.hdf')

		shutil.copy(options.input,hdfmodel)
		options.input = hdfmodel 

	return options

'''
def clip3D( vol, size ):
	
	volxc = vol['nx']/2
	volyc = vol['ny']/2
	volzc = vol['nz']/2
	
	Rvol =  Region( (2*volxc - size)/2, (2*volyc - size)/2, (2*volzc - size)/2, size , size , size)
	if Rvol:
		vol.clip_inplace( Rvol )
	else:
		print("ERROR!: Empty Region to clip", Rvol)
	#vol.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return vol


def clip2D( img, size ):
	
	imgxc = img['nx']/2
	imgyc = img['ny']/2
	#imgzc = img['nz']/2
	
	Rimg =  Region( (2*imgxc - size)/2, (2*imgyc - size)/2, 0, size , size , 1)
	if Rimg:
		img.clip_inplace( Rimg )
	else:
		print("ERROR!: Empty Region to clip", Rimg)
		
	#img.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return img
'''

'''
====================
RANDOMIZER - Takes a file in .hdf format and generates a set of n particles randomly rotated and translated
====================
'''
def randomizer(options, model, tag):
	
	import random
	
	#print "I am inside the RANDOMIZER"
	
	if options.verbose:
		print("You have requested to generate %d particles with random orientations and translations" %(options.nptcls))
	
	randptcls = {}
	
	#randstackname = options.path + '/' + options.input.replace('.hdf','_randst' + tag + '_n' + str(options.nptcls).zfill(len(str(options.nptcls))) + '.hdf').split('/')[-1]
	#if options.output:
	#	randstackname = options.path + '/' + options.output.replace('.hdf','_randst' + tag + '_n' + str(options.nptcls).zfill(len(str(options.nptcls))) + '.hdf').split('/')[-1]
	
	randstackname = options.path + '/randstack.hdf'
	
	print("###############\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n#################\nThe stackname inside RANDOMIZER, is", randstackname)
	print("--saverandstack is", options.saverandstack)
	print("#####################################\n\n\n\n\n\n\n\n\n\n\n\n\n")
	
	orientations = {}
	randomangles = False	
	randomtrans = False 
	
	if not options.notrandomize:
		
		randomangles =True

		orientations = { 0:Transform() }	#the first particle's orientation is not randomized
		palts=[]

		print("\nGenerating random orientation within asymmetric unit %s" %(options.sym))
		sym = Symmetries.get( options.sym )
		
		#orients = sym.gen_orientations("eman",{"n": options.nptcls,"random_phi":1,"inc_mirror":1})
		
		orients = sym.gen_orientations("rand",{"n": options.nptcls,"phitoo":1,"inc_mirror":1})

		ii=1	#start at index 1 since the first particle's orientation is not randomized
		for o in orients:
			orientations.update( { ii:o } )
			ii+=1

		if options.preferredtop and not options.preferredside:
			#paltstop = preferredalt( options, mu=180,sigma=45, nptcls=options.nptcls )
			ntop = int(round(old_div(options.nptcls,2.0)))
			nbottom = options.nptcls -ntop 
			#palts = paltstop
			palts = numpy.append( preferredalt( options, mu=180,sigma=45, nptcls=ntop ), preferredalt( options, mu=0,sigma=45, nptcls=nbottom ) )
			print("\nreturned palts",palts)
		
		if options.preferredside and not options.preferredtop:
			paltsside = preferredalt( options, mu=90, sigma=22.5, nptcls=options.nptcls ) 
			palts = paltsside

		if options.preferredside and options.preferredtop:
			ntop = int(round(old_div(options.nptcls,4.0)))
			nbottom = ntop
			nside = options.nptcls -ntop -nbottom  	
			palts = numpy.append( preferredalt( options, mu=180,sigma=45, nptcls=ntop ), preferredalt( options, mu=0,sigma=45, nptcls=nbottom ), preferredalt( options, mu=90,sigma=45, nptcls=nside ) )


		palts = numpy.array(palts)

		if options.preferredtop or options.preferredside: #reassign altitude if preferred orientation was specified and new altitudes were generated successfully above

			if palts.any():	

				for i in range( options.nptcls ):

					random_transform = orientations[i]
													
					rots=random_transform.get_rotation()
					
					az=rots['az']
					phi=rots['phi']
					alt=palts[i]
					
					random_transform=Transform({'type':'eman','az':az,'alt':alt,'phi':phi})
					
					orientations.update({ i:random_transform })

	elif options.notrandomize:
		for i in range( options.nptcls ):
			orientations.update({ i:Transform() })

	if options.trange or options.txrange or options.tyrange or options.tzrange :
		
		randomtrans = True
		for i in range( options.nptcls ): 

			randtx = randty = randtz = 0

			if i > 0:	#if translations are also requested to be randomized, do so for all particles after the first
				if options.trange and not options.txrange:
					randtx = random.randrange(-1 * options.trange, options.trange+1)	#Generate random translations
				elif options.txrange:													#add +1 so that python includes the upper limit
					randtx = random.randrange(-1 * options.txrange, options.txrange+1)	
				
				if options.trange and not options.txrange:
					randty = random.randrange(-1 * options.trange, options.trange+1)
				elif options.tyrange:
					randty = random.randrange(-1 * options.tyrange, options.tyrange+1)
				
				if options.trange and not options.tzrange:
					randtz = random.randrange(-1 * options.trange, options.trange+1)
				elif options.tzrange:
					randtz = random.randrange(-1 * options.tzrange, options.tzrange+1)

				if randtx or randty or randtz:
					#random_transform.translate(randtx, randty, randtz)
					print("\nbefore translation orientations[i] is",  orientations[i])
					orientations[i].set_trans(randtx, randty, randtz)
					print("\ntranslated orientations[i] is", orientations[i])
					#orientations.update({ i : newt  })
			
			else:
				pass #first particle never randomized
	
	outtransform = Transform()
	b = model.copy()
	b['sptsim_randT'] = outtransform
	b['xform.align3d'] = Transform()
	
	jsTransformParamsPath = options.path + '/sptsim.json'
	jsA = js_open_dict( jsTransformParamsPath )
	
	jsTransformParamsPath_solution = options.path + '/sptsim_solution.json'
	jsAS = js_open_dict( jsTransformParamsPath_solution )
	
	for i in range( options.nptcls ):
		
		xformslabel = 'subtomo_' + str( i ).zfill( len( str( options.nptcls ) ) )			

		if options.verbose:
			print("\n(e2spt_simulation.py) generating particle #%d" %( i ))
		
		if randomangles or randomtrans:
			
			b = model.copy()

			outtransform = orientations[i]
			jsA.setval( xformslabel, [ outtransform ] )

			outtransformi = outtransform.inverse()
			jsAS.setval( xformslabel, [ outtransformi ] )

			print("\nouttransform", outtransform)
			b.transform(outtransform)
		
			#transforms.append(random_transform)		

			b['sptsim_randT'] = outtransform
			b['xform.align3d'] = Transform()			#This parameter should be set to the identity transform since it can be used later to determine whether
													#alignment programs can successfuly "undo" the random transformation stored in spt_randT
			#b.process_inplace('normalize.edgemean')
			if options.saverandstack:	

				#print "The stackname to use is", stackname
				#randstackname = options.path + '/' + stackname.split('/')[-1]
				#print "\n\n\n\n\nI will save randstack! Using THE PATH in e2spt_simulation and stackname both of which together are", randstackname
				#print "options path received is", options.path
				#print "Whereas stackname is", stackname
				#print "Therefore stackname.split('/')[-1] is", stackname.split('/')[-1]
				
				
				b['origin_x'] = 0									#Make sure the origin is set to zero, to avoid display issues with Chimera
				b['origin_y'] = 0

				#if dimension == 3:
				b['origin_z'] = 0
				
				b.write_image(randstackname,i)
				print("\n(e2spt_simulation.py) saving random orientations stack. particle %d written to %s" % ( i, randstackname ))

			else:
				pass #stack of particles in random orientations not saved
		else:
			pass #angles and translations were not randomized; the stack to be return will contain copies of the model (might be modified by noise/CTF later) 

		randptcls.update({i:b})
		
		if options.verbose:
			print("\n(e2spt_simulation.py) applied transform", random_transform)


	jsA.close()
	jsAS.close()


	azs=[]
	alts=[]
	phis=[]
	xs=[]
	ys=[]
	zs=[]

	#if len(transforms) > 2:
	if len( orientations ) > 2:
		linestrans = []
		if randomangles or randomtrans:
			
			for i in orientations:
				
				t = orientations[i]

				print("\n t to get rotations and translations from is",t)
				if randomangles:
					rots=t.get_rotation()
					
					az=rots['az']
					azs.append(az)

					alt=rots['alt']
					alts.append(alt)
					
					phi=rots['phi']
					phis.append(phi)
					
				if randomtrans:
					trans=t.get_trans()
					
					x=trans[0]
					xs.append(x)

					y=trans[1]
					ys.append(y)

					z=trans[2]
					zs.append(z)

					line2writetrans = str(x) + '\t' + str(y) + '\t' + str(z) + '\n'
					linestrans.append(line2writetrans)

			

		if randomangles:
			textwriter(options, azs,'az')
			plotvals( options, azs, 'az' )

			textwriter(options, alts,'alt')
			plotvals( options, alts, 'alt' )

			textwriter(options, phis,'phi')
			plotvals( options, phis, 'phi' )

		if randomtrans:
			textwriter(options, xs,'x')
			plotvals( options, xs, 'x' )

			textwriter(options, ys,'y')
			plotvals( options,  ys,'y' )

			textwriter(options, zs,'z')
			plotvals( options, zs,'z' )

			if linestrans:
				with open(options.path + '/x_y_z_translations.txt','w') as f:
					print("\n\n\n\n\n\nwwwwwwwwww writing translations text file")
					f.writelines(linestrans)

	return randptcls,randstackname


def preferredalt( options, mu=0, sigma=1, nptcls=3 ):
	s = numpy.random.normal(mu, sigma, nptcls)
	return s


def plotvals( options, vals, tag ):
	import matplotlib.pyplot as plt

	sigmavals= numpy.std(vals)
	meanvals = numpy.mean(vals)

	cuberoot = numpy.power(len(vals),old_div(1.0,3.0))
	width = old_div((3.5*sigmavals),cuberoot)
	
	#print "Therefore, according to Scott's normal reference rule, width = (3.5*std)/cuberoot(n), the width of the histogram bins will be", width
	
	minvals = min(vals)
	maxvals = max(vals)
	
	#if options.bins:
	#	calcbins = options.bins

	if 'az' in tag or 'phi' in tag:
		minvals = 0
		maxvals = 360
	elif 'alt' in tag:
		minvals = 0
		maxvals = 180
	elif 'x' in tag or 'y' in tag or 'z' in tag:
		pass

	calcbins = int(round( old_div((maxvals - minvals ), width) ))

	#count, bins, ignored = plt.hist(vals, 30, normed=True)
	ignored = plt.hist(vals, calcbins)
		
	#plt.plot(bins, 1/(sigmavals * numpy.sqrt(2 * numpy.pi)) * numpy.exp( - (calcbins - meanvals)**2 / (2 * sigmavals**2) ), linewidth=2, color='r')

	plt.title( tag + ' distribution' )
	plt.ylabel("n")
	plt.xlabel(tag)
	
	plt.savefig(options.path + '/' + tag + '.png')
	plt.clf()

	return


def textwriter(options,data,tag):
	
	#if options.path not in name:
	name = options.path + '/' + tag + '.txt'
	
	print("I am in the text writer for this file", name)
	
	lines=[]
	
	for i in range(len(data)):
			
		line2write = str(i) + '\t' + str(data[i]) + '\n'
		#print "THe line to write is"
		lines.append(line2write)
	
	with open(name,'w') as f:
		f.writelines(lines)

	return


def orthostack(options, model):
	orthoptcls={}
	for i in range(4):
		b = model.copy()
		
		if i==0:
			b['sptsim_randT'] = Transform()
			b['xform.align3d'] = Transform()
			orthoptcls.update({2:b})
		
		#if i == 1:
		#	b.rotate(90,0,0)
		#	t = Transform({'type':'eman','az':0,'alt':90,'phi':0})
		#	b['sptsim_randT'] = t
		#	b['xform.align3d'] = Transform()				
		#	randptcls.append(b)
		
		if i == 1:
			b.rotate(0,90,0)
			t = Transform({'type':'eman','az':0,'alt':90,'phi':0})
			b['sptsim_randT'] = t
			b['xform.align3d'] = Transform()				
			orthoptcls.update({1:b})
			
		if i == 2:
			b.rotate(0,90,90)
			t = Transform({'type':'eman','az':0,'alt':90,'phi':90})
			b['sptsim_randT'] = t
			b['xform.align3d'] = Transform()
			orthoptcls.update({0:b})
		
		#orthostackname = options.path + '/orthostack.hdf'
		
	return orthoptcls
	

'''
====================
SUBTOMOSIM takes a set of particles in .hdf format and generates a simulated sub-tomogram for each, using user-defined parameters for tomographic simulation.
The function generates projections for each particle with a user-defined tilt step and missing wedge size (or data collection range),
adds noise and ctf to each projection, randomizes the position of each particle "in the ice" and changes the defocus for each picture in the tilt series accordingly,
and recounstructs a new 3D volume from the simulated tilt series.
====================
'''	
def subtomosim(options,ptcls,outname,dimension):
	#print "INSIDE SUBTOMOSIM"
	
	print("\n\n\n\n\n(e2spt_simulation) Outname received in subtomosim", outname)
	print("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n\n\n\n\n\n")	
	'''
	Initialize parallelism if being used
	'''
	if options.parallel :
		print("\n\n(e2spt_simulation.py) INITIALIZING PARALLELISM, for this outname (stack, or reference)", outname)
		print("\n\n")
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
	
	if options.verbose:
		
		print("(e2spt_simulation) There are these many slices to produce to simulate each subtomogram", options.nslices)
	
	print("\n\n\n\n\n\n\n\n\n\n\nIn subtomosim function before task, size is", ptcls[0]['nx'],ptcls[0]['ny'],ptcls[0]['nz'])
	
	tasks=[]
	#>>for i in range(len(ptcls)):
	
	
	tangles = genangles( options )
	
	for i in ptcls:	
		#if options.parallel:			
		task=SubtomoSimTask(ptcls[i],i,options,outname,tangles)
		tasks.append(task)
	
	tids=etc.send_tasks(tasks)
	
	results = get_results(etc,tids,options)
	
	#pn = 0
	
	#print "I have these many particles from results", len(results)
	#print "From outname", outname
	ii=0
	#>>for pn in range(len(results)):
	finaloutname = outname
	for result in results:
		
		key = list(result.keys())[0]
		
		#if options.path not in outname:
		#	finaloutname = options.path + '/' + outname
		
		#print "\n\n\n\n\n\nn\\n\n\nn\\nn\n\n\n\n\n THe final outname rigth before writing is", finaloutname
		#print "because options.path is", options.path
		#print "and outname is", outname
		#print "And the particle is", results[pn]
		#print "And its type is", type( results[pn])
		#print "And its index is", pn
		
		if i==0:
			print("\n\n(subtomosim) The size of the final particle is",result[key]['nx'],result[key]['ny'],result[key]['nz'])

		result[key]['origin_x'] = 0									#Make sure the origin is set to zero, to avoid display issues with Chimera
		result[key]['origin_y'] = 0
		result[key]['origin_z'] = 0
		
		#if options.tiltaxis == 'y':
		#	tt = Transform({'type':'eman','az':90,'alt':90,'phi':-90})
		#	result[key].transform( tt )
				
		result[key]['xform.align3d'] = Transform()
		#finaloutname = finaloutname.replace('_preproc','')
		result[key].write_image(finaloutname,key)
		#pn+=1
		ii+=1
	
	#if options.clip and dimension == 3:
	#	box = int(options.clip)
	#	print "\nActually, because of clip$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n options.clip is", options.clip
	#	print "Therefore, box of final simulated stack is", box
		
	#	clipcmdf = 'e2proc3d.py ' + finaloutname + ' ' + finaloutname + ' --clip=' + str(options.clip)
	#	print "with cmmand", clipcmdf
	#	
	#	p=subprocess.Popen( clipcmdf, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#	text=p.communicate()	
	#	p.stdout.close()
	
	
	if options.simref and not options.saveorthostack:
		simrefname = options.path + '/simmodel.hdf'
		
		simrefcmd = 'e2proc3d.py ' + finaloutname + ' ' + simrefname + ' --first=0 --last=0'
		p=subprocess.Popen( simrefcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
	
	#if options.writetomogram:
	#	tomogramsim(options,results)

	
	return
	
	
def genangles( options ):

	angles = []
	
	if options.tiltangles:
		anglesfile = open(options.tiltangles,'r')				#Open tilt angles file
		alines = anglesfile.readlines()							#Read its lines
		anglesfile.close()										#Close the file
		
		for line in alines:
			angle = line.replace('\n','').replace(' ','')
			angles.append( float(angle) )
	
	else:
		lower_bound = -1 * options.tiltrange
		alt = lower_bound
		upper_bound = options.tiltrange
		nslices = options.nslices
		tiltstep = round(old_div(( float(upper_bound) - float(lower_bound) ), float(nslices - 1)),2)	
	
		lines = []
	
		for j in range( nslices ):
			realalt = alt + j*tiltstep
			angles.append( realalt )
			if options.savetlt:
				line = str( realalt ) +'\n'
				lines.append( line )
	
		if options.savetlt:
			tiltfile = options.path + '/tiltangles.tlt'
			print("\n\n\n\nPath to save angles in is", tiltfile)
		
			f = open( tiltfile, 'w' )
			f.writelines( lines )
			f.close()
	
	return	angles
	
	
class SubtomoSimTask(JSTask):
	"""This is a task object for the parallelism system. It is responsible for generating a simulated subtomogram."""	
	
	def __init__(self,image,ptclnum,options,outname,tiltangles):
	
		"""fixedimage and image may be actual EMData objects, or ["cache",path,number]
		label is a descriptive string, not actually used in processing
		ptcl is not used in executing the task, but is for reference
		other parameters match command-line options from e2spt_classaverage.py
		Rather than being a string specifying an aligner, 'align' may be passed in as a Transform object, representing a starting orientation for refinement"""
		#data={}
		#data={"fixedimage":fixedimage,"image":image}
		data={"image":image}
		
		JSTask.__init__(self,"SubTomoSim",data,{},"")

		self.classoptions={"options":options,"ptclnum":ptclnum,"outname":outname,"tiltangles":tiltangles}
	
	def execute(self,callback=None):
		"""This simulates a subtomogram and saves projections before and after adding noise if instructed to."""
		classoptions=self.classoptions
		options=self.classoptions['options']
		outname = self.classoptions['outname']
		tiltangles = self.classoptions['tiltangles']
		tiltangles.sort()
		
		i = self.classoptions['ptclnum']
		
		image = self.data['image']
		
		print("\n\n(SubtomoSimTask) Size of the particle for simulation is", image['nx'],image['ny'],image['nz'])
		
		dimension = 3
		if int(image['nz']) < 2:
			dimension = 2
		
		outname = self.classoptions['outname']
		
		if options.verbose:
			print("Generating projections for particle #", i)

		#apix = ptcls[i]['apix_x']
		apix = image['apix_x']
		if options.apix:
			apix = options.apix
		
		lower_bound = -1 * options.tiltrange
		upper_bound = options.tiltrange
	
		nslices = options.nslices
		tiltstep = round(old_div(( float(upper_bound) - float(lower_bound) ), float(nslices - 1)),2)	
	
		#extraslices = 0
		#if options.fillwedge:
		#	extraslices = int( ( 90.0 - float( upper_bound ) ) / tiltstep )
			
			
		#nslices = int(round((upper_bound - lower_bound)/ options.tiltstep))	
		
		#print "\n\nBBBBBBBBBB\n\nThe apix of the simulated ptcl is", apix
		#print "\n\nBBBBBBBBBB\n\n"
	
		px = 0
		if options.applyctf and not options.set2tiltaxis:
			'''
			Random distance 'px' of the particle to the tilt axis at tilt=0, used for CTF simulation.
			Units in microns since it will be used to simulate the effects of having a "defocus gradient"
			with tilt angle, in microns.
			The center of a particle cannot be right at the edge of the tomogram; it has to be
			at least ptcl_size/2 away from it.
			'''
			px = random.uniform(-1* options.gridholesize/2.0 + old_div((apix*image['nx']/2.0),10000), old_div(options.gridholesize,2.0) - old_div((apix*image['nx']/2.0),10000))			
			
		pz=0
		if options.icethickness > old_div(image['nx'],2.0):
			'''
			Beware, --icethickness supplied in microns
			'''
			coordz = random.randint(0 + old_div(image['nx'],2), int(round(old_div(options.icethickness*10000,apix) - old_div(image['nx'],2) )) )
		
			'''
			Calculate the particle's distance 'pz' in microns to the mid section of the tomogram in Z 
			(1/2 of 'nz' or --icethickness). If in-focus occurs at nz/2, then the relative position of 
			a particle to this plane (above or below in the ice) will affect its defocus.
			'''
		
			pz = old_div(coordz*apix,10000) - old_div(options.icethickness,2)	
		
		'''
		Calculate coordx from px since it's shifted by half the gridholesize (the position of the tilt axis)
		and convert to pixels.
		For coordy, generate a random coordinate afresh.
		For coordz, generate a random coordinate within --icethickness range
		'''
		coordx = int( round(  old_div(10000*(px + old_div(options.gridholesize,2)), apix)))
		coordy = random.randint(0 + old_div(image['nx'],2), int(round(old_div(options.gridholesize*10000,apix) - old_div(image['nx'],2) )) )									#random distance in Y of the particle's center from the bottom edge in the XY plane, at tilt=0
		coordz = int( round( old_div(image['nx'],2.0) ) )
		
		if options.set2tiltaxis:
			coordx = 0
			
		sptcoords = tuple([coordx, coordy, coordz])
		
		print("Spt coords are", sptcoords)	

		alt = lower_bound
		raw_projections = []
		ctfed_projections = []
	
		randT = image['sptsim_randT']
	
		print("\n\nWill process particle i", i)
		print("Which was been assigned coordinates")
		print("Nslices are", nslices)
		print("Lower bound is", lower_bound)
	
		finalprjsRAW = finalprjsED = ''

		print("The 3d image is", image)
		print("Its size is",image['nx'],image['ny'],image['nz'])
			
		#tiltangles = []
		#for j in range( nslices ):						#Extra 'noise' slices are 0 if --fillwedge is off. Calculated above if on.
		

		prjindx=0
		for realalt in tiltangles:	
			#realalt = alt + j*tiltstep
		
			#tiltangles.append(realalt)
		
			#print "Real alt is", realalt
			#print "Iterating over slices. Am in tiltstep, slice, alt", tiltstep,j,realalt
		
		
			#Generate the projection orientation for each picture in the tilt series
			
			#t = Transform({'type':'eman','az':90,'alt':realalt,'phi':-90})		#Correct
			t = Transform({'type':'eman','az':-90,'alt':realalt,'phi':90})	#Alternative correct, -90 +90 are interchangeable for phi and az
			#trecon = Transform({'type':'eman','az':0,'alt':realalt,'phi':0}) #Garbage
			if options.tiltaxis == 'x':
				t = Transform({'type':'eman','az':0,'alt':realalt,'phi':0})				
		
			if dimension == 2:
				
				t = Transform({'type':'eman','az':realalt,'alt':0,'phi':0})
				#print "\n\nUSING 2D transform!!!!", t
				#sys.exit()
			
			'''			
			Calculate the defocus shift 'dz' per tilt image, per particle, depending on the 
			particle's position relative to the tilt axis (located at half of the Y dimension
			of an image frame; in this case, 1/2 of --gridholesize), and to the tomograms midsection
			in Z (in this case 1/2 of --icethickness or 'nz').
			There's a dz component coming from px, 'dzx'. 
			There's a dz component coming from pz, 'dzz'.
			
			For particles left of the tilt axis dzx is positive with negative tilt angles (i.e., defocus increases; the particle gets more defocused) 
			For particles right of the tilt axis dzx is negative with negative tilt angles (i.e., defocus decreases; the particle is less defocused) 
			For particles left of the tilt axis dzx is negative with positive tilt angles (i.e., defocus decreases; the particle gets less defocused) 
			For particles right of the tilt axis dzx is positive with positive tilt angles (i.e., defocus increases; the particle gets more defocused) 
			
			Particles above midZ will always be less defocused than if they were exactly at midZ, but this -dzz factor decreases with tilt angle.
				Therefore, for positive pz, dzz is negative (less defocused). 
			Particles above midZ will always be more defocused than if they were exactly at midZ, but this +dzz factor decreases with tilt angle.
				Therefore, for negative pz, dzz is positive (more defocused).
			'''
			dzx = px * numpy.sin( math.radians(realalt) )	
			dzz = -1 * pz * numpy.cos( math.radians(realalt) )			
			
			dz = dzx + dzz
			
			#print "Px is", px
			#print "And angle is", realalt
			#print "Therefore sin(alt) is", numpy.sin(realalt)
			print("And thus dz is", dz)
			defocus = options.defocus + dz
			print("So the final defocus to use is", defocus)
			
			
			
			#prj = image.process("misc.directional_sum",{"axis":"z"})
			print("\nprojecting from",t,realalt)
	
			prj = image.project("standard",t)
			
			print("projection done")
			
			'''
			if options.fillwedge and j > nslices:
	 			
	 			print "\n(e2spt_simulation.py) I'm adding an extra slice with just noise" + str(j-nslices) + '/' + str(extraslices)
				nx = image['nx']
				ny = image['ny']
			
				noise = test_image(1,size=(nx,ny))			
				noise2 = noise.process("filter.lowpass.gauss",{"cutoff_abs":.25})
				noise.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.75})
				
				finalnoise = ( noise*3 + noise2*3 )
			
				prj = finalnoise
			'''	
			
			terror=txerror=tyerror=0
			
			if options.terror and not options.txerror:
				txerror = random.randint( -1*options.terror, options.terror )
			elif options.txerror:
				txerror = random.randint( -1*options.txerror, options.txerror )
				
			if options.terror and not options.tyerror:
				tyerror = random.randint( -1*options.terror, options.terror )
			elif options.tyerror:
				tyerror = random.randint( -1*options.tyerror, options.tyerror )
			
			if txerror or tyerror:
				prj.translate(txerror,tyerror,0)
			
			prj.set_attr('xform.projection',t)
			prj['apix_x']=apix
			prj['apix_y']=apix
			prj['spt_tiltangle']=realalt
			prj['spt_tiltaxis']=options.tiltaxis
			prj['spt_txerror'] = txerror
			prj['spt_tyerror'] = tyerror
			prj['ptcl_source_coord']=sptcoords
			prj['spt_prj_indx'] = prjindx
		
			#print "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\nThe size of the prj is", prj['nx']
			#print "\n"
			
			prj.process_inplace('normalize.edgemean')
		
			if options.saveprjs:
				finalprjsRAW = outname.replace('.hdf', '_ptcl' + str(i).zfill(len(str(nslices))) + '_prjsRAW.hdf')
				#if options.path + '/' in outname:
				#	finalprjsRAW = outname.replace('.hdf', '_ptcl' + str(i).zfill(len(str(nslices))) + '_prjsRAW.hdf')
				
				finalprjsRAW = finalprjsRAW.replace('_preproc','')	
				prj.write_image( finalprjsRAW , prjindx )					#Write projections stack for particle i
				
			raw_projections.append(prj)
			
			if options.invert:
				prj.mult(-1)
			
			prj_r = prj
			
			ctffactor=1.0
			if options.applyctf:
				ctffactor=2.0
			
			if options.snr and options.snr != 0.0 and options.snr != '0.0' and options.snr != '0' and dimension == 3:
				prj_r = noiseit( prj_r, options, nslices, outname, i, ctffactor ) 
				
			prj_fft = prj.do_fft()
			
			
			#print "Sizes of prj and prj_ftt are", prj['nx'],prj['ny'],prj_fft['nx'],prj_fft['ny']
		
											#Reverse the contrast, as in "authentic" cryoEM data		
			if options.applyctf:
				ctf = EMAN2Ctf()
				ctf.from_dict({ 'defocus': defocus, 'bfactor': options.bfactor ,'ampcont': options.ampcont ,'apix':apix, 'voltage':options.voltage, 'cs':options.cs })	
				prj_ctf = prj_fft.copy()	
				ctf.compute_2d_complex(prj_ctf,Ctf.CtfType.CTF_AMP)
				prj_fft.mult(prj_ctf)
		
				prj_r = prj_fft.do_ift()							#Go back to real space
				prj_r['ctf'] = ctf
		
				if options.snr and options.snr != 0.0 and options.snr != '0.0' and options.snr != '0' and dimension == 3:
					prj_r = noiseit( prj_r, options, nslices, outname, i, ctffactor )

			ctfed_projections.append(prj_r)
			#print "Appended ctfed prj in slice j", j
		
			print("options.applyctf", options.applyctf)
			print("should save edited prjs...")
			if options.saveprjs and (options.applyctf or options.snr):
				finalprjsED = outname.replace('.hdf', '_ptcl' + str(i).zfill(len(str(nslices))) + '_prjsEDITED.hdf')
				#if options.path + '/' in outname:
				#	finalprjsED = outname.replace('.hdf', '_ptcl' + str(i).zfill(len(str(nslices))) + '_prjsEDITED.hdf') 
				
				finalprjsED = finalprjsED.replace('_preproc','')
				prj_r.write_image( finalprjsED , prjindx)	
				print("wrote edited prj to %s, indx %d" %( finalprjsED, prjindx ))
	
			prjindx += 1

		#print "The box for IMAGE is with image.get_xsize", image.get_xsize()
		#print "Whereas with image['nx']", image['nx']

		box = max(int(image['nx']),int(image['ny']),int(image['nz']))
		
		#print "!!!!!!!!!!!!!!!Therefore, box is", box
		
		mode='gauss_2'
		if options.reconstructor:
			if len(options.reconstructor) > 1:
				if 'mode' in options.reconstructor[-1]:
					mode=options.reconstructor[-1]['mode']
					
					print("\nThe reconstructor mode has been changed from default to", mode)
					#sys.exit()
					
		
		r = Reconstructors.get(options.reconstructor[0],{'size':(box,box,box),'sym':'c1','verbose':True,'mode':mode})
		
		if dimension == 2:
			print("Boxsize to set up 2D reconstructor is", box,box)
			r = Reconstructors.get(options.reconstructor[0],{'size':(box,box,1),'sym':'c1','verbose':True,'mode':mode})

		#
		#r = Reconstructors.get(options.reconstructor[0],options.reconstructor[1])
		r.setup()
	
		#print "There are these many projections to add to backprojection after all processing", len(ctfed_projections)
	
		k=0
		for p in ctfed_projections:
			#print "Adding projection k", k
			#print "Whose min and max are", p['minimum'], p['maximum']
			#print "The size of the prj to insert is", p['nx']
			pm = r.preprocess_slice(p,p['xform.projection'])
			r.insert_slice(pm,pm['xform.projection'],1.0)
			k+=1
		
		#print "\n\n!!!!!!Will reconstruct the volume for particle i",i
	
		rec = r.finish(True)
		
		print("\n(e2spt_simulation) I have finished simulating particle number", i)
		print("\n")
		#print "The mean of the reconstructed particle is", rec['mean']
		#mname = parameters['model'].split('/')[-1].split('.')[0]
		#name = 'rec_' + mname + '#' + str(i).zfill(len(str(len(particles)))) + '.hdf'
	
		rec['apix_x']=apix
		rec['apix_y']=apix
		rec['apix_z']=apix
		rec['sptsim_randT'] = randT
		rec['origin_x']=0
		rec['origin_y']=0
		rec['origin_z']=0
	
		print("sptcoords for header are", sptcoords)
		rec['ptcl_source_coord']=sptcoords
		
		rec['spt_tiltangles'] = tiltangles
		rec['spt_tiltaxis'] = options.tiltaxis
		
		#print "The apix of rec is", rec['apix_x']
		
		print("\nThe outname to write the particle i", i) 
		print("is", outname)
		print("\n\n")
		
		#finaloutname = options.path + '/' + outname
		#finaloutname.replace('_preproc','')
		#print "is, finaloutname", finaloutname

		print("rec to return is", rec)
			
		#rec.write_image(finaloutname,i)

		#if options.tomogramoutput:
		#	py = random.uniform(0 + image['nx']/2, options.gridholesize - image['nx']/2)									#random distance in Y of the particle's center from the bottom edge in the XY plane, at tilt=0
		#	pz = random.uniform(0 + image['nx']/2, options.icethickness - image['nx']/2) 		
		#	return {'ptcl':rec,'px':px,'py':py,'pz':pz}
		#	
		#else:
		
		#if options.clip:
		#	box = int(options.clip)
		#	print "\nActually, because of clip$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n options.clip is", options.clip
		#	print "Therefore, box is", box
			
		#print "The boxsize to use is", options.clip
			
		#if finalprjsRAW and options.clip and int(options.clip) > int(box):
		#	print "\n\n\n\n\n\n\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@prjs raw need to be clipped! to a size of", options.clip
			
		#	finalprjsRAWclipped = finalprjsRAW.replace('.hdf','_clipped.hdf')  
		#	clipcmd1 = 'e2proc2d.py ' + finalprjsRAW + ' ' + finalprjsRAWclipped + ' --clip=' + str(options.clip) + ',' + str(options.clip) + ' && mv ' + finalprjsRAWclipped + ' ' + finalprjsRAW
		#	print "with cmmand", clipcmd1
			
		#	
		#	p=subprocess.Popen( clipcmd1, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		#	text=p.communicate()	
		#	p.stdout.close()
			
		#if finalprjsED and options.clip and int(options.clip) > int(box):
		#	print "\n\n\n\n\n\n\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@prjs ED need to be clipped! to a size of", options.clip
			
		#	#finalprjsEDclipped = finalprjsED.replace('.hdf','_clipped.hdf')  
			
		#	clipcmd2 = 'e2proc2d.py ' + finalprjsED + ' ' + finalprjsEDclipped + ' --clip=' + str(options.clip) + ',' + str(options.clip) + ' && mv ' + finalprjsEDclipped + ' ' + finalprjsED 
		#	print "with cmmand", clipcmd2
			
			
			
		#	p=subprocess.Popen( clipcmd2, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		#	text=p.communicate()	
		#	p.stdout.close()
		
		#box = image.get_xsize()
		
		print("\n\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\nThe final boxsize of rec is", rec['nx'], rec['ny'], rec['nz']) 
		print("and box is", box)
		return { classoptions['ptclnum']:rec }


jsonclasses["SubtomoSimTask"]=SubtomoSimTask.from_jsondict


def noiseit( prj_r, options, nslices, outname, i, ctffactor ):
	
	noise = EMData()
	nx=prj_r['nx']
	ny=prj_r['ny']

	#print "I will make noise"
	#noise = 'noise string'
	#print "Noise is", noise
	
	noise = test_image(1,size=(nx,ny))
	#print "noise now is img", noise

	#induce spatial correlations
	#noise2 = noise.process("filter.lowpass.gauss",{"cutoff_abs":.25})
	#noise.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.75})
	
	#noise = ( noise*3 + noise2*3 )
	
	if options.savenoise:
		noiseStackName = outname.replace('.hdf', '_ptcl' + str(i).zfill(len(str(nslices))) + '_NOISE.hdf')
		noise.write_image( noiseStackName, -1 )
	
	noise *= old_div(float( options.snr ), ctffactor)

	#if noise:
	#print "I will add noise"
	#noise.process_inplace('normalize')
	
	#prj_r.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})
	#prj_r.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.75})
	
	fractionationfactor = old_div(61.0,nslices)		#At snr = 10, simulated subtomograms look like empirical ones for +-60 deg data collection range
											#using 2 deg tilt step. If 61 slices go into each subtomo, then the fractionation factor
											#Will be 1. Otherwise, if nslices is > 61 the signal in each slice will be diluted.
											#If nslices < 1, the signal in each slice will be enhanced. In the end, regardless of the nslices value, 
											#subtomograms will always have the same amount of signal instead of signal depending on number of images.
	prj_r.mult( fractionationfactor )
	prj_r.add(noise)
	
	
	
	#elif options.snr:
	#	print "WARNING: You specified snr but there's no noise to add, apparently!"
	return prj_r


def get_results(etc,tids,options):
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
				r=etc.get_results(tidsleft[i])				#Results for a completed task
				
				if r:
					#print "r is", r
					ptcl=r[0].classoptions["ptclnum"]		#Get the particle number from the task rather than trying to work back to it
					results[ptcl] = r[1]					
				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if options.verbose:
			print("  %d tasks, %d complete, %d waiting to start        \r"%(len(tids),ncomplete,nwait))
			sys.stdout.flush()
	
		if len(tidsleft)==0: break
		
	return results


if __name__ == '__main__':
	main()


'''
NOTES
			#prj_r.process_inplace('math.addnoise',{'noise':options.snr})
				
				#prj_n = EMData(nx,ny)
				#for i in range(options.snr):				
				#	prj_n.process_inplace(options.noiseproc[0],options.noiseproc[1])
				#	prj_r = prj_r + prj_n
				#prj_n.write_image('NOISE_ptcl' + str(i).zfill(len(str(nslices))) + '.hdf',j)

'''
