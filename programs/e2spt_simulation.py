#!/usr/bin/env python

'''
====================
Author: Jesus Galaz - 2011, Last update: July/08/2013
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

from optparse import OptionParser
from EMAN2 import *
from sys import argv
import EMAN2
import heapq
import operator
import random
import numpy
import math
from EMAN2jsondb import JSTask,jsonclasses

def main():

	usage = """e2spt_simulation.py <options>.
	 This program DEPENDS on e2spt_classaverage.py from which it imports the sptmakepath function.
	 The options should be supplied in "--option=value", replacing "option" for a valid option name, and "value" for an acceptable value for that option. This program produces simulated sub volumes in random orientations from a given PDB or EM file. The output is ALWAYS in HDF format, since it's the only format supported by E2SPT programs.
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	
	#parser.add_argument("--interpolator",default='',help="""What interpolation scheme 
	#	to use for reconstruction. Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 
	#	'gauss_5_slow', 'gypergeom_5', experimental""")
	
	parser.add_argument("--path",type=str,default=None,help="Directory to store results in. The default is a numbered series of directories containing the prefix 'sptsim'; for example, sptsim_02 will be the directory by default if 'sptsim_01' already exists.")
	#parser.add_argument("--output",type=str,default=None,help="Name of the output stack for the simulated subtomograms.")
	parser.add_argument("--randstack",type=str,default=None,help="If you already have a stack of particles (presumably in random orientations) you can supply it here.")
	parser.add_argument("--tomogramoutput",type=str,default=None,help="This will generate a simulated tilt series and tomogram containing the entire set of subtomograms.")

	parser.add_argument("--input", type=str, help="""The name of the input volume from which simulated subtomograms will be generated. 
							The output will be in HDF format, since volume stack support is required. The input CAN be PDB, MRC or and HDF stack. 
							If the input file is PDB or MRC, a version of the supplied model will be written out in HDF format.
							If the input file is a stack, simulatd subvolumes will be generated from each model in the stack and written to different output stacks.
							For example, if the input file contains models A and B, two output stacks with simulated subvolumes will be generated.""", default=None)
				
	parser.add_argument("--filter",type=str,help="""A filter (as in a processor from e2proc3d.py) apply to the model before generating simulated particles from it.
							Type 'e2help.py processors' at the command line and find the options availbale from the processors list)""",default=None)
	
	parser.add_argument("--shrink", type=int,default=0,help="Optionally shrink the input volume before the simulation if you want binned/down-sampled subtomograms.")
	parser.add_argument("--verbose", "-v", type=int, dest="verbose", action="store", metavar="n", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_argument("--nptcls", type=int,default=10,help="Number of simulated subtomograms tu generate per referece.")
	parser.add_argument("--txrange", type=int,default=None,help="""Maximum number of pixels to randomly translate each subtomogram in X. The random translation will be picked between -txrange and +txrange. 
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--tyrange", type=int,default=None,help="""Maximum number of pixels to randomly translate each subtomogram in Y. The random translation will be picked between -tyrange and +tyrange.
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--tzrange", type=int,default=None,help="""Maximum number of pixels to randomly translate each subtomogram in Z. The random translation will be picked between -tzrange and +tzrange.
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--transrange", type=int,default=0,help="""Maximum number of pixels to randomly translate each subtomogram in all X, Y and Z. 
									The random translation will be picked between -transrage and +transrange; --txrange, --tyrange and --tzrange overwrite --transrange for each specified direction.""")
	
	#parser.add_argument("--tiltstep", type=float,default=5.0,help="Degrees between each image in the simulated tilt series for each subtomogram.")
	parser.add_argument("--tiltrange", type=float,default=60,help="""Maximum angular value at which the highest tilt picture will be simulated. Projections will be simulated from -tiltrange to +titlrange. 
									For example, if simulating a tilt series collected from -60 to 60 degrees, enter a --tiltrange value of 60. 
									Note that this parameter will determine the size of the missing wedge.""")
	parser.add_argument("--applyctf", action="store_true",default=False,help="If on, it applies ctf to the projections in the simulated tilt series based on defocus, cs, and voltage parameters.")
	
	parser.add_argument("--saveorthostack", action="store_true",default=False,help="If on, --nptcls is ignored and you get 3 subtomograms (simulated from the model supplied) which are orthogonal to each other.")

	parser.add_argument("--defocus", type=float,default=3.0,help="Intended defocus at the tilt axis (in microns) for the simulated tilt series.")
	parser.add_argument("--voltage", type=int,default=200,help="Voltage of the microscope, used to simulate the ctf added to the subtomograms.")
	parser.add_argument("--cs", type=float,default=2.1,help="Cs of the microscope, used to simulate the ctf added to the subtomograms.")

	parser.add_argument("--gridholesize", type=float,default=0.5,help="""Size of the carbon hole in micrometers for the simulated grid (this will determine the shifts in defocus for each particle at 
									each tilt step, depending on the position of the particle respect to the tilt axis; the tilt axis by convention goes parallel to Y through the middle of the tomogram.
									of the subvolumes within the tomogram are assigned randomly.""")
	parser.add_argument("--icethickness", type=float,default=100,help="If --tomogramoutput is supplied, this will define the size of Z in PIXELS for the simulated tomogram.")

	parser.add_argument("--saverandstack", action="store_true",default=False,help="""Save 
		the stack of randomly oriented particles, before subtomogram simulation (before the 
		missing wedge and noise are added).""")

	parser.add_argument("--nosim", action="store_true",default=False,help="""If on
		the program will generate stacks of "perfect particles" in different random 
		orientations, but with no missing wedge, no noise, no ctf parameters, etc.
		The output randstack.hdf will be identical to simptcls.hdf""") 
	
	parser.add_argument("--saveprjs", action="store_true",default=False,help="Save the projections (the 'tilt series') for each simulated subtomogram.")

	parser.add_argument("--reconstructor", type=str,default="fourier",help="""The reconstructor 
		to use to reconstruct the tilt series into a tomogram. Type 'e2help.py reconstructors' 
		at the command line to see all options and parameters available.
		To specify the interpolation scheme for the fourier reconstruction, specify 'mode'.
		Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 
		'gauss_5_slow', 'gypergeom_5', 'experimental'.
		For example --reconstructor=fourier:mode=gauss_5 """)									
											
	parser.add_argument("--pad", type=float,default=0.0,help="""If on, it will increase the 
		box size of the model BEFORE generating projections and doing 3D reconstruction of 
		simulated sutomograms. Make sure to supply --finalboxsize to clip the simulated 
		subtomograms back to a reasonable box size.""")								
	#parser.add_argument("--noiseproc",type=str,help="A noise processor to be applied to the individual projections of the simulated tilt series",default=None)
	
	parser.add_argument("--finalboxsize", type=str,default='',help="""The final box size 
		to clip the subtomograms to. Enter an integer, or a positive number followed by 'x' to
		be used as an 'expansion factor'. For example, if the model in --input has an original
		boxsize of 64 and you enter --finalboxsize=1.5x, then the finalboxsize of the
		simulated subtomograms will be 1.5*64, that is, 96.
		If you provide --shrink, for example, --shrink=2, then the finalboxsize of the 
		subtomograms will be 64/2 * 1.5 = 48.""")								

	parser.add_argument("--snr",type=float,help="Weighing noise factor for noise added to the image.",default=0)
	#parser.add_argument("--addnoise",action="store_true",default=False,help="If on, it adds random noise to the particles")
	
	parser.add_argument("--sym",type=str,default='c1',help="If your particle is symmetrical, you should randomize it's orientation withing the asymmetric unit only. Thus, provide the symmetry.")

	parser.add_argument("--notrandomize",action="store_true",default=False,help="This will prevent the simulated particles from being rotated and translated into random orientations.")
	parser.add_argument("--simref",action="store_true",default=False,help="This will make a simulated particle in the same orientation as the original input (or reference).")
	parser.add_argument("--negativecontrast",action="store_true",default=False,help="This will make the simulated particles be like real EM data before contrast reversal. Otherwise, 'white protein' (positive density values) will be used.")

	parser.add_argument("--nslices", type=int,default=61,help="""This will determine the tilt step between slices, depending on tiltrange.
									For example, to simulate a 2 deg tilt step supply --nslices=61 --tiltrange=60.
									Recall that --tiltrange goes from - to + the supplied value, and that there is a central slice or projection at 0 deg,
									for symmetrical tilt series.""")	

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--parallel",  help="Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel", default='thread:1')
	
	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)

	print "e2spt_simulation received these many slices", options.nslices

	'''
	Parse the options
	'''
	if options.filter:
		options.filter = parsemodopt(options.filter)

	if options.reconstructor:
		options.reconstructor = parsemodopt(options.reconstructor)
	

	if options.input and options.randstack:
		print """\n(e2spt_simulation)WARNING: No point in supplying --input and --randstack simultaneously.
		They are mutually exclusive. If --randstack is provided, --input is ignored,
		because --randstack becomes --input.\n
		--input was""",options.input
		print '--randstack was', options.randstack
	
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'sptsim')
	
	rootpath = os.getcwd()
	
	#if rootpath not in options.path:
	#	options.path = rootpath + '/' + options.path
		
	randptcls = {}
	
	if options.randstack:
		print "\n\nI will not generate a randstack but will read it from", options.randstack

		randstackbase = os.path.basename( options.randstack )
		randstackcopy = options.path + '/' + randstackbase

		copycmd = 'e2proc3d.py ' + options.randstack + ' ' + randstackcopy
		
		p=subprocess.Popen( copycmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
		
		
		nr=EMUtil.get_image_count(randstackcopy)
		print "There are these many particles in the randstack", nr							#ATTENTION: Randstack might still need box size changes and padding...
		
		options.input = randstackcopy
		
		if options.filter or float(options.pad) > 1.0 or int(options.shrink) > 1:
			ret = preprocess(options,randstackcopy) 
			randstackcopy = ret[-1]

		print "Randstackcopy to read particles from is", randstackcopy
		
		#stackname = rootpath + '/' + os.path.basename(options.randstack)
		
		simptclsname = options.path + '/simptcls.hdf'
		if options.nosim:
			
			clip2cmd = 'e2proc3d.py ' + randstackcopy + ' ' + randstackcopy + ' --clip=' + str(options.finalboxsize)
			p=subprocess.Popen( clip2cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()
	
			os.system('cp ' + options.input + ' ' + simptclsname)
		
			simrefname = options.path + '/simmodel.hdf'
			
			simrefcmd = 'e2proc3d.py ' + randstackcopy + ' ' + simrefname + ' --first=0 --last=0'
			p=subprocess.Popen( simrefcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()
			
			
		else:
			for np in range(nr):
				a=EMData(randstackcopy,np)
				#>>randptcls.append(a)
				randptcls.update({np:a})
			
			subtomosim(options, randptcls, simptclsname )
		
		
		#if options.simref and options.input and not options.nosim:
		#	model = EMData(options.input,0)
		#	
		#	name = options.path + '/' + options.input.replace('.hdf','_SIM.hdf').split('/')[-1]
		#	model['sptsim_randT'] = Transform()
		#	model['xform.align3d'] = Transform()
		#	ret = subtomosim(options,[model],name)
		#
		#	if ret == 1:
		#		#os.system('e2proc3d.py ' + name + ' ' + name + ' --clip=' + str(options.finalboxsize) + ' --first=' + str(i) + ' --last=' + str(i))
		#		
		#		clip3cmd = 'e2proc3d.py ' + name + ' ' + name + ' --clip=' + str(options.finalboxsize) + ' --first=' + str(0) + ' --last=' + str(0)
		#		p=subprocess.Popen( clip3cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		#		text=p.communicate()	
		#		p.stdout.close()	
	
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
				print "These many particles will be simulated, for each of the supplied references/models", options.nptcls
				print "There are these many references/models", nrefs
		
			originalpath = options.path
			kkk=0
	
			originalinput = options.input
		
			for i in range(nrefs):
				if options.verbose:
					print "\n\nGenerating simulated subtomograms for reference number", kkk
					
				if nrefs>1:
					modelfilename = originalinput.split('/')[-1].replace('.hdf','_model' + str(i).zfill(2) + '.hdf')

					options.path = originalpath + '/model' + str(i).zfill(2)

					os.system('mkdir ' + options.path)
					#cmd = 'e2proc3d.py '  + options.input + ' ' + options.path + '/' + modelfilename + ' --first=' + str(i) + ' --last=' + str(i) + ' --append'
	
					os.system('e2proc3d.py '  + originalinput + ' ' + options.path + '/' + modelfilename + ' --first=' + str(i) + ' --last=' + str(i) + ' --append')
					print "This is the command to create the model"
					print 'e2proc3d.py '  + originalinput + ' ' + options.path + '/' + modelfilename + ' --first=' + str(i) + ' --last=' + str(i) + ' --append'
	
					options.input = options.path + '/' + modelfilename
					tag = str(i).zfill(len(str(nrefs)))
				
				
				modelhdr = EMData(options.input,0,True)
	
				'''
				Make model's box cubical if it isn't
				'''
			
				if modelhdr['nx'] != modelhdr['ny'] or modelhdr['nx'] != modelhdr['nz'] or modelhdr['ny'] != modelhdr['nz']:	
					sptmakecube( options )
					modelhdr = EMData(options.input,0,True)	
			
				'''
				Preprocess model if necessary
				'''
				if options.filter or options.pad or options.shrink:
					ret = preprocess(options,options.input)
					options = ret[0]
					
				newsize = max( int(modelhdr['nx']), int(modelhdr['ny']), int(modelhdr['nz']) )
				if options.finalboxsize:
					if 'x' not in options.finalboxsize:
						print "\nThere is no x factor in finalboxsize"
						finalboxsize = int( float(options.finalboxsize) )
					else:
						options.finalboxsize = int( round( float( options.finalboxsize.replace('x','') ) * newsize ) )	
						
						print "Therefore the expanded boxisze is", options.finalboxsize
						print "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n\n\n\n\n\n"
	
					#os.system('e2proc3d.py ' + options.input + ' ' + options.input + ' --clip=' + str(options.finalboxsize) + ' --first=' + str(i) + ' --last=' + str(i))
					
					clip2cmd = 'e2proc3d.py ' + options.input + ' ' + options.input + ' --clip=' + str(options.finalboxsize) + ' --first=' + str(i) + ' --last=' + str(i)
					p=subprocess.Popen( clip2cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					text=p.communicate()	
					p.stdout.close()
					
				else:
					options.finalboxsize = max( int(modelhdr['nx']), int(modelhdr['ny']), int(modelhdr['nz']) )

				model = EMData(options.input,0)
							
				retrand = randomizer(options, model,tag)
				randptcls = retrand[0]
				randstackname = retrand[-1]
				
				simptclsname = options.path + '/simptcls.hdf'
				
				print "\n\n\n\n\n\(e2spt_simulation) before subtomosim, simptclsname is", simptclsname
				print "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS\n\n\n\n\n\n\n\n"
				if not options.nosim:
					subtomosim(options,randptcls,simptclsname)		
				else:
					os.system('cp ' + randstackname + ' ' + simptclsname)
					
					simrefname = options.path + '/simmodel.hdf'
			
					simrefcmd = 'e2proc3d.py ' + options.input + ' ' + simrefname + ' --first=0 --last=0'
					p=subprocess.Popen( simrefcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					text=p.communicate()	
					p.stdout.close()
					
				if options.saveorthostack:
					orthoptcls = orthostack( options, model )
					
					if not options.nosim:
						simorthoptclsname = options.path + '/orthoptcls.hdf'
						subtomosim(options,orthoptcls,simorthoptclsname)		

				#if ret == 1:
				#	os.system('e2proc3d.py ' + options.input + ' ' + options.input + ' --clip=' + str(options.finalboxsize) + ' --first=' + str(i) + ' --last=' + str(i))	

				#if options.simref and not options.nosim:
				#	#name = options.input.replace('.hdf','_SIM.hdf').split('/')[-1]
				#	
				#	simmodelname = options.path + '/simmodel.hdf'
				#	model['sptsim_randT'] = Transform()
				#	model['xform.align3d'] = Transform()
				#	ret = subtomosim(options,[model],simmodelname)
				#
				#	if ret == 1:
				#		#os.system('e2proc3d.py ' + name + ' ' + name + ' --clip=' + str(options.finalboxsize) + ' --first=' + str(i) + ' --last=' + str(i))
				#		
				#		clip3cmd = 'e2proc3d.py ' + simmodelname + ' ' + simmodelname + ' --clip=' + str(options.finalboxsize) + ' --first=' + str(i) + ' --last=' + str(i)
				#		p=subprocess.Popen( clip3cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				#		text=p.communicate()	
				#		p.stdout.close()
				kkk+=1
	
	E2end(logger)
					
	return()




def sptfixformat( options ):
	check=0	
	
	if '.pdb' in options.input:
		pdbmodel = options.input
		os.system('cp ' + pdbmodel + ' ' + options.path)
		pdbmodel = options.path + '/' + pdbmodel.split('/')[-1]
		mrcmodel = pdbmodel.replace('.pdb','.mrc')
		os.system('e2pdb2mrc.py ' + pdbmodel + ' ' + mrcmodel + ' && rm ' + pdbmodel)
		options.input = mrcmodel
		check=1

	if '.mrc' in options.input:
		mrcmodel = options.input
		if check==0:
			os.system('cp ' + mrcmodel + ' ' + options.path)
			mrcmodel = options.path + '/' + mrcmodel.split('/')[-1]
		hdfmodel = mrcmodel.replace('.mrc','.hdf')
		os.system('e2proc3d.py ' + options.input + ' ' + hdfmodel + ' && rm ' + mrcmodel)
		options.input = hdfmodel
		check=1
			
	if '.hdf' in options.input:
		hdfmodel = options.input
		if check == 0:
			os.system('cp ' + hdfmodel + ' ' + options.path)
			hdfmodel = options.path + '/' + hdfmodel.split('/')[-1]
			options.input = hdfmodel
		workname = hdfmodel.replace('.hdf','_sptsimMODEL.hdf')
		#if nrefs > 1:
		#	workname = hdfmodel.replace('.hdf','_sptsimMODELS.hdf')

		os.system('cp ' + hdfmodel + ' ' + workname + ' && rm ' + hdfmodel)
		options.input = workname

	return options


def sptmakecube( options ):
	
	newsize = max( int(modelhdr['nx']), int(modelhdr['ny']), int(modelhdr['nz']) )
	print "\n\n\n\nNEWSIZE will be", newsize

	fixboxcmd='e2proc3d.py ' + options.input + ' ' + options.input + ' --clip=' + str(newsize)
	p=subprocess.Popen( fixboxcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
			
	return


def preprocess(options,stack):
		
	preprocessed = stack.replace('.hdf','_preproc.hdf')
	
	os.system('e2proc3d.py	' + stack + ' ' + preprocessed)
	
	hdr = EMData(stack,0,True)
	newsize=hdr['nx']
	nf = EMUtil.get_image_count(stack)
	
	if options.shrink and int(options.shrink) > 1:
		newsize = newsize/options.shrink	
		if newsize % 2:
			newsize += 1

		#os.system('e2proc3d.py ' + options.input + ' ' + options.input + ' --process=math.meanshrink:n=' + str(options.shrink))
		
		shrinkcmd='e2proc3d.py ' + preprocessed + ' ' + preprocessed + ' --process=math.meanshrink:n=' + str(options.shrink)
		p=subprocess.Popen( shrinkcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()

	if options.pad and float(options.pad) > 1.0:
		newsize *= options.pad
		newsize = int( round( newsize ) )
		
		clipcmd = 'e2proc3d.py ' + preprocessed + ' ' + preprocessed + ' --clip=' + str(newsize)
		p=subprocess.Popen( clipcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
	
		fixcmd = 'e2fixheaderparam.py --input=' + preprocessed + ' --stem=origin --stemval=0.0'
		p=subprocess.Popen( fixcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
		
		options.input = preprocessed
		
	if options.filter:
		pass #Filter not working for now
		
	

	return [options,preprocessed]
	

'''
====================
RANDOMIZER - Takes a file in .hdf format and generates a set of n particles randomly rotated and translated, where n is
the size of the set and is defined by the user
====================
'''
def randomizer(options, model, tag):
	
	#print "I am inside the RANDOMIZER"
	
	if options.verbose:
		print "You have requested to generate %d particles with random orientations and translations" %(options.nptcls)
	
	randptcls = {}
	
	#randstackname = options.path + '/' + options.input.replace('.hdf','_randst' + tag + '_n' + str(options.nptcls).zfill(len(str(options.nptcls))) + '.hdf').split('/')[-1]
	#if options.output:
	#	randstackname = options.path + '/' + options.output.replace('.hdf','_randst' + tag + '_n' + str(options.nptcls).zfill(len(str(options.nptcls))) + '.hdf').split('/')[-1]
	
	randstackname = options.path + '/randstack.hdf'
	
	print "###############\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n#################\nThe stackname inside RANDOMIZER, is", randstackname
	print "--saverandstack is", options.saverandstack
	print "#####################################\n\n\n\n\n\n\n\n\n\n\n\n\n"
	for i in range(options.nptcls):
		if options.verbose:
			print "I will generate particle #", i

		b = model.copy()
		
		random_transform = Transform()	
		
		if not options.notrandomize:
			if i > 0:
				rand_orient = OrientGens.get("rand",{"n":1, "phitoo":1})						#Generate a random orientation (randomizes all 3 euler angles)
				c1_sym = Symmetries.get("c1")
				random_transform = rand_orient.gen_orientations(c1_sym)[0]
			
				randtx = randty = randtz = 0	
				if options.transrange:
					randtx = random.randrange(-1 * options.transrange, options.transrange)			#Generate random translations
					randty = random.randrange(-1 * options.transrange, options.transrange)
					randtz = random.randrange(-1 * options.transrange, options.transrange)

				if options.txrange:
					randtx = random.randrange(-1 * options.txrange, options.txrange)	
				if options.tyrange:
					randty = random.randrange(-1 * options.tyrange, options.tyrange)
				if options.tzrange:
					randtz = random.randrange(-1 * options.tzrange, options.tzrange)
			
				if randtx or randty or randtz:
					random_transform.translate(randtx, randty, randtz)

			b.transform(random_transform)		

		b['sptsim_randT'] = random_transform
		b['xform.align3d'] = Transform()							#This parameter should be set to the identity transform since it can be used later to determine whether
													#alignment programs can "undo" the random rotation in spt_randT accurately or not
		if options.saverandstack:	

			#print "The stackname to use is", stackname
			#randstackname = options.path + '/' + stackname.split('/')[-1]
			#print "\n\n\n\n\nI will save randstack! Using THE PATH in e2spt_simulation and stackname both of which together are", randstackname
			#print "options path received is", options.path
			#print "Whereas stackname is", stackname
			#print "Therefore stackname.split('/')[-1] is", stackname.split('/')[-1]
			
			
			b['origin_x'] = 0									#Make sure the origin is set to zero, to avoid display issues with Chimera
			b['origin_y'] = 0
			b['origin_z'] = 0
			
			b.write_image(randstackname,i)
			print "Actually, particle written to", randstackname

		#>>randptcls.append(b)

			randptcls.update({i:b})
		if options.verbose:
			print "The random transform applied to it was", random_transform
			
		if options.finalboxsize:
			if int(b['nx']) != int(options.finalboxsize) or int(b['ny']) != int(options.finalboxsize) or int(b['nz']) != int(options.finalboxsize):
				
				
				clipcmdf = 'e2proc3d.py ' + randstackname + ' ' + randstackname + ' --clip=' + str(options.finalboxsize)
		
				p=subprocess.Popen( clipcmdf, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				text=p.communicate()	
				p.stdout.close()
	

	return(randptcls,randstackname)

	
def orthostack(options, model):
	orthoptcls={}
	for i in range(4):
		b = model.copy()
		
		if i==0:
			b['sptsim_randT'] = Transform()
			b['xform.align3d'] = Transform()
			orthoptcls.append(b)
		
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
			orthoptcls.append(b)
			
		if i == 2:
			b.rotate(0,90,90)
			t = Transform({'type':'eman','az':0,'alt':90,'phi':90})
			b['sptsim_randT'] = t
			b['xform.align3d'] = Transform()
			orthoptcls.update({i:b})
		
		#orthostackname = options.path + '/orthostack.hdf'
		
	return( orthoptcls )
	
'''
====================
SUBTOMOSIM takes a set of particles in .hdf format and generates a simulated sub-tomogram for each, using user-defined parameters for tomographic simulation.
The function generates projections for each particle with a user-defined tilt step and missing wedge size (or data collection range),
adds noise and ctf to each projection, randomizes the position of each particle "in the ice" and changes the defocus for each picture in the tilt series accordingly,
and recounstructs a new 3D volume from the simulated tilt series.
====================
'''	
def subtomosim(options,ptcls,outname):
	#print "INSIDE SUBTOMOSIM"
	
	print "\n\n\n\n\n(e2spt_simulation) Outname received in subtomosim", outname
	print "OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO\n\n\n\n\n\n"	
	'''
	Initialize parallelism if being used
	'''
	if options.parallel :
		print "\n\n(e2spt_simulation.py) INITIALIZING PARALLELISM, for this outname (stack, or reference)", outname
		print "\n\n"
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel)
	
	if options.verbose:
		
		print "(e2spt_simulation) There are these many slices to produce to simulate each subtomogram", options.nslices
	
	#outname = stackname.replace('.hdf','_ptcls.hdf')
		
	#
	#if options.output:
	#	outname = options.path + '/' + options.output.replace('.hdf','_ptcls.hdf')
	
		
	#if len(ptcls) == 1 and '_SIM.hdf' in stackname:
	#	outname = stackname.split('/')[-1]
	
	#tomogramdata=[]
	
	#print "\n\n\n%%%%%%%%%%%%%%%%The number of particles are", len(ptcls)
	#print "\n\n\n"
	
	tasks=[]
	#>>for i in range(len(ptcls)):
	for i in ptcls:	
		#if options.parallel:			
		task=SubtomoSimTask(ptcls[i],i,options,outname)
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
		
		key = result.keys()[0]
		
		#if options.path not in outname:
		#	finaloutname = options.path + '/' + outname
		
		#print "\n\n\n\n\n\nn\\n\n\nn\\nn\n\n\n\n\n THe final outname rigth before writing is", finaloutname
		#print "because options.path is", options.path
		#print "and outname is", outname
		#print "And the particle is", results[pn]
		#print "And its type is", type( results[pn])
		#print "And its index is", pn
		
		if i==0:
			#>>print "\n\n(subtomosim) The size of the final particle is",results[pn]['nx'],results[pn]['ny'],results[pn]['nz']
			print "\n\n(subtomosim) The size of the final particle is",result[key]['nx'],result[key]['ny'],result[key]['nz']

		result[key]['origin_x'] = 0									#Make sure the origin is set to zero, to avoid display issues with Chimera
		result[key]['origin_y'] = 0
		result[key]['origin_z'] = 0
		
		#finaloutname = finaloutname.replace('_preproc','')
		result[key].write_image(finaloutname,key)
		#pn+=1
		ii+=1
	
	if options.finalboxsize:
		box = int(options.finalboxsize)
		print "\nActually, because of finalboxsize$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n options.finalboxsize is", options.finalboxsize
		print "Therefore, box of final simulated stack is", box
		
		clipcmdf = 'e2proc3d.py ' + finaloutname + ' ' + finaloutname + ' --clip=' + str(options.finalboxsize)
		print "with cmmand", clipcmdf
		
		p=subprocess.Popen( clipcmdf, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
	
	
	if options.simref and not options.saveorthostack:
		simrefname = options.path + '/simmodel.hdf'
		
		simrefcmd = 'e2proc3d.py ' + finaloutname + ' ' + simrefname + ' --first=0 --last=0'
		p=subprocess.Popen( simrefcmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text=p.communicate()	
		p.stdout.close()
	
	if options.tomogramoutput:
		tomogramsim(options,results)

	
	return
	
	
	
class SubtomoSimTask(JSTask):
	"""This is a task object for the parallelism system. It is responsible for generating a simulated subtomogram."""	
	
	def __init__(self,image,ptclnum,options,outname):
	
		"""fixedimage and image may be actual EMData objects, or ["cache",path,number]
		label is a descriptive string, not actually used in processing
		ptcl is not used in executing the task, but is for reference
		other parameters match command-line options from e2spt_classaverage.py
		Rather than being a string specifying an aligner, 'align' may be passed in as a Transform object, representing a starting orientation for refinement"""
		#data={}
		#data={"fixedimage":fixedimage,"image":image}
		data={"image":image}
		
		JSTask.__init__(self,"SubTomoSim",data,{},"")

		self.classoptions={"options":options,"ptclnum":ptclnum,"outname":outname}
	
	def execute(self,callback=None):
		"""This simulates a subtomogram and saves projections before and after adding noise if instructed to."""
		classoptions=self.classoptions
		options=self.classoptions['options']
		outname = self.classoptions['outname']
		
		i = self.classoptions['ptclnum']
		
		image = self.data['image']
		
		print "\n\n(SubtomoSimTask) Size of the particle for simulation is", image['nx'],image['ny'],image['nz']
		
		outname = self.classoptions['outname']
		
		if options.verbose:
			print "Generating projections for particle #", i

		#apix = ptcls[i]['apix_x']
		apix = image['apix_x']
		
		lower_bound = -1 * options.tiltrange
		upper_bound = options.tiltrange
	
		nslices = options.nslices
		tiltstep = round(( float(upper_bound) - float(lower_bound) )/ float(nslices - 1),2)	
	
		#nslices = int(round((upper_bound - lower_bound)/ options.tiltstep))	
		
		#print "\n\nBBBBBBBBBB\n\nThe apix of the simulated ptcl is", apix
		#print "\n\nBBBBBBBBBB\n\n"
	
		px = random.uniform(-1* options.gridholesize/2 + image['nx']/2, options.gridholesize/2 - image['nx']/2)			#random distance in X of the particle's center from the tilt axis, at tilt=0
																																#The center of a particle cannot be right at the edge of the tomogram; it has to be
																																#at least ptcl_size/2 away from it
		alt = lower_bound
		raw_projections = []
		ctfed_projections = []
	
		randT = image['sptsim_randT']
	
		print "\n\nWill process particle i", i
		print "Nslices are", nslices
		print "Lower bound is", lower_bound
	
		finalprjsRAW = finalprjsED = ''

		print "The 3d image is", image
		print "Its size is",image['nx'],image['ny'],image['nz']
			
		for j in range(nslices):
			realalt = alt + j*tiltstep
		
			#print "Real alt is", realalt
			#print "Iterating over slices. Am in tiltstep, slice, alt", tiltstep,j,realalt
		
			t = Transform({'type':'eman','az':90,'alt':realalt,'phi':0})				#Generate the projection orientation for each picture in the tilt series
		
			dz = -1 * px * numpy.sin(realalt)							#Calculate the defocus shift per picture, per particle, depending on the 
																	#particle's position relative to the tilt axis. For particles left of the tilt axis,
																	#px is negative. With negative alt [left end of the ice down, right end up], 
																	#dz should be negative.
			defocus = options.defocus + dz
			
			
			
			#prj = image.process("misc.directional_sum",{"axis":"z"})
				
			prj = image.project("standard",t)
			prj.set_attr('xform.projection',t)
			prj['apix_x']=apix
			prj['apix_y']=apix
			prj['sptsim_tiltangle']=realalt
		
			#print "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\nThe size of the prj is", prj['nx']
			#print "\n"
			
			prj.process_inplace('normalize')
		
			if options.saveprjs:
				finalprjsRAW = outname.replace('.hdf', '_ptcl' + str(i).zfill(len(str(nslices))) + '_prjsRAW.hdf')
				#if options.path + '/' in outname:
				#	finalprjsRAW = outname.replace('.hdf', '_ptcl' + str(i).zfill(len(str(nslices))) + '_prjsRAW.hdf')
				
				finalprjsRAW = finalprjsRAW.replace('_preproc','')	
				prj.write_image( finalprjsRAW , -1)					#Write projections stack for particle i
				
			
			raw_projections.append(prj)
				
			prj_fft = prj.do_fft()
			
			#print "Sizes of prj and prj_ftt are", prj['nx'],prj['ny'],prj_fft['nx'],prj_fft['ny']
		
			if options.negativecontrast:
				prj_fft.mult(-1)								#Reverse the contrast, as in "authentic" cryoEM data		
		
			if options.applyctf:
				ctf = EMAN2Ctf()
				ctf.from_dict({'defocus':options.defocus,'bfactor':100,'ampcont':0.05,'apix':apix,'voltage':options.voltage,'cs':options.cs})	
				prj_ctf = prj_fft.copy()	
				ctf.compute_2d_complex(prj_ctf,Ctf.CtfType.CTF_AMP)
				prj_fft.mult(prj_ctf)
		
			prj_r = prj_fft.do_ift()							#Go back to real space
			noise = ''
		
			if options.snr and options.snr != 0.0 and options.snr != '0.0' and options.snr != '0':
				nx=prj_r['nx']
				ny=prj_r['ny']
			
				#print "I will make noise"
			
				#noise = 'noise string'
				#print "Noise is", noise
				noise = test_image(1,size=(nx,ny))
				#print "noise now is img", noise
			
				noise2 = noise.process("filter.lowpass.gauss",{"cutoff_abs":.25})
				noise.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.75})
				noise = ( noise*3 + noise2*3 ) * int(options.snr)
			
				if noise:
					#print "I will add noise"
					#noise.process_inplace('normalize')
					
					prj_r.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})
					prj_r.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.75})
					
					fractionationfactor = 61.0/nslices		#At snr = 10, simulated subtomograms look like empirical ones for +-60 deg data collection range
															#using 2 deg tilt step. If 61 slices go into each subtomo, then fractionation factor
															#Will be 1. If nslices is > 61 the signal in each slice will be diluted.
															#If nslices < 1, the signal in each slice will be enhanced. In the end, regardless of the nslices value, 
															#subtomograms will always have the same amount of signal.
					prj_r.mult( fractionationfactor )
					prj_r.add(noise)
			
				elif options.snr:
					print "WARNING: You specified snr but there's no noise to add, apparently!"

			ctfed_projections.append(prj_r)
			#print "Appended ctfed prj in slice j", j
		
			if options.saveprjs and (options.applyctf or options.snr):
				finalprjsED = outname.replace('.hdf', '_ptcl' + str(i).zfill(len(str(nslices))) + '_prjsEDITED.hdf')
				#if options.path + '/' in outname:
				#	finalprjsED = outname.replace('.hdf', '_ptcl' + str(i).zfill(len(str(nslices))) + '_prjsEDITED.hdf') 
				
				finalprjsED = finalprjsED.replace('_preproc','')
				prj_r.write_image( finalprjsED , -1)	
	
	
	
	
		print "The box for IMAGE is with image.get_xsize", image.get_xsize()
		print "Whereas with image['nx']", image['nx']

		box = max(int(image['nx']),int(image['ny']),int(image['nz']))
		
		print "!!!!!!!!!!!!!!!Therefore, box is", box
		
		mode='gauss_2'
		if options.reconstructor:
			if len(options.reconstructor) > 1:
				if 'mode' in options.reconstructor[-1]:
					mode=options.reconstructor[-1]['mode']
					
					print "\nThe reconstructor mode has been changed from default to", mode
					#sys.exit()
					
		r = Reconstructors.get(options.reconstructor[0],{'size':(box,box,box),'sym':'c1','verbose':True,'mode':'gauss_2'})
		#r = Reconstructors.get(options.reconstructor[0],options.reconstructor[1])
		r.setup()
	
		#print "There are these many projections to add to backprojection after all processing", len(ctfed_projections)
	
		k=0
		for p in ctfed_projections:
			#print "Adding projection k", k
			#print "Whose min and max are", p['minimum'], p['maximum']
			#print "The size of the prj to insert is", p['nx']
			p = r.preprocess_slice(p,p['xform.projection'])
			r.insert_slice(p,p['xform.projection'],1.0)
			k+=1
		
		#print "\n\n!!!!!!Will reconstruct the volume for particle i",i
	
		rec = r.finish(True)
		print "\n(e2spt_simulation) I have finished simulating particle number", i
		print "\n"
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

		#print "The apix of rec is", rec['apix_x']
		
		print "\nThe outname to write the particle i", i 
		print "is", outname
		print "\n\n"
		
		#finaloutname = options.path + '/' + outname
		#finaloutname.replace('_preproc','')
		#print "is, finaloutname", finaloutname

		print "rec to return is", rec
			
		#rec.write_image(finaloutname,i)

		#if options.tomogramoutput:
		#	py = random.uniform(0 + image['nx']/2, options.gridholesize - image['nx']/2)									#random distance in Y of the particle's center from the bottom edge in the XY plane, at tilt=0
		#	pz = random.uniform(0 + image['nx']/2, options.icethickness - image['nx']/2) 		
		#	return {'ptcl':rec,'px':px,'py':py,'pz':pz}
		#	
		#else:
		
		if options.finalboxsize:
			box = int(options.finalboxsize)
			print "\nActually, because of finalboxsize$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n options.finalboxsize is", options.finalboxsize
			print "Therefore, box is", box
			
		#print "The boxsize to use is", options.finalboxsize
			
		if finalprjsRAW and options.finalboxsize and float(options.pad) > 1.0:
			print "\n\n\n\n\n\n\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@prjs raw need to be clipped! to a size of", options.finalboxsize
			
			finalprjsRAWclipped = finalprjsRAW.replace('.hdf','_clipped.hdf')  
			clipcmd1 = 'e2proc2d.py ' + finalprjsRAW + ' ' + finalprjsRAWclipped + ' --clip=' + str(options.finalboxsize) + ',' + str(options.finalboxsize) + ' && mv ' + finalprjsRAWclipped + ' ' + finalprjsRAW
			print "with cmmand", clipcmd1
			
			
			p=subprocess.Popen( clipcmd1, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()
			
		if finalprjsED and options.finalboxsize and float(options.pad) > 1.0:
			print "\n\n\n\n\n\n\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@prjs ED need to be clipped! to a size of", options.finalboxsize
			
			finalprjsEDclipped = finalprjsED.replace('.hdf','_clipped.hdf')  
			
			clipcmd2 = 'e2proc2d.py ' + finalprjsED + ' ' + finalprjsEDclipped + ' --clip=' + str(options.finalboxsize) + ',' + str(options.finalboxsize) + ' && mv ' + finalprjsEDclipped + ' ' + finalprjsED 
			print "with cmmand", clipcmd2
			
			
			
			p=subprocess.Popen( clipcmd2, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text=p.communicate()	
			p.stdout.close()
		
		#box = image.get_xsize()
		
		print "\n\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\nThe final boxsize of rec is", rec['nx'], box
		return { classoptions['ptclnum']:rec }


jsonclasses["SubtomoSimTask"]=SubtomoSimTask.from_jsondict


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
				r=etc.get_results(tidsleft[i])					# results for a completed task
				
				if r:
					#print "r is", r
					ptcl=r[0].classoptions["ptclnum"]		# get the particle number from the task rather than trying to work back to it
					results[ptcl] = r[1]						# this will be a list of (qual,Transform)
				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if options.verbose:
			print "  %d tasks, %d complete, %d waiting to start        \r"%(len(tids),ncomplete,nwait)
			sys.stdout.flush()
	
		if len(tidsleft)==0: break
		
	return results



def tomogramsim(options,tomogramdata):
	'''
	Transform gridholesize and icethickness to pixels
	'''

	modelhdr = EMData(options.input,0,True)
	
	if options.gridholesize:
		options.gridholesize = int( options.gridholesize * math.pow(10,6) / model['apix_x'] )

	if options.icethickness:
		options.icethickness = int( options.icethickness * math.pow(10,0) / model['apix_x'] )
	
	tomox=options.gridholesize 						
	tomoy=options.gridholesize
	tomoz=options.icethickness
	print "Gridholesize and icethickness are"
	print options.gridholesize 
	print options.icethickness
	#T=EMData(tomox,tomoy,tomoz)
	print "The size of the tomogam to create is", tomox,tomoy,tomoz

	ptcls=[]
	for w in tomogramdata:
		ptcl=w['ptcl']
		x=w['px']+options.gridholesize/2
		y=w['py']
		z=w['pz']

		r=Region( (x-tomox)/2, (y-tomoy)/2, (z-tomoz)/2, tomox,tomoy,tomoz )
		ptclr=ptcl.get_clip(r)
		ptcls.append(ptclr)
	
	T=sum(ptcls)/len(ptcls)	
	T.process_inplace('math.addnoise',{'noise':1})
	tomogrampath = options.path + '/' + options.tomogramoutput
	T.write_image(tomogrampath,0)
	
	return


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
