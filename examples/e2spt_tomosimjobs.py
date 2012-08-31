#!/usr/bin/env python

#
# Author: Jesus Galaz, 02/21/2012 - using code and concepts drawn from Jesus Galaz's scripts
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
	
	parser.add_argument("--snrlowerlimit", type=int,default=0,help="Minimum weight for noise compared to singal.")
	parser.add_argument("--snrupperlimit", type=int,default=1,help="Maximum weight for noise compared to singal.")

	parser.add_argument("--snrchange", type=int,default=1,help="""Step to vary snr from one run to another. 
										For example, if this parameter is set to 2, snr will be tested at 1,3,5... up to --snrupperlimit.""")

	parser.add_argument("--tiltrangelowerlimit", type=int,default=60,help="""Minimum value for imaging range (at a value of 90, there's no missing wedge. 
											60 would mean the data will be simulated, as if it came from a tomogram that was reconstructed from a til series
											collected from -60deg to 60deg)""")
	parser.add_argument("--tiltrangeupperlimit", type=int,default=61,help="""Maximum value to simulate the imaging range. Simulations will be run starting at --tiltrangelowerlimit, and will increase
											by --tiltrangestep, until the last simulation is done at --tiltrangeupperlimit.""")
	parser.add_argument("--tiltrangechange", type=int,default=1,help="Amount (in degrees) to decrease the size of the missing wedge from one run to another.")	
	
	
	parser.add_argument("--tiltsteplowerlimit", type=int,default=1,help="""Within each tiltrange simulated, you can simulate individual pictures taken with different tilt steps.
										For example, if you collect images from -60deg to 60deg with a 2deg tilt step, the tilt series will have 61 images.
										If, on the other hand, the tilt step was 4deg, the tilt series will only have 31 images.
										--tiltstepupperlimit is the largest step size you want to simulate.""")
	parser.add_argument("--tiltstepupperlimit", type=int,default=2,help="""Within each tiltrange simulated, you can simulate individual pictures taken with different tilt steps.
										For example, if you collect images from -60deg to 60deg with a 2deg tilt step, the tilt series will have 61 images.
										If, on the other hand, the tilt step was 4deg, the tilt series will only have 31 images.
										--tiltstepupperlimit is the largest step size you want to simulate.""")
	parser.add_argument("--tiltstepchange", type=int,default=1,help="""Increase in size of tilt step from one run to another. 
											Jobs will be run using --tiltstepstep as the first value, and then adding that value on subsequent runs until --tiltsteplimit is reached""")	
	
	"""
	Parameters to be passed on to e2spt_simulation.py
	"""
	parser.add_argument("--input", type=str, help="""The name of the input volume from which simulated subtomograms will be generated. 
							The output will be in HDF format, since volume stack support is required. The input CAN be PDB, MRC or and HDF stack. 
							If the input file is PDB or MRC, a version of the supplied model will be written out in HDF format.
							If the input file is a stack, simulatd subvolumes will be generated from each model in the stack and written to different output stacks.
							For example, if the input file contains models A and B, two output stacks with simulated subvolumes will be generated.""", default=None)
				
	parser.add_argument("--filter",type=str,help="""A filter (as in a processor from e2proc3d.py) apply to the model before generating simulated particles from it.
							Type 'e2help.py processors' at the command line and find the options availbale from the processors list)""",default=None)
	
	parser.add_argument("--shrink", type=int,default=1,help="Optionally shrink the input volume before the simulation if you want binned/down-sampled subtomograms.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_argument("--nptcls", type=int,default=10,help="Number of simulated subtomograms tu generate per referece.")
	parser.add_argument("--tx", type=int,default=None,help="""Maximum number of pixels to randomly translate each subtomogram in X. The random translation will be picked between -txrange and +txrange. 
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--ty", type=int,default=None,help="""Maximum number of pixels to randomly translate each subtomogram in Y. The random translation will be picked between -tyrange and +tyrange.
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--tz", type=int,default=None,help="""Maximum number of pixels to randomly translate each subtomogram in Z. The random translation will be picked between -tzrange and +tzrange.
								     Default value is set by --transrange, but --txrange will overwrite it if specified.""")
	parser.add_argument("--transrange", type=int,default=4,help="""Maximum number of pixels to randomly translate each subtomogram in all X, Y and Z. 
									The random translation will be picked between -transrage and +transrange; --txrange, --tyrange and --tzrange overwrite --transrange for each specified direction.""")
	
	parser.add_argument("--applyctf", action="store_true",default=False,help="If on, it applies ctf to the projections in the simulated tilt series based on defocus, cs, and voltage parameters.")

	parser.add_argument("--defocus", type=int,default=3,help="Intended defocus at the tilt axis (in microns) for the simulated tilt series.")
	parser.add_argument("--voltage", type=int,default=200,help="Voltage of the microscope, used to simulate the ctf added to the subtomograms.")
	parser.add_argument("--cs", type=int,default=2,help="Cs of the microscope, used to simulate the ctf added to the subtomograms.")

	parser.add_argument("--gridholesize", type=int,default=2,help="""Size of the carbon hole for the simulated grid (this will determine shifts in defocus for each particle at 
									each tilt step, depending on the position of the particle respect to the tilt axis, which is assigned randomly.""")
	parser.add_argument("--saverandstack", action="store_true",default=False,help="Save the stack of randomly oriented particles, before subtomogram simulation (before the missing wedge and noise are added).")
	parser.add_argument("--saveprjs", action="store_true",default=False,help="Save the projections (the 'tilt series') for each simulated subtomogram.")

	parser.add_argument("--reconstructor", type=str,default="fourier",help="""The reconstructor to use to reconstruct the tilt series into a tomogram. Type 'e2help.py reconstructors' at the command line
											to see all options and parameters available.""")

	parser.add_argument("--pad", type=int,default=0,help="""If on, it will increase the box size of the model BEFORE generating projections and doing 3D reconstruction of simulated sutomograms.""")								
	
	parser.add_argument("--finalboxsize", type=int,default=0,help="""The final box size to clip the subtomograms to.""")								

	parser.add_argument("--snr",type=int,help="Weighing noise factor for noise added to the image. Only words if --addnoise is on.",default=5)
	parser.add_argument("--addnoise",action="store_true",default=False,help="If on, it adds random noise to the particles")
	
	parser.add_argument("--sym",type=str,default='c1',help="If your particle is symmetrical, you should randomize it's orientation withing the asymmetric unit only. Thus, provide the symmetry.")

	parser.add_argument("--notrandomize",action="store_true",default=False,help="This will prevent the simulated particles from being rotated and translated into random orientations.")
	parser.add_argument("--simref",action="store_true",default=False,help="This will make a simulated particle in the same orientation as the original input (or reference).")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--negativecontrast",action="store_true",default=False,help="This will make the simulated particles be like real EM data before contrast reversal. Otherwise, 'white protein' (positive density values) will be used.")

	parser.add_argument("--testalignment",action="store_true",default=False,help="This will run e2spt_classaverage.py to test the alignment of the particles against the simulated reference.")

	parser.add_argument("--quicktest",action="store_true",default=False,help="This will run e2spt_classaverage.py with minimal parameters to quickly test the program.")

	(options, args) = parser.parse_args()	
	
	logger = E2init(sys.argv, options.ppid)
	
	snrl = options.snrlowerlimit
	snru = options.snrupperlimit
	snrch = options.snrchange

	tiltrangel = options.tiltrangelowerlimit
	tiltrangeu = options.tiltrangeupperlimit
	tiltrangech = options.tiltrangechange

	tiltstepl = options.tiltsteplowerlimit
	tiltstepu = options.tiltstepupperlimit
	tiltstepch = options.tiltstepchange
	
	if tiltstepl == 0:
		print "ERROR! You cannot start with a tilt step of 0. The minimum tiltstep is 1, thus, the lower limit for this parameter, tiltsteplowerlimit, must be at least 1."
	
	'''
	Make the directory where to create the database where the results will be stored
	'''
	
	#if options.path and ("/" in options.path or "#" in options.path) :
	#	print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
	#	sys.exit(1)
		
	#if options.path and options.path[:4].lower()!="bdb:": 
	#	options.path="bdb:"+options.path

	#if not options.path: 
	#	#options.path="bdb:"+numbered_path("sptavsa",True)
	#	options.path = "sptsim_01"
	
	
	if options.path and ("/" in options.path or "#" in options.path) :
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. "
		sys.exit(1)

	if not options.path: 
		#options.path="bdb:"+numbered_path("sptavsa",True)
		options.path = "sptsim_01"
	
	files=os.listdir(os.getcwd())
	print "right before while loop"
	while options.path in files:
		print "in while loop, options.path is", options.path
		#path = options.path
		if '_' not in options.path:
			print "I will add the number"
			options.path = options.path + '_00'
		else:
			jobtag=''
			components=options.path.split('_')
			if components[-1].isdigit():
				components[-1] = str(int(components[-1])+1).zfill(2)
			else:
				components.append('00')
						
			options.path = '_'.join(components)
			#options.path = path
			print "The new options.path is", options.path

	if options.path not in files:
		
		print "I will make the path", options.path
		os.system('mkdir ' + options.path)
	
	#if options.testalignment:
	#	resultsdir = 'results_ali_errors'
	#	os.system('cd ' + options.path + ' && mkdir ' + resultsdir)
	
	nrefs = EMUtil.get_image_count(options.input)
		
	tiltrange=tiltrangel
	while tiltrange <tiltrangeu:
		print "tiltrage is", tiltrange
		
		tiltstep=tiltstepl		
		while tiltstep < tiltstepu:
			print "Tilt step is", tiltstep
				
			snr=snrl
			while snr < snru:
				print "Snr is", snr
								
				rootpath = os.getcwd()

				for d in range(nrefs):
					
					modeltag = ''
					subpath = 'TR' + str(tiltrange).zfill(2) + '_TS' + str(tiltstep).zfill(2) + '_SNR' + str(snr).zfill(2)
					
					inputdata = options.input
					
					if nrefs > 1:
						modeltag = 'model' + str(d).zfill(2)
						subpath += '_' + modeltag
						
						model = EMData(options.input,d)
						newname = options.path + '/' + inputdata.split('/')[-1].replace('.hdf','_' + modeltag + '.hdf')
						model.write_image(newname,0)
						
						inputdata = newname.split('/')[-1]
					
					subtomos =  subpath + '.hdf'

					jobcmd = 'e2spt_simulation.py --input=' + inputdata + ' --output=' + subtomos + ' --snr=' + str(snr) + ' --nptcls=' + str(options.nptcls) + ' --tiltstep=' + str(tiltstep) + ' --tiltrange=' + str(tiltrange) + ' --transrange=' + str(options.transrange) + ' --pad=' + str(options.pad) + ' --shrink=' + str(options.shrink) + ' --finalboxsize=' + str(options.finalboxsize)

					if options.simref:
						jobcmd += ' --simref'
					if options.addnoise:
						jobcmd += ' --addnoise'
					if options.saveprjs:
						jobcmd += ' --saveprjs'
					if options.negativecontrast:
						jobcmd += ' --negativecontrast'

					jobcmd += ' --path=' + subpath				

					cmd = 'cd ' + options.path + ' && ' + jobcmd

					resultsfiles=[]

					if options.testalignment:

						#modeldir = ''
						#if nrefs > 1:
						#	modeldir = '/model' + str(d).zfill(2)

						cmd = cmd + ' && cd ' + rootpath + '/' + options.path + '/' + subpath

						#subtomos = options.input.split('/')[-1].replace('.hdf','_sptsimMODEL_randst_n' + str(options.nptcls) + '_' + subpath + '_subtomos.hdf')

						print "\n\nSubtomos name will be\n", subtomos

						ref = inputdata.split('/')[-1].replace('.hdf','_sptsimMODEL_SIM.hdf')

						#if nrefs > 1:
						#	modd = 'model' + str(d).zfill(2)
						#	ref = options.input.split('/')[-1].replace('.hdf','_sptsimMODELS_' + modd + '.hdf')

						output=subtomos.replace('.hdf', '_avg.hdf')
						print "\n\n$$$$$$$$$$$$$$$$$$$$$$\nRef name is\n$$$$$$$$$$$$$$$$$$$\n", ref

						print "\n\noutput name is\n", output

						alipath=output.replace('_avg.hdf','_ali')
						print "\n\nAlipath for results will be\n", alipath

						alicmd = " && e2spt_classaverage.py --path=" + alipath + " --input=" + subtomos.replace('.hdf','_ptcls.hdf') + " --output=" + output + " --ref=" + ref + " --npeakstorefine=4 -v 0 --mask=mask.sharp:outer_radius=-4 --lowpass=filter.lowpass.gauss:cutoff_freq=.02 --align=rotate_translate_3d:search=" + str(options.transrange) + ":delta=12:dphi=12:verbose=0 --parallel=thread:8 --ralign=refine_3d_grid:delta=12:range=12:search=2 --averager=mean.tomo --aligncmp=ccc.tomo --raligncmp=ccc.tomo --shrink=2 --savesteps --saveali --normproc=normalize.mask"

						if options.quicktest:
							alicmd = " && e2spt_classaverage.py --path=" + alipath + " --input=" + subtomos.replace('.hdf','_ptcls.hdf') + " --output=" + output + " --ref=" + ref + " -v 0 --mask=mask.sharp:outer_radius=-4 --lowpass=filter.lowpass.gauss:cutoff_freq=.02 --align=rotate_symmetry_3d:sym=c1:verbose=0 --parallel=thread:8 --ralign=None --averager=mean.tomo --aligncmp=ccc.tomo --raligncmp=ccc.tomo --shrink=3 --savesteps --saveali --normproc=normalize.mask"

						aliptcls = output.replace('_avg.hdf','_ptcls_ali.hdf')

						print "\n\aliptcls name is\n", aliptcls

						extractcmd = " && cd " + alipath + " && e2proc3d.py bdb:class_ptcl " + aliptcls

						resultsfile=aliptcls.replace('_ptcls_ali.hdf','_ali_error.txt')

						solutioncmd = " && e2spt_transformdistance.py --input=" + aliptcls + ' --output=' + resultsfile

						rfilecmd =  ' && mv ' + resultsfile + ' ' +  rootpath + '/' + options.path

						cmd = cmd + alicmd + extractcmd + solutioncmd + rfilecmd

					print "The command to execute is \n\n %s \n\n" %(cmd)

					os.system(cmd)

				snr += snrch
				
			tiltstep += tiltstepch
		
		tiltrange += tiltrangech
	
	resultsdir = 'results_ali_error'
	if nrefs >1:
		for i in range(nrefs):
		
			modname = 'model' + str(i).zfill(2)
			print "\nI will make this moddir", modname
			cmdd='cd ' + options.path + ' && mkdir ' + modname + ' && mv *' + modname + '* ' + modname
			print "\nBy executing this command", cmdd
			os.system(cmdd)

		for i in range(nrefs):
			modname = 'model' + str(i).zfill(2)
			
			resultsdir = rootpath + '/' + options.path + '/' + modname + '/results_ali_error_' + modname
			print "\n\n\n\n*******************\nResults dir is\n", resultsdir
			
			os.system('mkdir ' +  resultsdir + ' && cd ' + options.path + '/' + modname + ' && mv *error.txt ' + resultsdir)
	
	else:
		os.system('cd ' + options.path + ' && mkdir ' + resultsdir + ' && mv *error.txt ' + resultsdir)
			 	
	E2end(logger)

	return()

if __name__ == "__main__":
    main()
