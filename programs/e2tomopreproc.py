#!/usr/bin/env python																																																																																																																																																																																																																																																																																																																			#!/usr/bin/python2.7

#====================
#Author: Jesus Galaz-Montoya January/21/2016 , Last update: January/21/2016
#====================
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


from optparse import OptionParser

from EMAN2 import *
from EMAN2jsondb import JSTask,jsonclasses

import sys
from sys import argv

from shutil import copyfile

def main():

	usage = """e2tomopreproc.py <imgs> <options> . 
	This program takes a tiltseries ('.st' or '.ali' file from IMOD) and applies preprocessing operations to them, such as lowpass, highpass, masking, etc.
	The options should be supplied in "--option=value" format, replacing "option" for a valid option name, and "value" for an acceptable value for that option. 
	"""
			
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	
	parser.add_argument("--path",type=str,default='',help="""Directory to store results in. 
		The default is a numbered series of directories containing the prefix 'tomopreproc';
		for example, tomopreproc_02 will be the directory by default if 'tomopreproc_01' 
		already exists.""")
	
	parser.add_pos_argument(name="stack_files",default="",help="Stacks or images to process.")
	
	parser.add_argument("--input",type=str,default='',help=""""tiltseries to process. redundant with --tiltseries, or with providing images as arguments (separated by a space: e2tomopreproc.py stack1.hdf stack2.hdf), but --input takes precedence.""")
	
	parser.add_argument("--tiltseries",type=str,default='',help=""""tiltseries to process. redundant with --input""")

	parser.add_argument("--tltfile",type=str,default='',help="""".tlt file containing the tilt angles for --tiltseries""")
	
	parser.add_argument("--outmode", type=str, default='', help="""All EMAN2 programs write images with 4-byte floating point values when possible by default. This allows specifying an alternate format when supported: float, int8, int16, int32, uint8, uint16, uint32. Values are rescaled to fill MIN-MAX range.""")
	
	parser.add_argument("--dontcleanup", action='store_true', default=False, help="""If specified, intermediate files will be kept.""")
	
	parser.add_argument("--clip",type=str,default='',help="""Default=None. This resizes the 2-D images in the tilt series. If one number is provided, then x and y dimensions will be made the same. To specify both dimensions, supply two numbers, --clip=x,y. Clipping will be about the center of the image.""")
			
	#parser.add_argument("--apix",type=float,default=0.0,help="""True apix of images to be written on final stack.""")
	
	parser.add_argument("--shrink", type=float,default=0.0,help="""Default=0.0 (no shrinking). Can use decimal numbers, larger than 1.0. Optionally shrink the images by this factor. Uses processor math.fft.resample.""")
		
	parser.add_argument("--threshold",type=str,default='',help="""Default=None. A threshold processor applied to each image.""")

	parser.add_argument("--erasegold",action='store_true',default='',help="""Default=False. Runs erase_gold.py on the stack.""")
	
	parser.add_argument("--mask",type=str,default='', help="""Default=None. Masking processor applied to each image.""")
	
	parser.add_argument("--maskbyangle",action='store_true',default=False,help="""Default=False. Requires --tltfile. This will mask out from tilted images the info that isn't present at the 0 tilt angle. It uses the tomo.tiltedgemask processor (type 'e2help.py processors' at the commandline to read a description of the processor and its parameters). Provide --maskbyanglefalloff and --maskbyanglesigma to modify the default parameters.""")
	
	parser.add_argument("--maskbyanglefalloff", type=int, default=4,help="""Default=4. Number of pixels over which --maskbyangle will fall off to zero.""")
	
	parser.add_argument("--maskbyanglesigma", type=float, default=2.0,help="""Default=2.0. Number of sigmas for the width of the gaussian fall off in --maskbyangle and --maskbyanglefalloff""")
	
	parser.add_argument("--normproc",type=str, default='',help="""Default=None (see 'e2help.py processors -v 10' at the command line). Normalization processor applied to each image.""")
	
	parser.add_argument("--normalizeimod",action='store_true',default=False,help="""Default=False. This will apply 'newstack -float 2' to the input stack. Requires IMOD.""")
	
	parser.add_argument("--preprocess",type=str,default='',help="""Any processor (see 'e2help.py processors -v 10' at the command line) to be applied to each image.""")
	
	parser.add_argument("--lowpassfrac",type=float,default=0.0,help="""Default=0.0 (not used). Fraction of Nyquist to lowpass at. The processor used is filter.lowpass.tanh""")
	
	parser.add_argument("--highpasspix",type=int,default=0,help="""Default=0 (not used). Number of Fourier pixels to apply highpass filter at. The processor used is filter.highpass.gauss.""")
	
	parser.add_argument("--parallel",type=str, default="thread:1", help="""default=thread:1. Parallelism. See http://blake.bcm.edu/emanwiki/EMAN2/Parallel""")
	
	parser.add_argument("--prenadminite",type=int, default=0, help="""Default=0. Requires IMOD to be installed. Used to apply prenad filtering to a tiltseries. This is the --minite parameter in IMOD's preNAD program (minimum number of iterations).""")
	
	parser.add_argument("--prenadmaxite",type=int, default=0, help="""Default=0. Requires IMOD to be installed. Used to apply prenad filtering to a tiltseries. This is the --maxite parameter in IMOD's preNAD program (maximum number of iterations).""")
	
	parser.add_argument("--prenadsigma",type=int, default=0, help="""Default=0. Requires IMOD to be installed. Used to apply prenad filtering to a tiltseries. This is the --sigma parameter in IMOD's preNAD program (initial sigma for 'smoothing structure tensor').""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")
	
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()	
	

	logger = E2init(sys.argv, options.ppid)
	print "\n(e2tomopreproc)(main) started log"	
	
	from e2spt_classaverage import sptmakepath
	
	options = sptmakepath(options,'tomopreproc')
	
	#print "args are",args

	infiles = []
	if not options.input:
		#try:
		#	infiles.append( sys.argv[1] )
		#except:
		if options.tiltseries:
			infiles.append( options.tiltseries )
		else:
			if args:
				print "copying args to infiles"
				infiles = list(args)
				print "infiles are", infiles
			else:
				print "\n(e2tomopreproc)(main) ERROR: must provide input files as arguments or via the --input or --tiltseries parameters."


	if infiles:
		print "\n(e2tomopreproc)(main) identified --input", options.input
		#print " .ali in options.input[:-4]", '.ali' in options.input[-4:]
		#print "options.input[-4] is", options.input[-4:]
		
		for infile in infiles:
			if '.ali' in infile[-4:] or '.st' in infile[-3:] or '.mrc' in infile[-4:] or '.mrcs' in infile[-5:] or '.hdf' in infile[-4:]:
				pass
			else:
				print "\n(e2tomopreproc)(main) ERROR: invalid image extension %s for image %s. Extension must be .st, .ali, .hdf, .mrc or .mrcs" %(options.input.split('.')[-1], infile)
				sys.exit(1)
	else:
		print "\n(e2tomopreproc)(main) ERROR: no images found/provided"
		sys.exit(1)
		
	originalextension = infiles[0].split('.')[-1]
	
	angles = {}
	if options.maskbyangle or (options.prenadminite and options.prenadmaxite and options.prenadsigma):
	
		if not options.tltfile:
			print "\n(e2tomopreproc)(main) ERROR: --maskbyangle and --prenad parameters require --tltfile"
			sys.exit(1)
		
		else:
			f = open( options.tltfile, 'r' )
			lines = f.readlines()
			print "\nnumber of lines read from --tltfile", len(lines)
			f.close()
			#print "lines in tlt file are", lines
			k=0
			for line in lines:
				line = line.replace('\t','').replace('\n','')
	
				if line:
					angle = float(line)
					angles.update( { k:angle } )
					if options.verbose:
						print "appending angle", angle
					k+=1
			if len(angles) < 2:
				print "\nERROR: something went terribly wrong with parsing the --tltlfile. This program does not work on single images"
				sys.exit()

		if len(angles) < 2:
			print "\nERROR: (second angle check) something went terribly wrong with parsing the --tltlfile. This program does not work on single images"
			sys.exit()
				
	
	
	
	
	print "\n(e2spt_preproc)(main) - INITIALIZING PARALLELISM!\n"

	from EMAN2PAR import EMTaskCustomer
	etc=EMTaskCustomer(options.parallel)
	pclist=[options.input]

	etc.precache(pclist)
	print "\n(e2spt_preproc)(main) - precaching --input"

	tasks=[]
	results=[]
	
	mrcstacks = []
	print "there are these many infiles to loop over", len(infiles)



	if options.lowpassfrac:
		hdr = EMData( infiles[0], 0, True )
		apix = hdr['apix_x']
		print "\n(e2spt_preproc)(main) apix is",apix
		nyquist = 2.0 * apix
		print "\n(e2spt_preproc)(main) therefore nyquist resolution is", nyquist
		print
		lowpassres = nyquist/options.lowpassfrac
		
		options.lowpassfrac = 1.0/(lowpassres)
		if float(options.shrink) > 1.0:
			options.lowpassfrac /= float(options.shrink)
			
			print "there's shrinking", options.shrink
			lowpassres = nyquist/options.lowpassfrac

		print "\n(e2spt_preproc)(main) and final lowpass frequency is", options.lowpassfrac

		print "corresponding to lowpassres of",lowpassres

	for infile in infiles:
	
		mrcstack = options.path + '/' + infile
		print "infile is", infile
		print "infile[-5:] is ", infile[-5:]
		if '.hdf' in infile[-5:]:
			print "replacing .hdf extension"
			mrcstack = options.path + '/' + infile.replace('.hdf','.mrc')
	
		if '.mrcs' in infile[-5:]:
			print "replacing .mrcs extension"
			mrcstack = options.path + '/' + infile.replace('.mrcs','.mrc')
	
		if '.st' in infile[-5:]:
			print "replacing .st extension"
			mrcstack = options.path + '/' + infile.replace('.st','.mrc')	

		if '.ali' in infile[-5:]:
			print "replacing .ali extension"
			mrcstack = options.path + '/' + infile.replace('.ali','.mrc')
			
		if '.tif' in infile[-5:]:
			print "replacing .ali extension"
			mrcstack = options.path + '/' + infile.replace('.tif','.mrc')
	
		#go = 0
		#if go:
		print "mrcstack is",mrcstack
		
		#outname = outname.replace('.mrc','.mrcs')
	
		mrcstacks.append( mrcstack )
		
		go = 0
		if options.maskbyangle:
			outname = mrcstack.replace('.mrc','_UNSTACKED.mrc')
			print "therefore, outname is", outname
	
			cmd = 'e2proc2d.py ' + infile + ' ' + outname + ' --unstacking --threed2twod'

			#from shutil import copyfile
			#copyfile(options.input, outname)
			#print "copied input to", outname

			if options.outmode:
				cmd += ' --outmode=' + options.outmode

			if options.verbose:
				cmd += ' --verbose=' + str(options.verbose)
				print "\ncommand to unstack original input tiltseries is", cmd	

			print "\n(e2tomopreproc)(main) unstacking command is", cmd

			p = subprocess.Popen( cmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			#p = subprocess.Popen( cmd , shell=True, stdout=subprocess.PIPE)

			text = p.communicate()	
			#p.stdout.close()

			p.wait()
		
			if p.returncode == 0:
				go = 1
		else:
			go = 1
	
		
		if go:

			imgs = []
			if options.maskbyangle:
				c = os.getcwd() + '/' + options.path 
				findir = os.listdir( os.getcwd() + '/' + options.path )

				print "\n(e2tomopreproc)(main) directory to look for images is", c	
				for f in findir:
					#if '.mrcs' in f:
					if "_UNSTACKED" in f:
						imgs.append( options.path + '/' +f )

				kk=0
				imgs.sort()
				print "\n(e2spt_preproc)(main) found these many images", len( imgs )		

				for img in imgs:
					#task=None

					#if options.maskbyangle:
					outimage = img.replace('.mrc','_preproc.mrc')
					task = TomoPreproc2DTask( img, options, angles[kk], outimage )
					tasks.append(task)
					kk+=1
			else:
				outimage = options.path + '/' + infile.replace('.mrc','_preproc.mrcs')
				task = TomoPreproc2DTask( infile, options, 0, outimage )
				tasks.append(task)
				
					
			#else:
			#	newmrcs = mrcstack.replace('.mrc','.mrcs')
			#	print "copying file %s to %s" %(infile,newmrcs)
			#	copyfile( infile, newmrcs  )
			#	imgs.append( newmrcs )
			
			

			

				
				#print "and the final lowpass frequency will be", options.lowpassfrac

			

			
	tids = etc.send_tasks(tasks)
	if options.verbose: 
		print "\n(e2spt_preproc)(main) preprocessing %d tasks queued" % (len(tids)) 

	results = get_results( etc, tids, options )

	print "\n(e2tomopreproc)(main) these many images have been processsed",len(results)

	
	imgspreproc = []
	findir = os.listdir( os.getcwd() + '/' + options.path )
	
	#for mrcstack in mrcstacks:


	for f in findir:
		if "_preproc.mrc" in f:
			print "found preprocessed image", f
			imgspreproc.append( options.path + '/' + f )
		else:
			print "this file is NOT a preprocessed image", f

	imgspreproc.sort()

	print "\n(e2tomopreproc)(main) these many preprocessed images loaded", len(imgspreproc)
	
	finalfiles=[]
	
	if options.maskbyangle:
		
		outfile = mrcstack.replace('.mrc','.mrcs')
		print "for RESTACKING"
		print "\n\n\noutfile is", outfile

		for f in imgspreproc:
			print "appending image %s to outfile %s" %(f,outfile)			
			cmd = 'e2proc2d.py ' + f + ' ' + outfile
			if options.outmode:
				cmd += ' --outmode=' + options.outmode

			if options.verbose:
				cmd += ' --verbose ' + str(options.verbose)

			print "\ncmd is with .mrcs outputformat is", cmd
			print "becauase outfile is",outfile	
			p = subprocess.Popen( cmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text = p.communicate()	
			p.stdout.close()		
	
		finaloutput = outfile.replace('.mrcs', '.' + originalextension)
		os.rename( outfile, finaloutput )
		
		finalfiles.append( finaloutput )
	else:
		finalfiles = list( imgspreproc )
	
	
	for finalf in finalfiles:
		if not options.tltfile:
			break
	
		if options.normalizeimod:
			try:
				cmd = 'newstack ' + finalf + ' ' + finalf + ' --float 2'
				print "normalizeimod cmd is", cmd
				p = subprocess.Popen( cmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				text = p.communicate()	
				p.wait()
			except:
				print "\nERROR: --normalizeimod skipped. Doesn't seem like IMOD is installed on this machine"		

		if not options.dontcleanup and options.maskbyangle:
			purge( options.path, '_preproc.mrc')
			purge( options.path, '_UNSTACKED')	
			purge( options.path, '~')
		
		if options.tltfile:
			if options.prenadminite or options.prenadmaxite or options.prenadsigma:

				if options.prenadminite and options.prenadmaxite and options.prenadsigma:
					cmd = 'preNAD -input ' + finalf + ' -output ' + finalf.replace('.'+originalextension, '_prenad.' + originalextension) + ' -minite ' + str(options.prenadminite) + ' -maxite ' + str(options.prenadmaxite) + ' -sigma ' + str(options.prenadsigma) + ' -angles ' + options.tltfile 
					if options.verbose:
						print "\n(e2tomopreproc)(main) prenad cmd to run is", cmd
					try:
						p = subprocess.Popen( cmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
						text = p.communicate()	
						p.wait()
					except:
						print "\nERROR: check that a version of IMOD containing the preNAD program is correctly installed on this machine"

				else:
					if options.prenadminite:
						if not options.prenadmaxite:
							print "\nERROR: --prenadmaxite required with --prenadminite"
						if not options.prenadsigma:
							print "\nERROR: --prenadsigma required with --prenadminite"

					if options.prenadmaxite:
						if not options.prenadminite:
							print "\nERROR: --prenadminite required with --prenadmaxite"
						if not options.prenadsigma:
							print "\nERROR: --prenadsigma required with --prenadmaxite"

					if options.prenadsigma:
						if not options.prenadminite:
							print "\nERROR: --prenadminite required with --prenadsigma"
						if not options.prenadmaxite:
							print "\nERROR: --prenadmaxite required with --prenadsigma"
					
		
	E2end(logger)	
	return()
	

def purge( dir, stem ):
	import os, re
	
	for f in os.listdir( dir ):
		if re.search( stem, f ):
			os.remove( os.path.join( dir, f) )
	return


'''
CLASS TO PARALLELIZE PREPROCESSING STEPS
'''
class TomoPreproc2DTask(JSTask):
	'''This is a task object for the parallelism system. It is responsible for preprocessing a 2-D image in a tiltseries, with a variety of options'''


	def __init__(self, image, options, angle, outimage):
	
		data={"image":image}
		
		JSTask.__init__(self,"TomoPreproc2d",data,{},"")

		self.classoptions={"options":options, "angle":angle,"outimage":outimage }
	
	def execute(self,callback=None):
		
		options = self.classoptions['options']
		angle = self.classoptions['angle']
		outimage = self.classoptions['outimage']
		
		image = self.data["image"]
		
		if options.verbose:
			print "\n(e2tomopreproc)(class), image is", image
		
		hdr=EMData(image,0,True)
		nx=hdr['nx']
		ny=hdr['ny']
		nz=hdr['nz']
		print "original nx,ny,nz are",nx,ny,nz
		
		#cmd = 'e2proc2d.py ' + image + ' ' + image.replace('.mrc','_preproc.mrc')

		#if nz < 2:
		#cmd = 'e2proc2d.py ' + image + ' ' + image.replace('.mrcs','_preproc.mrcs')
		#
		cmd = 'e2proc2d.py ' + image + ' ' + outimage

		if '.mrcs' in outimage:
			cmd += ' --threed2twod'
		#if nz > 1:
		#		cmd += ' --threed2threed'
		#print "cmd is", cmd
		
		
		if options.outmode:
			print "adding outmode"
			cmd += ' --outmode ' + options.outmode
	
		if options.normproc:
			print "adding normproc"
			cmd += ' --process ' + options.normproc
									
		if options.threshold:
			print "adding threshold"
			cmd += ' --process ' + options.threshold	
		
		if options.lowpassfrac:
			print "adding lowpassfrac"
			cmd += ' --process filter.lowpass.tanh:cutoff_freq=' + str(options.lowpassfrac)	

		if options.highpasspix:
			print "adding highpasspix"
			cmd += ' --process filter.highpass.gauss:cutoff_pixels=' + str(options.highpasspix)	
		
		if options.preprocess:
			print "adding preprocess"
			cmd += ' --process ' + options.preprocess		
		
		if options.maskbyangle:	
			print "adding maskbyangle"
			cmd += ' --process tomo.tiltedgemask:angle=' + str(angle) + ':gauss_falloff=' + str(options.maskbyanglefalloff) + ':gauss_sigma=' + str(options.maskbyanglesigma)
	
		if options.mask:
			print "adding mask"
			cmd += ' --process ' + options.mask	
		
		if options.clip:
			print "adding clip"
			clips = options.clip.split(',')
			print "clips", clips
		
			clipx = clips[0]
			cmd += ' --clip ' + clipx
			if len(	clips ) > 1:
				clipy = clips[1]
				cmd += ',' + clipy
			else:
				cmd += ',' + clipx
				
			
		if options.shrink:
			print "adding shrink"
			cmd += ' --process math.fft.resample:n=' + str(options.shrink)
		
		if options.verbose:
			print "\n(e2tomopreproc)(class) cmd", cmd
		
		os.system( cmd )
		
		#p = subprocess.Popen( cmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		#text = p.communicate()	
		#p.stdout.close()

		
		return


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
					#ptcl=r[0].classoptions["ptclnum"]		#Get the particle number from the task rather than trying to work back to it
					#results[ptcl] = r[1]
					results[i] = 1					
				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if options.verbose:
			print ("\n%d tasks, %d complete, %d waiting to start \r" % (len(tids),ncomplete,nwait))
			sys.stdout.flush()
	
		if len(tidsleft)==0: break
		
	return results


	
if __name__ == "__main__":
    main()
    