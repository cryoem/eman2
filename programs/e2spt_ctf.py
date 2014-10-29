#!/usr/bin/python2.7

#
# Author: Jesus Galaz, 11/01/2012; last update 26/oct/2014
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

import os
from EMAN2 import *
#from time import time


#import matplotlib
#matplotlib.use('Agg',warn=False)		 

#import matplotlib.pyplot as plt

import sys
import numpy		 
import math
	 
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """
		WARNING: Experimental program. Processing subtiltseries is enabled (Oct/2014),
		but other functionality might be incomplete.
		
		This program can do the following:
		1) Flip phases of the entire image for all images in a tilt series using the same 
		parameters for all images (approximate solution that will only work for low-tilt
		and/or small images) or
		2) Using different parameters for each image (supply --ctfparamfile or
		--defocilist in that case)
		
		The same can be done for subtiltseries extracted with e2spt_subtilt.py,
		except that for subtiltseries the correction is tweaked on a per-particle basis and
		is thus as accurate as possible.
		
		3) The program can also correct frames (in tilt series or subtiltseries) on a 
		strip-by-strip basis."""
		
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--tiltseries", type=str, default='', help="""Aligned tilt series.
		File format must be MRC and must end in .mrc or .st or .ali""")
	
	parser.add_argument("--imagestem",type=str,default='',help="""If the images to apply ctf
		correction on are already unstacked and are individual mrc files, supply a common
		string to all of them.""")
	
	parser.add_argument("--subtiltsdir",type=str,default='',help="""Provide a directory
		containing individual stacks, where each stack is a 'mini tilt series' or a 'subtilt series'
		for single particles. Then, each image for each particle in the dir will be
		phase-phlipped using the ctf parameters you provide.
		If each image in the subtilt series is at a different defocus, then the parameters
		should be provided through --ctfparamsfile, whith a different defocus value per row.
		(There should be as many rows as images in each subtiltseries).""")
		
	parser.add_argument("--infodir",type=str,default='',help="""Folder typically produced
		by e2evalimage.py containing info.json files, one per tilt image in a tilt series.
		Each .json file should contain the fitted ctf and all associated parameters for each tilt image.""")
	
	parser.add_argument("--output", type=str, default='',help="""Name for the tilt series saved as an 
		.hdf stack; also, this name will be used as the stem for all other files produced.""")

	parser.add_argument("--path",type=str,default='sptctf',help="""Directory to store results in. 
		The default is a numbered series of directories containing the prefix 'sptctf'; 
		for example, sptctf_02 will be the directory by default if 'sptctf_01' already exists.""")
	
	parser.add_argument("--reconstructor", type=str,default="fourier",help="""The reconstructor 
		to use to reconstruct the tilt series into a tomogram. Type 'e2help.py reconstructors' 
		at the command line to see all options and parameters available.
		To specify the interpolation scheme for the fourier reconstruction, specify 'mode'.
		Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 
		'gauss_5_slow', 'gypergeom_5', 'experimental'.
		For example --reconstructor=fourier:mode=gauss_5 """)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	#parser.add_argument("--path",type=str,default='',help="""Directory to store results in. 
	#	The default is a numbered series of directories containing the prefix 'sptctf'; 
	#	for example, sptctf_02 will be the directory by default if 'sptctf_01' already exists.""")
		
	#parser.add_argument("--tilesize",type=int,default=512,help="""This option will divide each image
	#	in the tilt series into squared tiles of the specified side length, and all the tiles 
	#	in a strip parallel to the tilt axis will be grouped in a single .hdf stack""")
		
	#parser.add_argument("--dontunstack",action='store_true',default=False,help=""""Supply this
	#	parameter when the tilt series is already unstacked. In this case, you should supply
	#	a common string to all the images in the tilt series via --images""")
	
	parser.add_argument("--pad2d", type=float,default=0.0,help="""Padding factor to zero-pad
		the 2d images in the tilt series prior to reconstruction.
		(The final reconstructed subvolumes will be cropped to the original size).""")

	parser.add_argument("--pad3d", type=float,default=0.0,help="""Padding factor to zero-pad
		the reconstruction volume. (The final reconstructed subvolumes will be cropped to 
		the original size).""")	
		
	parser.add_argument("--save3d",action='store_true',default=False,help="""If on, the CTF
		corrected subtiltseries will be reconstrcuted into subvolumes.
		Options --reconstructor, --pad2d, --pad3d are used if --save3d is on.""")
	
	parser.add_argument("--icethickness", type=int,default=0,help="""This corresponds
		to the Z dimension in pixels of the reconstructed raw tomogram (uncropped), at the same binning
		(sampling) as the provided tiltseries, images or subtiltseries.
		This value MUST be provided, except if --subtiltsdir is given.
		""")
	
	parser.add_argument("--icethicknessauto",action='store_true',default=False,help="""
		If --subtiltsdir is provided (and if --icethickness is *not* provided), the thickness of the 
		specimen in Z will be calculated by computing the difference between the largest 
		and the smallest Z coordinate found in the header of the subtiltseries.""")
	
	parser.add_argument("--framexsize",type=int,default=0,help="""This correspond to the X
		size in pixes of the images/frames in the raw tilt series; that is, the size of the entire frame
		along the X axis (perpendicular to the direction of the tilt axis in the aligned tilt series).
		It is used to calculate the distance of each particle (subtiltseries) to the tilt axis.
		This MUST be provided if --subtiltsdir is provided.
		Othwerwise, it will be read from the header of the images provided.""")
	
	parser.add_argument("--phaseflipwhole", action='store_true',default=False,help="""This 
		will perform phase flipping on the entire image for each image in an aligned tilt 
		series using the CTF parameters supplied.""")
	
	parser.add_argument("--phaseflipstrips",action='store_true',default=False,help="""This will
		perform phase flipping on images of an aligned tilt series on a strip-by-strip basis,
		assuming the supplied ctf parameters correspond to the proper values at the tilt axis,
		either the same values for all images (--defocus,--ampcont,--cs,--apix,--voltage,--bfactor)
		or a different set for each (--ctfparamsfile), taking into account the tilt angle for 
		each image (--tltfile), which should be supplied through an IMOD-like .tlt file.""")
		
	parser.add_argument("--tltfile",default='',type=str,help="""File containing a list of 
		tilt angles corresponding to the tilt angles of images 0 to n of an aligned
		tilt series""")
	
	parser.add_argument("--ctfparamsfile",type=str,default='',help="""This should be a text file
		with ctf parameters in the following format;
		defocus=value voltage=value cs=value apix=value bfactor=value ampcont=value
		A single space should separate each parameter from the next.
		Do not write any unit symbols for the values; just the numerical value.
		Defocus should be in microns, voltage in kV, apix in angstroms per pixel, and ampcont (amplitude contrast)
		should be a decimal; for example, 0.1 for 10 percent amplitude contrast.
		IF you want to use DIFFERENT PARAMETERS PER IMAGE, then the file must contain
		multiple rows with the different values.
		The first row will be used to phase flip the first image,
		the second row to phase flip the second, etc.""")
	
	parser.add_argument("--defocilist",type=str,default='',help='''Text file containing
		a single column of defocus values in microns. The file should have as many
		defocus values as images in the tiltseries or subtiltseries supplied.''')
	
	parser.add_argument("--defocus", type=float,default=0.0,help="""Default=0. 
		Target defocus at the tilt axis. In the absence of ctfparamsfile(s)
		this value will be assumed to be the defocus at the tilt axis for all tilt images.""")
	
	parser.add_argument("--voltage", type=int,default=200,help="""Default=200. Voltage of
		the microscope with which the images where collected. Supply it to replace the value
		in ctfparamsfile(s), or if ctfparamsfile(s) are lacking altogether.""")
	
	parser.add_argument("--cs", type=float,default=2.1,help="""Default=2.1. Cs of the microscope
		with which the images were collected. Supply it to replace the value in ctfparamsfile(s), 
		or if ctfparamsfile(s) are lacking altogether.""")
	
	parser.add_argument("--apix",type=float,default=0.0,help="""Default=whatever is on the header
		of the images. Sampling of the images in angstroms/pixel. 
		Supply --apix here to replace the value in ctfparamsfile(s), or if ctfparamsfile(s) 
		are lacking altogether.""")	
	
	parser.add_argument("--bfactor",type=int,default=1000,help="""Default=1000. Bfactor or
		temperature factor to use. Supply it to replace the value
		in ctfparamsfile(s), or if ctfparamsfile(s) are lacking altogether.""")
	
	parser.add_argument("--ampcont",type=float,default=0.05,help="""Default=0.05. Amplitude 
		contrast to use for CTF correction phase flipping. Supply it to replace the value
		in ctfparamsfile(s), or if ctfparamsfile(s) are lacking altogether.""")
	
	parser.add_argument("--nozcorrection",action='store_true',default=False,help="""If you 
		turn on this option and --subtiltsdir is provided, the position in Z of each subtomogram
		will not be considered for CTF correction""")
		
	(options, args) = parser.parse_args()
	
	#print "options are", options
	
	if options.tiltseries and options.subtiltsdir:
		print """ERROR: You either 1) supply a tiltseries with .mrc, .st, or .ali extension (in MRC format),
		or with .hdf extension (in HDF format), for an entire tomogram, 2) a stem (a 'string') common 
		to all individual .mrc or .hdf images corresponding to a tiltseries, OR 3) a directory 
		with subtiltseries in .hdf format for individual subtomograms. You cannot supply both
		--tiltseries, --subtiltsdir and --imagestem at the same time. Pick one."""
		sys.exit()
		
	if options.tiltseries and options.imagestem:
		print """ERROR: You either 1) supply a tiltseries with .mrc, .st, or .ali extension (in MRC format),
		or with .hdf extension (in HDF format), for an entire tomogram, 2) a stem (a 'string') common 
		to all individual .mrc or .hdf images corresponding to a tiltseries, OR 3) a directory 
		with subtiltseries in .hdf format for individual subtomograms. You cannot supply both
		--tiltseries, --subtiltsdir and --imagestem at the same time. Pick one."""
		sys.exit()
	
	if options.imagestem and options.subtiltsdir:
		print """ERROR: You either 1) supply a tiltseries with .mrc, .st, or .ali extension (in MRC format),
		or with .hdf extension (in HDF format), for an entire tomogram, 2) a stem (a 'string') common 
		to all individual .mrc or .hdf images corresponding to a tiltseries, OR 3) a directory 
		with subtiltseries in .hdf format for individual subtomograms. You cannot supply both
		--tiltseries, --subtiltsdir and --imagestem at the same time. Pick one."""
		sys.exit()
	
	if options.reconstructor and options.reconstructor != 'None' and options.reconstructor != 'none': 
		options.reconstructor=parsemodopt(options.reconstructor)
		
	'''
	Figure out apix and nimgs to process in the tiltseries (or in each subtiltseries).
	'''
	apix = 0.0
	nimgs = 0
	subtilts = []
	autoIcethickness=0
	nx=0
	if options.subtiltsdir:
		if not options.framexsize:
			print "ERROR: provide a value for options.framexsize if processing --subtiltsdir"
			sys.exit(1)
			
		findir = os.listdir( options.subtiltsdir )
		
		nimgs = 0
		zs = []
		for f in findir:
			if '.hdf' in f:
				stsfile = options.subtiltsdir + '/' + f
				subtilts.append( stsfile )
				print "\nFound subtiltseries", f
				
				
				stshdr = EMData( stsfile, 0, True )
				hdrcoords = stshdr['ptcl_source_coord']
				nx = stshdr['nx']
				print "hdrcoords are", hdrcoords
				z = hdrcoords[-1]
				zs.append( z )
		
		autoIcethickness = max( zs ) -  min( zs )
		
		icefile=options.path+'/autoicethickness.txt'
		os.system( 'touch ' + icefile )
		f=open(icefile,'w')
		line=[str(autoIcethickness)+'\n']
		f.writelines(line)
		f.close()
		
		nimgs = EMUtil.get_image_count( subtilts[0] )
		apix = EMData( subtilts[0], 0, True)['apix_x']
				
	elif options.imagestem:
		if not options.icethickness:
			print "ERROR: provide a value for options.icethickness if processing --imagestem or --tiltseries."
			sys.exit(1)
		
		findir = os.listdir( os.getcwd )
		
		imgs = []
		for f in findir:
			if '.hdf' in f or '.mrc' in f:
				imgs.append( f )
		
		nimgs = len( imgs )	
		apix = EMData( imgs[0], 0, True )['apix_x']
		framexsize = EMData( imgs[0], 0, True )['nx']
				
	elif options.tiltseries:
		if not options.icethickness:
			print "ERROR: provide a value for options.icethickness if processing --imagestem or --tiltseries."
			sys.exit(1)
		
		if '.mrc' in options.tiltseries or '.st' in options.tiltseries or '.ali' in options.tiltseries:
			nimgs =  EMData( options.tiltseries, 0, True )['nz']
		elif '.hdf' in options.tiltseries:
			nimgs = EMUtil.get_image_count( options.tiltseries )
			
		apix = EMData( options.tiltseries, 0, True )['apix_x']
		framexsize = EMData( imgs[0], 0, True )['nx']
		
	if options.apix:
		apix = options.apix
		
	if options.framexsize:
		framexsize = options.framexsize
	
	'''
	If no crashes till now, make the directory where to create the database where the results will be stored
	'''
	from e2spt_classaverage import sptmakepath
	options = sptmakepath (options, 'sptctf')
	
	icethickness=0
	if options.icethickness: 
		if int(options.icethickness) < int(nx):
			icethickness = int(nx)
		else:
			icethickness = options.icethickness
	
	elif options.icethicknessauto:
		icethickness = autoIcethickness
	
	angles = []
	if options.tltfile:
		angles = getangles( options )

	nangles = len( angles )
	
	if nangles != nimgs:
		print "ERROR: The number of angles %d does not coincide with number of images %d" % ( nangles, nimgs )
		sys.exit(1)
		
	'''
	Read/generate the CTF parameters to use
	'''
	ctfs = genctfparamlines( options, apix, nimgs, angles )
	

	'''
	Log current run of the program
	'''
	logger = E2init(sys.argv, options.ppid)
	
	
	'''
	If input consists of subtiltseries, adjust the defocus value in 'ctfs' based on the particle's
	X and Z coordinates in the tomogram
	'''
	if options.subtiltsdir:
		correctsubtilt( options, subtilts, angles, ctfs, apix, nangles, nimgs, framexsize, icethickness )
	
	'''
	#CTF correction for an entire tilt series is different than for a subtiltseries of a subtomogram,
	at least when using the same correction for the entire image. In the latter case (subtiltseries)
	the images are only accepted as a stack, and you want to account for differences in z-height.
	'''
	if not options.subtiltsdir:	
		
		imagefilenames = {}

		'''
		#If input consists of individual image files, find them and put them into an imagefilenames dictionary
		'''
		if options.imagestem:
			print """\nI will process all images in the current directory containing the following string""", options.imagestem
			findir=os.listdir(os.getcwd())
		
			for f in findir:
				if options.imagestem in f:
					imagestem = f.replace('.mrc','')
					imagefilenames.update({imagestem:f})
			print "\nFound these many tilt images",len(imagefilenames)
		
		'''
		#If input is a tiltseries, unstack it, then put the individual images into an imagefilenames dictionary
		'''
		if options.tiltseries:
			if '.st' in options.tiltseries or '.mrc' in options.tiltseries or '.ali' in options.tiltseries:
				print "\nI will process this tiltseries", options.tiltseries
				#if not options.dontunstack:
				cmd = 'cd ' + options.path + ' && e2spt_tilstacker.py --unstack --options.input=' + options.tiltseries
	
				print "\nCmd to extract tilts is", cmd		
				p = subprocess.Popen( cmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				text = p.communicate()	
				p.stdout.close()
			
				nimgs = EMData(options.tiltseries,0,True)['nz']
				print "\n(e2spt_ctf.py)(main) There are these many tilts in the tilt series", n
		
				print "\nCurrent dir is", os.getcwd()
				findir = os.listdir( os.getcwd() + '/' + options.path + '/sptstacker/' )
				#imagefilenames = []
				for f in findir:
					imagestem = f.replace('.mrc','')
					imagefilenames.update({imagestem:f})
		
				nfiles = len(imagefilenames)
				if nimgs != nfiles:
					print """\n(e2spt_ctf.py)(main) WARNING: It seems like not all the images
					in the tilt series were properly unstacked. nimages and nfiles are""",nimages,nfiles
	
			else:
				print "\nTilt series has to be in .st or .mrc or .ali file format"
				sys.exit()	
		
		
		'''
		#Verify that you have CTF parameters for each image, returned from genctfparamlines
		'''
		if len( paramlines ) != len( imagefilenames ):
			print """(e2spt_ctf.py)(main) ERROR: It seems like you have less parameter
				lines in the --ctfparamsfile than images in the tilt series.
				You need one line of parameters per image.
				To apply the same correction to all images, enter the parameters directly,
				through --defocus, --apix, --ampcont, --bfactor, --cs and --voltage"""
			sys.exit()
		
		
		pp=0
		'''
		#Once the images are organized, proceed to CTF correct them. Loop over all images.
		'''
		for imagestem in imagefilenames.keys():
			print "\nWorking on image",pp,imagefilenames[imagestem]
			ctf=None
			
			'''
			#The first scheme is to apply the same correction to each entire image. This will only work if
			#the images are "low tilt" (or if they're all taken at high mag), such that the defocus 
			#gradient results irrelevant.
			'''
			if options.phaseflipwhole:
				print "\nI will phaseflip these images ignoring the defocus gradient",len(imagefilenames)
			
				if ctf:
					print "\nLoading image", imagefilenames[imagestem]
					img = EMData( imagefilenames[imagestem] )
					
					originalx = img['nx']
					originaly = img['ny']
					
					recrop=0
					if img['nx'] != img['ny']:
						print "\n\nThe image is not square and therefore will be padded so that all 4 sides are equal to its largest dimension"
						size=max(img['nx'],img['ny'])
						print "The image will be padded into a square of size", size
						img=padder(options,img,size,size)
						img.write_image(options.path + '/' + imagefilenames[imagestem].replace('.mrc','_sq.mrc') )
						recrop=1
					
					print "The returned ctf to use is", ctf
					print "Of type", type(ctf)
			
					ret = phaseflipper(options,img,ctf)
					flippedtilt = ret[0]
					
					if recrop:
						print "The image will be recropped into its oiriginal size of", originalx,originaly
						flippedtilt=padder(options,img,originalx,originaly)
					
					print "\nPhase-flipped image"
		
					outflipimg = options.path + '/' + imagefilenames[imagestem].replace('.mrc','_flip.mrc')
			
					flippedtilt['spt_phaseflipped']='yes'
			
					flippedtilt.write_image(outflipimg,0)
					print "\nWrote flipped image to", outflipimg
			
					#cmdmrc = 'e2proc2d.py ' + outflipimg + ' ' + outflipimg.replace('.hdf','.mrc') + ' --mrc16bit'
			
					#p = subprocess.Popen( cmdmrc , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
					#text = p.communicate()	
					#p.stdout.close()
					#print "\nCoverted to mrc"
				else:
					print "ERROR: CTF for this image file is blank!",imagefilenames[imagestem]
					sys.exit()
			elif options.phaseflipstrips:
				pass
		

				
	
		'''
		Compile a CTF corrected tilt series if the input was a tilt series (opposed to individual images)
		'''	
		if options.tiltseries:
			outflipseries = options.path + '/tiltseriesflipped.mrc'
			if '.mrc' in options.tiltseries:
				outflipseries = options.tiltseries.replace('.mrc','_flip.mrc')

			if '.st' in options.tiltseries:
				outflipseries = options.tiltseries.replace('.st','_flip.st')

			if '.ali' in options.tiltseries:
				outflipseries = options.tiltseries.replace('.ali','_flip.ali')
			
			outflipseries = options.path + '/' + outflipseries
	
		cmdst = 'newstack ' + options.path + '/*_flip.mrc ' + outflipseries
		
		print "\nCreating flipped tilt series with this command", cmdst
		p = subprocess.Popen( cmdst , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		text = p.communicate()	
		p.stdout.close()
		
	"""
	'''
	Correct individual particle subtiltseries.
	'''
	if options.subtiltsdir:
		findir = os.listdir( options.subtiltsdir )
		
		subtomofiles = []
		for f in findir:
			if '.hdf' in f:
				subtomofiles.append( f )
		
		ntilts = EMUtil.get_image_count( subtomofiles[0] )
		
		'''
		Verify that you have parameters for each tilt in each subtomogram
		'''		
		if len( paramlines ) != ntilts:
			print '''ERROR: (e2spt_ctf.py)(main) ERROR: It seems like you have less parameter
				lines in the --ctfparamsfile than images in each subtilt series.
				You need one line of parameters per image in the subtilt series.
				To apply the same correction to all images, enter the parameters directly,
				through --defocus, --apix, --ampcont, --bfactor, --cs and --voltage.'''
			sys.exit()
		
				
		'''
		Proceed to CTF correction. Loop over all subtomograms.
		'''
		for stomof in subtomofiles:
			ntilts = EMUtil.get_image_count( stomof )
			
			'''
			Loop over all tilts
			'''
			for i in range( ntilts ):
		
				ctf = None
				pp = 0
				
				if options.phaseflipwhole:

					if options.ctfparamsfile:
						
						if 1:
							pass
					
						else:
							pline = paramlines[pp]
							ctf = ctfparamparser(pline)
							pp+=1
			
					elif options.infodir:
						infodir = options.indir
						
						infofiles = []
						findirinfo = os.listdir( infodir )
						
						for inf in findirinfo:
							if '_info.json' in inf:
								infofile = infodir + '/' + inf
								infofiles.append( infofile )
						
						if len(infofiles) != ntilts:
							print '''\nERROR: You must have as many info files as tilt images in each 
								subtilt series.'''
							print "\nNumber of infofiles is", len(infofiles)
							print "\nNumber of tilts in current subtiltseries is", ntilts
							sys.exit()
							
						else:	
							for infof in infofiles:
								infofstem = infof.split('/')[-1].split('_info')[0]
								if infostem == stem:
									#infofile = infodir + '/' + infof
									js = js_open_dict( infofile )
									ctf = js["ctf_frame"][1]
			
					elif pp==0:
						ctf=ctfparamparser(paramlines[0])
				
					img = EMData(f,ns)
					ret = phaseflipper(options,img,ctf)
					phfimg=ret[0]
				
				else:
					pass
	
		"""
		
	E2end(logger)
	
	return


def correctsubtilt( options, subtilts, angles, ctfs, apix, nangles, nimgs, framexsize, icethickness ):
	
	ii = 0
	globalAvgDefErrors=[]
	for sts in subtilts:
		imghdr = EMData( sts, 0, True )
		
		#print " icethickness is", icethickness
		#print "ice < int(imghdr['nx'])", int(icethickness) < int(imghdr['nx'])
		#print "ice and type are", type(options.icethickness), icethickness
		#print "imghdrnx and type are",type(imghdr['nx']),imghdr['nx']
		
		if icethickness and int(icethickness) < int(imghdr['nx']):
			print """\nThe ice must be thick enough to contain a layer of molecules and will be set
			to half the X and  Y sides of the images. It must be => than %d pixels.""" %( imghdr['nx'] )
			sys.exit()
			
		coords = imghdr['ptcl_source_coord']
		coordx = coords[0]
		coordz =  coords[-1]
		nx = imghdr['nx']
		
		if options.verbose:
			print "Fixing subtomogram", ii
		
		n = EMUtil.get_image_count( sts )
		if n != nangles:
			print "ERROR: The number of angles %d does not coincide with number of images %d" % ( nangles, nimgs )
			sys.exit(1)
		
		flippedsts = options.path + '/' + os.path.basename( sts ).replace('.hdf','_PHFLIP.hdf')
		phfimgs = []
		defocuserrors=[]
		for m in range( n ):
			img = EMData( sts, m )
			angle = angles[ m ]
			ctf = ctfs[ angle ]
			print "\n\nUncorrected defocus is", ctf.defocus
			
			px = ( coordx - framexsize/2.0 ) * apix/10000
			dzx = px * numpy.sin( math.radians( angle ) )
			
			
			print"\nangle is", angle
			
			newdefocus = ctf.defocus + dzx 
			print "First corrected defocus is", newdefocus 
			
			pz = ( coordz - icethickness/2.0 ) * apix/10000	
			dzz = -1 * pz * numpy.cos( math.radians( angle ) )
				
			if not options.nozcorrection and icethickness and icethickness > nx:
					
				newdefocus += dzz
	
				print "Second corrected defocus is", newdefocus
			else:
				print "\n!!!!!!!!!\ndid NOT correct dzz", dzz
				
			
			finalctf = EMAN2Ctf()
			finalctf.from_dict({ 'defocus':newdefocus, 'bfactor':ctf.bfactor, 'ampcont':ctf.ampcont, 'apix':ctf.apix, 'voltage':ctf.voltage, 'cs':ctf.cs })	
			
			try:
				actualctf = img['ctf']
				print "\nactual defocus is", actualctf.defocus
				defocuserror = actualctf.defocus - newdefocus
				print "\nTherefore, defocus error is", defocuserror
				defocuserrors.append( math.fabs(defocuserror) )			#The average error needs to sum all positive errors
				
			except:
				pass
			
			ret = phaseflipper( options,img,finalctf )
			imgflipped = ret[0]
			imgflipped['ctf'] = finalctf
			
			#print "Flipped outstack to write is", flippedsts
			#print "imgflipped and type are", imgflipped, type(imgflipped)
			#print "index to write is", m
			
			imgflipped.write_image( flippedsts, m )	
			phfimgs.append( imgflipped )
		
		sts3d = options.path + '/' + os.path.basename( sts ).replace('.hdf','_PHFLIP3D.hdf')
		
		rec = reconstruct3d( options, phfimgs, apix )
		
		if defocuserrors:
			defocuserrorsAvg=sum(defocuserrors)/len(defocuserrors)
			rec['spt_avgDefocusError']=defocuserrorsAvg
			
			globalAvgDefErrors.append(defocuserrorsAvg)
		else:
			print "Defocus errors is empty!", defocuserrors
			sys.exit()
		
		if options.save3d:
			rec.write_image( sts3d , 0 )
			
		ii+=1
	
	lines=[]	
	if globalAvgDefErrors:
		globalAvgDefError=sum(globalAvgDefErrors)/len(globalAvgDefErrors)
		lines.append('Global average error = '+str(globalAvgDefError)+'\n')
		for error in globalAvgDefErrors:
			line = str(error)+'\n'
			lines.append(line)
		
		defErrorsFile=options.path+'/defocusErrorAvg.txt'
		g=open(defErrorsFile,'w')
		#line=[str(globalAvgDefError)+'\n']
		g.writelines(lines)
		g.close()
	
	return


def reconstruct3d( options, phfimgs, apix ):
	
	box = phfimgs[0]['nx']
	
	originalboxsize = box
	
	if options.pad3d:
		if options.pad2d:
			if options.pad3d > options.pad2d:
				box = box*options.pad3d
			else:
				box = box*options.pad2d
		else:
			box = box*options.pad3d			
	elif options.pad2d:
		box = box*options.pad2d
	
	mode='gauss_2'
	if options.reconstructor:
		if len(options.reconstructor) > 1:
			if 'mode' in options.reconstructor[-1]:
				mode=options.reconstructor[-1]['mode']
				
				print "\nThe reconstructor mode has been changed from default to", mode
				#sys.exit()
						
	r = Reconstructors.get(options.reconstructor[0],{'size':(box,box,box),'sym':'c1','verbose':True,'mode':mode})
	r.setup()

	k=0
	weight = 1.0
	for p in phfimgs:
		#print "Adding projection k", k
		#print "Whose min and max are", p['minimum'], p['maximum']
		#print "The size of the prj to insert is", p['nx']
		if options.pad2d:
			p = clip2D( p, box )
		
		pm = r.preprocess_slice(p,p['xform.projection'])
		r.insert_slice(pm,pm['xform.projection'],weight)
		k+=1
	
	rec = r.finish(True)

	rec['apix_x']=apix
	rec['apix_y']=apix
	rec['apix_z']=apix
	rec['origin_x']=0
	rec['origin_y']=0
	rec['origin_z']=0

	'''
	Preserve SPT parameters in header
	'''
	names = phfimgs[0].get_attr_dict()
	for name in names:
		if 'spt_' in name or 'tomogram' in name or 'ptcl_source_coord' in name or 'spt' in name:
			rec[ name ] = names[ name ]

	return rec


def clip2D( img, size ):
	
	imgxc = img['nx']/2
	imgyc = img['ny']/2
	
	Rimg =  Region( (2*imgxc - size)/2, (2*imgyc - size)/2, 0, size , size , 1)
	img.clip_inplace( Rimg )
	
	return img
	

def genctfparamlines( options, apix, nimgs, angles ):

	ctfs={}
	
	'''
	#Determine where to get ctf parameters from
	'''
	if options.ctfparamsfile:
		print "\nI'll parse the ctfparamsfile", options.ctfparamsfile
		g = open(options.ctfparamsfile,'r')
		initiallines = g.readlines()
		g.close()
	
		print "\nCTF will be derived from --ctfparamsfile",options.ctfparamsfile
		
		if len( initiallines ) != nimgs:
			print """ERROR: The number of lines in the file provided through 
				--ctfparamsfile should match the number of images in the tiltseries
				or in each subtiltseries provided.
				To use the same parameters for all images provide them
				explicitly through --cs, --bfactor,--voltage, --defocus, --ampcont and, optionally --apix
				(this apix will be read from the header if not provided, so make sure it's correct). """
			sys.exit(1)
		
		kk=0
		for line in initiallines:
			if len(line) > 25 and len(line.split(' ')) == 6:
				ctf = ctfparamparser( line )
				angle = angles[ kk ]
				ctfs.update( { angle:ctf } )
				kk+=1
		
	elif options.infodir:
		print "\nCTF will be read from info files in --infodir",options.infodir
		
		findirinfo = os.listdir( options.infodir )
		
		kk=0
		ctflines = []
		for inf in findirinfo:
			if '_info.json' in inf:
				infofile = options.infodir + '/' + inf
				#infofstem = inf.split('_info')[0]
				#infofiles.update( { infofstem:infofile } )
				
				js = js_open_dict( infofile )
				ctf = js["ctf_frame"][1]
				js.close()
				
				line = 'defocus=' + str( ctf['defocus'] ) + 'ampcont=' + str( ctf['ampcont'] ) + 'voltage=' + str( ctf['voltage'] ) + 'cs=' + str( ctf['cs'] ) + 'apix=' + apix + 'bfactor=' + str( ctf['bfactor'] )
				ctflines.append( line )
				angle = angles[ kk ]
				ctfs.update( { angle:ctf } )
				kk+=1

		if len( ctflines ) != nimgs:
			print """ERROR: The number of _info.json files inside the directory provided
				through --infodir should match the number of images in the tiltseries
				or in each subtiltseries provided.
				To use the same parameters for all images provide them
				explicitly through --cs, --bfactor,--voltage, --defocus, --ampcont and, optionally --apix
				(this apix will be read from the header if not provided, so make sure it's correct). """
			sys.exit(1)
	
	elif options.defocilist:
		print "\nI'll parse the defocilist %s and fill the other CTF parameters with default values" % ( options.defocilist )
		g = open(options.defocilist,'r')
		defoci = g.readlines()
		g.close()
		
		if len( defoci ) != nimgs:
			print """ERROR: The number lines in the file provided
				through --defocilist should match the number of images in the tiltseries
				or in each subtiltseries provided.
				To use the same parameters for all images provide them
				explicitly through --cs, --bfactor,--voltage, --defocus, --ampcont and, optionally --apix
				(this apix will be read from the header if not provided, so make sure it's correct). 
				"""
			sys.exit(1)
		
		if options.voltage and options.cs and apix and options.bfactor and options.ampcont:
			print """\nExplicit parameters --cs,--apix,--bfactor,--voltage and --ampcont 
				will be used."""
				
			kk=0
			for d in defoci:
				de = d.replace('\n','').replace(' ','').replace('\t','')
				line = 'defocus=' + str( de ) + 'ampcont=' + str( options.ampcont ) + 'voltage=' + str( options.voltage ) + 'cs=' + str( options.cs ) + 'apix=' + str( apix ) + 'bfactor=' + str( options.bfactor )
				ctf = ctfparamparser( line ) 
				angle = angles[ kk ]
				ctfs.update( { angle:ctf } )
				kk+=1
		else:
			print """\nERROR: There's nothing to do. If you supply --defocilist, you also 
			have to provide the following 4 parameters:
			--ampcont,--cs,--voltage,--bfactor (and, optionally, --apix, if the images 
			don't have the correct apix in their headers.)"""
			sys.exit(1)
		
	else:
		if options.voltage and options.cs and options.defocus and apix and options.bfactor and options.ampcont:
			print """\nExplicit parameters --defocus,--cs,--apix,--bfactor,--voltage and 
				--ampcont will be used."""
			
			print '\ndefocus to set is', options.defocus	
			line = 'defocus=' + str( options.defocus ) + 'ampcont=' + str( options.ampcont ) + 'voltage=' + str( options.voltage ) + 'cs=' + str( options.cs ) + 'apix=' + str( apix ) + 'bfactor=' + str( options.bfactor )
			ctf = ctfparamparser( line )
			print '\nafter parsing and making ctf object, set defocus is', ctf.defocus
			for kk in range( nimgs ):
				angle = angles[ kk ]
				ctfs.update({ angle : ctf } )
		
		else:
			print """\nERROR: There's nothing to do. If you don't provide --infodir, 
			--ctfparamsfile or --defocilist, you have to provide the following 5 parameters:
			--defocus, --ampcont,--cs,--voltage,--bfactor (and, optionally, --apix, if the images 
			don't have the correct apix in their headers.)"""
			sys.exit(1)
	
	return ctfs


def getangles( options ):
	
	angles = []
	
	f = open( options.tltfile, 'r' )
	lines = f.readlines()
	f.close()
	
	for line in lines:
		line = line.replace('\t','').replace('\n','')
	
		if line:
			angles.append( float(line) )
	
	if options.verbose > 9:
		print "\n(e2spt_ctf.py)(getangles) angles are", angles

	return angles


def padder(options,img,sizex,sizey):
	
	xc = img['nx']/2
	yc = img['ny']/2
	
	print "\nThe center of the image to pad is at", xc,yc
	
	if int(sizex) % 2 != 0:
		sizex += 1
		
	if int(sizey) % 2 != 0:
		sizey += 1
	
	print "\nThe size to pad to is", sizex,sizey		
	#if side % 8 and parameters['box_mult_of_8']=='yes':
	#	factor = int(side/8)*8 + 8
	#	print "\nKKKKKK\nKKKKKKK\nKKKKKKK The box has been changed to the closest (and larger) power of 2, see", side
			
	r = Region((2*xc - sizex)/2,(2*yc - sizey)/2, sizex, sizey)
	
	print "Therefore,region is",r
	imgp = img.get_clip(r)
	
	return imgp
	
	
def ctfparamparser( pline ):
	
	defocus = pline.replace('\n',' ').split("defocus=")[-1].split(' ')[-1]
	ampcont = pline.replace('\n',' ').split("ampcont=")[-1].split(' ')[-1]
	voltage = pline.replace('\n',' ').split("voltage=")[-1].split(' ')[-1]
	cs = pline.replace('\n',' ').split("cs=")[-1].split(' ')[-1]
	apix = pline.replace('\n',' ').split("apix=")[-1].split(' ')[-1]
	bfactor = pline.replace('\n',' ').split("bfactor=")[-1].split(' ')[-1]

	params = {'ampcont':ampcont,'apix':apix,'bfactor':bfactor,'cs':cs,'defocus':defocus,'voltage':voltage}
	print "\n(e2spt_ctf.py)(ctfparamparser) The parsed parameters are"
	for key in params.keys():
		print key + '=' + params[key] 
	
	ctf = EMAN2Ctf()
	#ctf.from_dict({'defocus':params['defocus'],'bfactor':params['bfactor'],'ampcont':params['ampcont'],'apix':params['apix'],'voltage':params['voltage'],'cs':params['cs']})	
	ctf.from_dict(params)
	
	return ctf

	
def phaseflipper(options,img,ctf):	

	#prj=EMData(imgfile,0)
	prj=img.copy() 
	
	#maskradius = min(prj['nx'],prj['ny'])/2 - 200
	
	#mask=EMData(prj['nx'],prj['ny'])
	#mask.to_one()
	
	#mask.process_inplace('mask.sharp',{'outer_radius':maskradius})
	
	#prj.process_inplace('normalize.mask',{'mask':mask})
	
	#prj.process_inplace('normalize.edgemean')
	
	
	#apix = prj['apix_x']
	#if options.apix:
	#	apix = options.apix
	
	prj_fft = prj.do_fft()

	#ctf = EMAN2Ctf()
	#ctf.from_dict({'defocus':params['defocus'],'bfactor':params['bfactor'],'ampcont':params['ampcont'],'apix':params['apix'],'voltage':params['voltage'],'cs':params['cs']})	
	#ctf.from_dict(params)

	flipim = prj_fft.copy()	
	ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)

	prj_fft.mult(flipim)

	intermediate = prj_fft.copy()

	prj_flipped=prj_fft.do_ift()
	
	return prj_flipped, intermediate
	
	
def tiler(options):	
	n = EMUtil.get_image_count(stack)
	print "The number of slices is", n
	testimg = EMData(stack,0)
	nx = testimg['nx']
	print "The dimenions are", nx
	
	actualsize = nx / options.gridsize
	strips = []
	for i in range(n):
		img = EMData(stack,i)
		nx = img['nx']
		ny = img['ny']
		y1=0
		y2=ny
		while j <= nx:
			Region=(x1,y1,x2,y2)
			img.get_clip(R)
			
		
	return
	
	

	
		
			 

if '__main__' == __name__:
	main()