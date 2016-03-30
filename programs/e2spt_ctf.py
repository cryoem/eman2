#!/usr/bin/python2.7

#
# Author: Jesus Galaz, 11/01/2012; last update 17/mar/2016
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
import collections
import sys
import numpy		 
import math
import e2ctf
	 
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
	
	

	parser.add_argument("--tiltseries", type=str, default='', help="""Aligned tilt series. File format must be MRC and must have .mrc or .st or .ali extension.""")
	
	parser.add_argument("--exclude",type=str,default='',help="""Comma-separated list of image indexes in the --tiltseries to exclude from CTF fitting. For example, --exclude 0,3,4,6,7.""")
	
	parser.add_argument("--skipstripping", type=str, default='', help="""Default=None. Comma-separated list of image indexes to exclude from strip-based fitting (in this case, only global defocus tiling the entire image wil be measured).""")
	
	parser.add_argument("--imagestem",type=str,default='',help="""Default=None. If the images to apply ctf correction on are already unstacked and are individual mrc files, supply a common string to all of them.""")
	
	parser.add_argument("--invert",action='store_true',default=False,help='''Invert the contrast of the output data, compared to the input data.''')
	
	parser.add_argument("--excludeedges",action='store_true',default=False,help='''Ignore 'excedent' (smaller than the width of a strip) at the edge of micrographs after dividing them into strips.''')
	
	parser.add_argument("--mintiles", type=int, default=0, help="""Minimum number of 'good tiles' in strip to consider it.""")
		
	parser.add_argument("--defocusvariationlimit",type=float,default=0.0,help="""default=0.0. Automatically calculated based on --apix and --voltage, but can be overwritten by the user. Total variation in defocus ('depth of focus' in micrometers) tolerated within a strip to still consider it a region of 'constant defocus'.""")
	
	parser.add_argument("--infodir",type=str,default='',help="""Folder typically produced by e2evalimage.py or previous runs of this program containing info.json files, one per tilt image in a tilt series. Each .json file should contain the fitted ctf and all associated parameters for each tilt image.""")
	
	parser.add_argument("--output", type=str, default='',help="""Filename for the output CTF-corrected tilt series.""")
		
	parser.add_argument("--subtiltsdir",type=str,default='',help="""Provide a directory containing individual stacks, where each stack is a 'mini tilt series' or a 'subtilt series' for single particles. Then, each image for each particle in the dir will be phase-phlipped using the ctf parameters you provide. If each image in the subtilt series is at a different defocus, then the parameters should be provided through --ctfparamsfile, whith a different defocus value per row. (There should be as many rows as images in each subtiltseries).""")
		
	parser.add_argument("--path",type=str,default='sptctf',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'sptctf'; for example, sptctf_02 will be the directory by default if 'sptctf_01' already exists.""")
	
	parser.add_argument("--reconstructor", type=str,default="fourier:mode=gauss_2",help="""Default=fourier:mode=gauss_2. The reconstructor to use to reconstruct the tilt series into a tomogram. Type 'e2help.py reconstructors' at the command line to see all options and parameters available. To specify the interpolation scheme for the fourier reconstruction, specify 'mode'. Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 'gauss_5_slow', 'gypergeom_5', 'experimental'. For example --reconstructor=fourier:mode=gauss_5 """)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	parser.add_argument("--maxangerror", type=float,default=180.0,help="""Default=180.0. If the calculated angle exceeds this value, the global defocus will be used.""")
	
	parser.add_argument("--pad2d", type=float,default=0.0,help="""Padding factor to zero-pad
		the 2d images in the tilt series prior to reconstruction.
		(The final reconstructed subvolumes will be cropped to the original size).""")

	parser.add_argument("--pad3d", type=float,default=0.0,help="""Padding factor to zero-pad
		the reconstruction volume. (The final reconstructed subvolumes will be cropped to 
		the original size).""")	
		
	parser.add_argument("--save3d",action='store_true',default=False,help="""If on, the CTF
		corrected subtiltseries will be reconstrcuted into subvolumes and save into a stack.
		Options --reconstructor, --pad2d, --pad3d are used if --save3d is on.""")
		
	parser.add_argument("--save2d",action='store_true',default=False,help="""If on, the CTF
		corrected subtiltseries will be saved as 2-D imag stacks [one per particle].""")
	
	parser.add_argument("--outputstem",type=str,default='',help="""Stem common to all
		output image stacks. For example, if --outputstem=myvirus and --save2d is provided, 
		the phase-flipped images for each subtiltseries wille be saved to myvirus_subtiltptclXXXX.hdf.
		If --save3d is provided, the stack of reconstructed subvolumes will be saved to myvirus_stack3d.hdf""")
	
	parser.add_argument("--savestriptiles",action='store_true',default=False,help="""Saves
		all tiles for all strips, for all images, in one stack per strip.""")
	parser.add_argument("--saveffts",action='store_true',default=False,help="""Saves
		ffts of each average of tiles per strip, for all images.""")
	
	parser.add_argument("--icethickness", type=int,default=0,help="""This corresponds
		to the Z dimension in pixels of the reconstructed raw tomogram (uncropped), at the same binning
		(sampling) as the provided tiltseries, images or subtiltseries.
		This value MUST be provided, only if --subtiltsdir is given.
		""")
	
	parser.add_argument("--autofit", action='store_true', default=False,help="""Runs automated
		CTF fitting on the input images, based on tiling.""")
	
	parser.add_argument("--firstfitglobal", action='store_true', default=False,help="""Default=False.
		Supplying this option will tile the entire image first (for each tilt angle) and find the average
		defocus. Then it will use that value to provide an educated 'guess' during stripe-by-stripe 
		fitting for each image.""")
	
	parser.add_argument("--tilesize",type=int,default=512,help="""Tile size to use for strips
		when --autofit is provided.""")
		
	parser.add_argument("--defocusmin",type=float,default=0.0,help=""" If --autofit, minimum autofit defocus. Default=0.0, not used. A value will be estimated based on tilt angle and distance from the tilt axis.""")
	
	parser.add_argument("--targetdefocus",type=float,default=0.0,help=""" Default=0.0, not used. Target defocus value at the tilt axis in micrometers for the images in the tiltseries.""")

	parser.add_argument("--defocusmax",type=float,default=0.0,help="""Default=0.0, not used. If --autofit, maximum autofit defocus. A value will be estimated based on tilt angle and distance from the tilt axis.""")
		
	parser.add_argument("--stripstep",type=int,default=0,help="""Default=0 (not used; no overlap). This will determine the
		amount of strips and the overlap between them for defocus estimation. For example, for a 4000x4000 pixels image, a tile size of
		400 would yield 20, not 10 strips, by default. If --stripstep=1 were provided, the
		image would be devided into 4000-400=3600 strips. The first strip would go from pixel
		0 to pixel 400, the second strip from pixel 1 to pixel 401, the third from pixel 2
		to 402, etc... up to the last strip going from pixel 3600 to 4000.""")
	
	parser.add_argument("--tileoverlap",action='store_true',default=False,help="""If provided, it will cause tiles to overlap by 50 percent in x and y for power spectrum computation (periodogram averaging).""")
	
	parser.add_argument('--subset', type=int, default=0, help='''Requires --subtiltsdir. Specify how many subtiltseries (or particles) to ctf correct. If you specify 10, the first 10 subtiltseires in --subtiltsdir will be corrected. 0 means "process all" because it makes no sense to process none''')

	
	parser.add_argument("--icethicknessauto",action='store_true',default=False,help="""
		If --subtiltsdir is provided (and if --icethickness is *not* provided), the thickness of the 
		specimen in Z will be calculated by computing the difference between the largest 
		and the smallest Z coordinate found in the header of the subtiltseries, plus the size of the specimen, calculated from --radius.""")
	
	parser.add_argument("--radius",type=int,default=0,help="""Radius of the particle in pixels.""")
	
	parser.add_argument("--framexsize",type=int,default=0,help="""This correspond to the X
		size in pixes of the images/frames in the raw tilt series; that is, the size of the entire frame
		along the X axis (perpendicular to the direction of the tilt axis in the aligned tilt series).
		It is used to calculate the distance of each particle (subtiltseries) to the tilt axis, since
		this will induce different shifts in defocus in 3-D for the actual particles. Particles
		right at the tilt axis don't move "up" or "down" as they are tilted.
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
	
	
	parser.add_argument("--prunetest",type=float,default=0.1,help="""Default=0.1.
		Decimal number that indicates the percentage of --tilesize (in terms of side length) 
		to tolerate of 'bad' values (i.e., empty regions of constant density) at the corners, 
		and still include the tile for CTF fitting. For example, if --tilesize=256, and
		--prunetest=0.1, a box of ~25-26 pixels each corner of every tile will be analyzed
		and if the standard deviation of any of the corners is 0, the tile will be excluded.
		To turn off this option supply --prunetest=-1.0. The program automatically adjusts 
		things so that the minimum size of regions at the corners to check will be 4x4 pixels.""")
		
	parser.add_argument("--tltfile",default='',type=str,help="""File containing a list of 
		tilt angles corresponding to the tilt angles of images 0 to n of an aligned
		tilt series""")
	
	parser.add_argument("--coords",default='',type=str,help="""text file containing x y z (or just z) coordinates for the particles, used to calculate icethickness if --icethicknessauto is specified for ctf fitting. NOT needed if --subtiltsdir is provided for ctf correction.""")
	
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
		
	parser.add_argument("--defocustop",action='store_true',default=False,help="""Assumes the signal for defocus measurement (e.g., carbon film) is at the top layer of the tomogram.""")
	
	parser.add_argument("--defocusbottom",action='store_true',default=False,help="""Assumes the signal for defocus measurement (e.g., carbon film) is at the top layer of the tomogram.""")
	
	(options, args) = parser.parse_args()
	
	#print "options are", options
	
	
	errordetector(options)
	
	
	'''
	Log current run of the program
	'''
	logger = E2init(sys.argv, options.ppid)
	
	
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
	
	'''
	If no crashes till now, make the directory where to create the database where the results will be stored
	'''
	from e2spt_classaverage import sptmakepath
	options = sptmakepath (options, 'sptctf')
	
	
	'''
	Store used parameters in a text file
	'''
	from e2spt_classaverage import writeParameters
	cmdwp = writeParameters(options,'e2spt_ctf.py', 'sptctf')
	
	xs = []
	ys = []
	zs = []
	ptclnx = 0
	
	if options.subtiltsdir:
		procsubtiltsdir( options )
				
	elif options.imagestem:
		#if not options.icethickness:
		#	print "ERROR: provide a value for options.icethickness if processing --imagestem or --tiltseries."
		#	sys.exit(1)
		
		findir = os.listdir( os.getcwd() )
		
		imgs = []
		for f in findir:
			if options.imagestem in f:
				imgs.append( f )
		
		if imgs:
			nimgs = len( imgs )	
			apix = EMData( imgs[0], 0, True )['apix_x']
			framexsize = EMData( imgs[0], 0, True )['nx']
		else:
			print "ERROR: no images found with stem", options.imagestem
			
	elif options.tiltseries:
		#if not options.icethickness:
		#	print "ERROR: provide a value for options.icethickness if processing --imagestem or --tiltseries."
		#	sys.exit(1)
		print "tiltseries is", options.tiltseries
		print "test to see if '.mrcs'", '.mrcs' in options.tiltseries
		
		if '.hdf' in options.tiltseries or '.mrcs' in options.tiltseries:
			print "series is .mrcs or .hdf"
			nimgs = EMUtil.get_image_count( options.tiltseries )
			print "nimgs is therefore", nimgs
		
		elif '.mrc' in options.tiltseries or '.st' in options.tiltseries or '.ali' in options.tiltseries:
			
			nimgs =  EMData( options.tiltseries, 0, True )['nz']
		
		if options.exclude:
			nimgsex = len(options.exclude.split(','))
			nimgs -= nimgsex
			print "but --exclude is %s and therefore %d images will be excluded and the final image count is %d" %( options.exclude, nimgsex, nimgs )
		
		hdr = EMData( options.tiltseries, 0, True )
		apix = hdr['apix_x']
		framexsize = hdr['nx']
		
	if options.apix:
		apix = options.apix
		
	if options.framexsize:
		framexsize = options.framexsize
	
	
	

	
	icethickness = 0
	if options.icethickness:
		if ptclnx:
			if int( options.icethickness ) < int(ptclnx):
				icethickness = int(ptclnx)
		else:
			icethickness = options.icethickness
	
	elif options.icethicknessauto:
		
		if options.coords:
			f=open(options.coords,'r')
			lines=f.readlines()
			zs = [ int(float( line.replace('\t',' ').split(' ')[-1].replace('\n','') )) for line in lines ]
			f.close()
		
		if options.subtiltsdir or options.coords:
			
			radius = ptclnx/2.0
			
			if options.radius:
				radius = options.radius
			elif not options.radius:
				print "\nWARNING: --icethicknessauto requires --radius. Since it wasn't provided, half the box size of the subtiltseries willbe assumed", radius
			
			maxz = max( zs ) 
			minz =  min( zs ) 
			print "maxz", maxz
			print "minz", minz
			print "radius", radius
			autoIcethickness = maxz - minz + 2*radius
			icethickness = autoIcethickness	
			print "autoicethickness therefore is", autoIcethickness
			
			icefile = options.path+'/autoicethickness.txt'
			f = open( icefile, 'w')
			
			line = [str(autoIcethickness)+'\n']
			f.writelines(line)
			f.close()
					
		elif not options.subtiltsdir and not options.coords:
			print "\nWARNING: --icethicknessauto requires --subtiltsdir or --coords"
			
	elif not options.icethicknessauto:
		print "WARNING: No icethickness provided, and --icethicknessauto is also turned off."
			
	angles = []
	if options.tltfile:
		angles = getangles( options )

	nangles = len( angles )

	#indxstoexclude=[int(i) for i in options.exclude.split(',')]

	#if nangles != nimgs-len(indxstoexclude):
	if nangles != nimgs:
		print "\nERROR: The number of angles %d does not coincide with number of images %d" % ( nangles, nimgs)
		sys.exit(1)
	else:
		pass
		#print "\nnumber of images to exclude",len( indxstoexclude )
	
	imagefilenames = {}
	
		
	'''
	#If input consists of individual image files, find them and put them into an imagefilenames dictionary
	'''
	if options.imagestem:
		print """\nprocessing all images in the current directory containing the following string""", options.imagestem
		findir=os.listdir(os.getcwd())
		
		
		kk=0
		for f in findir:
			if options.imagestem in f:
				imagestem = f.replace('.mrc','')
				imagefilenames.update({imagestem:[f,angles[kk]]})
				kk+=1
		print "\nFound these many tilt images",len(imagefilenames)
	
	'''
	#If input is a tiltseries, unstack it, then put the individual images into an imagefilenames dictionary
	'''
	if options.tiltseries:
		if '.st' in options.tiltseries or '.mrc' in options.tiltseries or '.ali' in options.tiltseries or '.hdf' in options.tiltseries:
			print "\nI will process this tiltseries", options.tiltseries
			
			nimgs = EMData(options.tiltseries,0,True)['nz']
			print "\n(e2spt_ctf.py)(main) There are these many tilts in the tiltseries", nimgs
			
	
			
			
			#if not options.dontunstack:
			
			if nimgs > 1:
				#cmd = 'e2spt_tiltstacker.py --unstack=' + options.tiltseries + ' --tltfile=' + options.tltfile
				
				outname = options.path + '/' + options.tiltseries.replace('.mrc','.hdf')
				outname = options.path + '/' + options.tiltseries.replace('.mrcs','.hdf')
				outname = options.path + '/' + options.tiltseries.replace('.st','.hdf')	
				outname = options.path + '/' + options.tiltseries.replace('.ali','.hdf')
	
				outname = outname.replace('.hdf','_UNSTACKED.hdf')
				
				cmdun = 'e2proc2d.py ' + options.tiltseries + ' ' + outname+ ' --unstacking '
				
				if options.exclude:
					excludefile = options.path + '/exclude.lst'
					f = open( excludefile,'w')
					excludelines = options.exclude.split(',')
					lines = [ line + '\n' for line in excludelines ]
					f.writelines( lines )
					f.close()
					
					cmdun +=  ' --exclude ' + excludefile
					
				
				runcmd( options, cmdun )
					
					
				#if options.outmode:
				#	cmdun += ' --outmode=' + options.outmode
	
			
	
				outnamestem = outname.replace('.hdf','')
				
				#print 'outnamestem is', outnamestem
				#cmdun += ' && mv ' + outnamestem + '* ' + options.path
				
				c = os.getcwd()
				findirroot = os.listdir( c )
				for f in findirroot:
					if outnamestem in f:
						os.rename( f, options.path + '/' + f )
				
				
				#print "\nCmd to extract tilts is", cmdun	
				#p = subprocess.Popen( cmdun , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				#text = p.communicate()	
				#p.stdout.close()
				
				#p = subprocess.Popen( cmd , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				#text = p.communicate()	
				#p.stdout.close()
	
				#print "\nCurrent dir is", os.getcwd()
				
				findir = os.listdir( options.path )
				
				#imagefilenames = []
				
				
				kk=0
				for f in findir:
					imagestem = ''
					if '.mrc' in f:
						imagestem = f.replace('.mrc','')
					if '.st' in f:
						imagestem = f.replace('.st','')
					if '.ali' in f:
						imagestem = f.replace('.ali','')
					if '.hdf' in f:
						imagestem = f.replace('.hdf','')
				
					if imagestem:
						imagefilenames.update({imagestem:[ options.path + '/' + f,angles[kk]]})
						kk+=1
	
				nfiles = len(imagefilenames)
				if nimgs != nfiles:
					print """\n(e2spt_ctf.py)(main) WARNING: It seems like not all the images
					in the tilt series were properly unstacked. nimages and nfiles are""",nimgs,nfiles
				
			else:
				f = options.tiltseries
				imagestem = ''
				if '.mrc' in f:
					imagestem = f.replace('.mrc','')
				if '.st' in f:
					imagestem = f.replace('.st','')
				if '.ali' in f:
					imagestem = f.replace('.ali','')
				if '.hdf' in f:
					imagestem = f.replace('.hdf','')
				
				if imagestem:
					imagefilenames.update({imagestem:[f,angles[0]]})

			print "\n\n\nThese are the imagefilenames to process",imagefilenames,len(imagefilenames)
		else:
			print "\nTilt series has to be in .st or .mrc or .ali extension"
			sys.exit()	
		
		
	'''
	Read or generate the CTF parameters to use
	'''
	
	ctfs = genctfparamlines( options, apix, nimgs, angles, imagefilenames, icethickness )
	
	if not options.defocilist:
		defocusesfile = options.path + '/defocuses.txt'
		os.system( 'touch ' + defocusesfile )
		fd = open(defocusesfile,'w')
		linesd = []
		for angle in angles:
			defocus = ctfs[angle].defocus
				
			lined = str(defocus) + '\n'
			linesd.append( lined )
			
		fd.writelines(linesd)
		fd.close()
	
	
	
	
	print "ctfs len is", len(ctfs)
	print "and ctfs are" 
	for ctf in ctfs:
		print ctf, ctfs[ctf]
	
	if options.tiltseries or options.imagestem:
		'''
		#Verify that you have CTF parameters for each image, returned from genctfparamlines
		'''
		
		if ctfs:
			if len( ctfs ) != len( imagefilenames ):
				print """(e2spt_ctf.py)(main) ERROR: It seems like you have fewer parameter
					lines in the --ctfparamsfile than images in the tilt series.
					You need one line of parameters per image.
					To apply the same correction to all images, enter the parameters directly,
					through --defocus, --apix, --ampcont, --bfactor, --cs and --voltage"""
				sys.exit(1)
		else:
			print """ERROR: There is no CTF information for any of the images. If you were
				using --autofit, this means autofitting failed for all images."""
			sys.exit(1)
				
		
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
			
			pp+=1

				
		if options.phaseflipstrips or options.phaseflipwhole:
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
				
				if '.hdf' in options.tiltseries:
					outflipseries = options.tiltseries.replace('.hdf','_flip.hdf')
			
				outflipseries = options.path + '/' + outflipseries
	
			cmdst = 'newstack ' + options.path + '/*_flip.mrc ' + outflipseries
		
			print "\nCreating flipped tilt series with this command", cmdst
			p = subprocess.Popen( cmdst , shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			text = p.communicate()	
			p.stdout.close()
		
	
	
	
	'''
	If input consists of subtiltseries, adjust the defocus value in 'ctfs' based on the particle's
	X and Z coordinates in the tomogram
	'''
	if options.subtiltsdir:
		if not options.save3d:
			if not options.save2d:
				print "ERROR: either --save2d or --save3d must be provided when --subtiltsdir is provided"
			elif options.save2d:
				pass
		elif options.save3d:
			pass
			
		if not options.save2d:
			if not options.save3d:
				print "ERROR: either --save2d or --save3d must be provided when --subtiltsdir is provided"
			elif options.save3d:
				pass
		elif options.save2d:
			pass
		
		
		maxz = max( zs ) 
		minz =  min( zs ) 
		
		correctsubtilt( options, subtilts, angles, ctfs, apix, nangles, nimgs, framexsize, icethickness, maxz, minz )
	
	'''
	#CTF correction for an entire tilt series is different than for a subtiltseries of a subtomogram,
	at least when using the same correction for the entire image. In the latter case (subtiltseries)
	the images are only accepted as a stack, and you want to account for differences in z-height.
	'''	
		
		
		
		
		
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
			print '''ERROR: (e2spt_ctf.py)(main) ERROR: It seems like you have fewer parameter
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


def errordetector( options ):
	
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
	
	
	return



def procsubtiltsdir():
	
	if not options.framexsize:
		print "ERROR: provide a value for options.framexsize if processing --subtiltsdir"
		sys.exit(1)
		
	findir = os.listdir( options.subtiltsdir )
	
	nimgs = 0
	
	linesxz = []
	linesyz = []
	maxfiles = options.subset
	fn = 0
	for f in findir:
		
		if options.subset:
			if fn == options.subset:
				break
				
		if '.hdf' in f:
			stsfile = options.subtiltsdir + '/' + f
			
			
			subtilts.append( stsfile )
			print "\nFound subtiltseries", f
			
			
			stshdr = EMData( stsfile, 0, True )
			hdrcoords = stshdr['ptcl_source_coord']
			ptclnx = stshdr['nx']
			print "hdrcoords are", hdrcoords
			x = hdrcoords[0]
			y = hdrcoords[1]
			z = hdrcoords[2]
			
			xs.append( x )
			ys.append( y )
			zs.append( z )
			
			linexz = str(x) + ' ' + str(z) + '\n'
			linesxz.append( linexz )
			
			lineyz = str(y) + ' ' + str(z) + '\n'
			linesyz.append( lineyz )
			
			fn+=1
					
	xzcoords = options.path + '/xz_distribution.txt'
	f = open( xzcoords, 'w')
	f.writelines( linesxz )
	f.close()
	
	xlabel = 'X (pixels)'
	zlabel = 'Z (pixels)'
	
	plotnamexz = options.path + '/xz_distribution.png'	
	title = 'X vs Z'
	
	if xs and zs:
		generalplotter( options, xs, zs, xlabel, zlabel, plotnamexz, title, flipyaxis=False, fit=True, ptclnx=0 )
	else:
		print "\nERROR: cannot plot x vs z values because arrays are empty", xs, zs
	
	
	yzcoords = options.path + '/yz_distribution.txt'
	f = open( yzcoords, 'w')
	f.writelines( linesyz )
	f.close()
	
	ylabel = 'Y (pixels)'

	plotnameyz = options.path + '/yz_distribution.png'	
	title = 'Y vs Z'
		
	if ys and zs:
		generalplotter( options, ys, zs, ylabel, zlabel, plotnameyz, title, flipyaxis=False, fit=True )
	else:
		print "\nERROR: cannot plot y vs z values because arrays are empty", ys, zs
	
	
	nimgs = EMUtil.get_image_count( subtilts[0] )
	apix = EMData( subtilts[0], 0, True)['apix_x']
		
	return



'''
c:function to run commands and the command line
'''
def runcmd( options, cmd ):
	if options.verbose > 9:
		print "\n(e2spt_autoboxer)(runcmd) Running command", cmd
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	if options.verbose > 9:
		print "\n(e2spt_autoboxer)(runcmd) Done"
	return


def correctsubtilt( options, subtilts, angles, ctfs, apix, nangles, nimgs, framexsize, icethickness, maxz, minz ):
	
	ii = 0
	globalAvgDefErrors=[]
	
	zfillfactor = len( str( len(subtilts) ) )
	
	#indxstoexclude=[int(i) for i in options.exclude.split(',')]
	
	#ss=0
	for sts in subtilts:
		imghdr = EMData( sts, 0, True )
		
		#print " icethickness is", icethickness
		#print "ice < int(imghdr['nx'])", int(icethickness) < int(imghdr['nx'])
		#print "ice and type are", type(options.icethickness), icethickness
		#print "imghdrnx and type are",type(imghdr['nx']),imghdr['nx']
		
		if icethickness and int(icethickness) < int(imghdr['nx']):
			print """\nWARNING: the ice must usually be thick enough to contain a layer of molecules and will be set
			to half the X and Y sides of the subtomograms. It must be >= than %d pixels.""" %( imghdr['nx'] )
			#sys.exit()
			
		coords = imghdr['ptcl_source_coord']
		coordx = coords[0]
		coordz =  coords[-1]
		nx = imghdr['nx']
		
		if options.verbose:
			print "Fixing subtomogram", ii
		
		n = EMUtil.get_image_count( sts )
		#print "\nprocessing subtiltseries %d, %s, with n=%d images in it, and will exclude %d"%(ii, sts,n,len(indxstoexclude))
		
		print "\nprocessing subtiltseries %d, %s, with n=%d images in it"%(ii, sts,n)

		if n!= nangles:
			print """WARNING: The number of angles %d does not coincide with number of images %d. 
			However, the actual angles being used should be obtained directly from the 2d image's header""" % ( nangles,  n )
			#sys.exit(1)
		
		flippedsts = options.path + '/' + os.path.basename( sts ).replace('.hdf','_phflip.hdf')
		if options.outputstem:
			flippedsts = options.path + '/' + options.outputstem + 'ptcl' + str(ii).zfill( zfillfactor ) + '_phflip.hdf'
		
		phfimgs = []
		defocuserrors=[]
		checkerrors = 0
		
		print "angles are", angles
		
		for m in range( n ):
		
			print "img number %d/%d" %(m,n)
			#if m not in indxstoexclude:
			#	print "m=%d is NOT in indxstoexclude", m
			img = EMData( sts, m )
			img['xform.align3d'] = Transform()
		
			angle = round(img['spt_tiltangle'],2)
		
			#angle2 = round(angles[ m ],2)
		
			#if angle != angle2:
			#	print "ERROR: The angle in the particle's header %.4f does not match the one in the angle's list"
		
			img['xform.projection'] = Transform({"type":"eman","az":90.0,"alt":float(angle),"phi":-90.0,'tx':0,'ty':0,'tz':0}) 
		
			#Multiple alt * -1 since EMAN2's convention for positive altitude is backwards from IMOD *NOT TRUE
			#,"weight":1.0}
		
			ctf = ctfs[ angle ]
			print "\n\nUncorrected defocus is", ctf.defocus
		
			'''
			For positive tilt angles (counter clockwise) the defocus decreases (the particle is more overfocus, less defocused) for positions px right of the tilt axis
			while defocus increases for particles left of the tilt axis (they are more defocused).
			For negative tilt angles (clockwise) the defocuses increases (the particle is more defocused)for px right of the tilt axis while
			defocus decreases (more overfocused, less defocused) for particles left of the tilt axis.
			'''
			px = ( coordx - framexsize/2.0 ) * apix/10000
			if px < 0:
				print "\npx (in microns) is left of the tilt axis", px
			elif px > 0:
				print "\npx (in microns) is right of the tilt axis", px
			elif px==0:
				print "\npx (in microns) is on the tilt axis", px
			
			dzx = -1 * px * numpy.sin( math.radians( angle ) )		#the -1 accounts for the fact that positive tilt angles are clockwise, negative counter clockwise
		
			if angle < 0.0:
				print "\ngiven a negative, CLOCKWISE tilt angle=%f, and coordx=%f pixels, px=%f microns, THEN dzx=%f microns" %( angle,coordx,px,dzx) 
			if angle > 0.0:
				print "\ngiven a positive, COUNTER CLOCKWISE tilt angle=%f, and coordx=%f pixels, px=%f microns, THEN dzx=%f microns" %( angle,coordx,px,dzx) 

		
			newdefocus = ctf.defocus + dzx 
			print "\ntherefore, for angle=%f, and defocus=%f, the first corrected defocus is NEWdefocus1=%f" % ( angle, ctf.defocus, newdefocus )
		
			pz = 0
			dzz = 0
			if options.defocustop:
				relativecoordz = coordz - maxz
				pz = relativecoordz * apix/10000
				#pz = ( coordz - 2.0*icethickness/2.0 ) * apix/10000
		
			elif options.defocusbottom:
				relativecoordz = coordz - minz
				pz = relativecoordz * apix/10000
				#pz = ( coordz + 2.0*icethickness/2.0 ) * apix/10000
		
			elif options.nozcorrection:
				pass

			else:														#assume the defocus signal comes from the middle
				middle = ( maxz + minz ) / 2.0
				relativecoordz = coordz - middle 
				pz = relativecoordz * apix/10000
				#pz = ( coordz - icethickness/2.0 ) * apix/10000
		
			if not options.nozcorrection:	
				dzz = -1 * pz * numpy.cos( math.radians( angle ) )		#for negative positions pz, particles are MORE defocus due to their depth on the ice;
				newdefocus += dzz										#or positive positions pz, particles are LESS defocused. 
																		#at tilt angle 0, the contribution of dzz is equal to pz (in magnitude) in microns, yet opposite in sign because defocus is defined as positive.
				print "dzz=%f applied to defocus=%f, therefore NEWdefocus2=%f" %(dzz,ctf.defocus,newdefocus)			#at tilt angle 90 (hypothetical) the contribution would be irrelevant, 0.
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
				checkerrors = 1
			except:
				pass
		
			ret = phaseflipper( options,img,finalctf )
			imgflipped = ret[0]
			imgflipped['ctf'] = finalctf
		
			#print "Flipped outstack to write is", flippedsts
			#print "imgflipped and type are", imgflipped, type(imgflipped)
			#print "index to write is", m
		
			print "received from flipper and will write to stack", imgflipped['minimum'],imgflipped['maximum'],imgflipped['sigma'],imgflipped['mean']
				
			if options.save2d:
			
				imgflipped.write_image( flippedsts, -1 )	
		
		
			phfimgs.append( imgflipped )
			#else:
			#	print "image %d being excluded because indx in indxstoexclude", m
		
		print "sending these many phfimgs", len(phfimgs)
		rec = reconstruct3d( options, phfimgs, apix )
		
		if defocuserrors:
			defocuserrorsAvg=sum(defocuserrors)/len(defocuserrors)
			rec['spt_avgDefocusError']=defocuserrorsAvg
			
			globalAvgDefErrors.append(defocuserrorsAvg)
		elif checkerrors:
			print "Defocus errors is empty!", defocuserrors
			sys.exit()
		
		if options.save3d:
			#sts3d = options.path + '/' + os.path.basename( sts ).replace('.hdf','_PHFLIP3D.hdf')
			stack3d = options.path + '/stack_phflip3d.hdf'
			if options.outputstem:
				stack3d = options.path + '/' + options.outputstem +'_phflip3d.hdf'
			
			if options.invert:
				stack3d.replace('.hdf','_inv.hdf')

			rec.write_image( stack3d , ii )
			
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
	
	print "phfimgs len is", len(phfimgs)
	print "phfimgs are", phfimgs
	box = phfimgs[-1]['nx']
	
	originalboxsize = box
	
	print "in reconstruct3d"
	
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
		print "--reconstructor, options.reconstructor is", options.reconstructor
		print "its len is", len(options.reconstructor)
		if len(options.reconstructor) > 1:
			if 'mode' in options.reconstructor[-1]:
				mode=options.reconstructor[-1]['mode']
				
				print "\nThe reconstructor mode has been changed from default to", mode
				#sys.exit()
	
	print "\Boxsize to reconstruct, after padding, is", box
	
	print "reconstructor is", options.reconstructor
	box = int(box)				
	r = Reconstructors.get(options.reconstructor[0],{'size':(box,box,box),'sym':'c1','verbose':True,'mode':mode})
	
	print "r is", r
	
	r.setup()

	k=0
	weight = 1.0
	for p in phfimgs:
		print "Adding projection k", k
		print "Whose min and max are", p['minimum'], p['maximum']
		print "The size of the prj to insert is", p['nx']
		
		pc=p.copy()
		
		if options.pad2d:
			pc = clip2D( pc, box )
		
		print "\n\n\n\nPPPPP the projection direction is", pc['xform.projection']
		
		pm = r.preprocess_slice(pc,pc['xform.projection'])
		r.insert_slice(pm,pm['xform.projection'],weight)
		k+=1
	
	rec = r.finish(True)

	rec['apix_x']=apix
	rec['apix_y']=apix
	rec['apix_z']=apix
	rec['origin_x']=0
	rec['origin_y']=0
	rec['origin_z']=0
	
	recfinal = clip3D( rec, originalboxsize )
	print "recfinal, min, max, sigma, mean", rec['minimum'],rec['maximum'],rec['sigma'],rec['mean']

	'''
	Preserve SPT parameters in header
	'''
	names = phfimgs[-1].get_attr_dict()
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
	

def clip3D( vol, sizex, sizey=0, sizez=0 ):
	
	if not sizey:
		sizey=sizex
	
	if not sizez:
		sizez=sizex
	
	volxc = vol['nx']/2
	volyc = vol['ny']/2
	volzc = vol['nz']/2
	
	Rvol =  Region( (2*volxc - sizex)/2, (2*volyc - sizey)/2, (2*volzc - sizez)/2, sizex , sizey , sizez)
	vol.clip_inplace( Rvol )
	#vol.process_inplace('mask.sharp',{'outer_radius':-1})
	
	return vol



def genctfparamlines( options, apix, nimgs, angles, imagefilenames, icethickness=0 ):
	
	#indxstoexclude = [int(i) for i in options.exclude.split(',')]
	
	print "\n\n\n\n\n\ne2spt_ctf (genctfparamlines)"
	print "received imagefilenames", imagefilenames, len(imagefilenames)
	
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
			#if kk not in indxstoexclude:
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
			#if kk not in indxstoexclude:
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
		#if kk not in indxstoexclude: 
		if len( defoci ) != nimgs:
			print """ERROR: The number lines in the file provided
				through --defocilist should match the number of images in the tiltseries
				or in each subtiltseries provided.
				To use the same parameters for all images provide them
				explicitly through --cs, --bfactor,--voltage, --defocus, --ampcont and, optionally --apix
				(this apix will be read from the header if not provided, so make sure it's correct). 
				"""
			sys.exit(1)
		else:
			print "same number of defoci %d as imgs %d" %(len(defoci),nimgs)
			if len(defoci) != len(angles):
				print "ERROR: number of defoci %d not the same as number of angles %d" %(len(defoci),len(angles))
				sys.exit()
		
		if options.voltage and options.cs and apix and options.bfactor and options.ampcont:
			print """\nExplicit parameters --cs,--apix,--bfactor,--voltage and --ampcont 
				will be used."""
				
			kk=0
			for d in defoci:
				#if kk not in indxstoexclude and kk < len(angles):
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
		
	elif options.autofit:
		print "\n\n\nautofitting using voltage=%.2f, cs=%.2f, apix=%.2f, ampcont=%.2f" %( float(options.voltage), float( options.cs), float( apix), float( options.ampcont))
		print "imagefilenames", imagefilenames
		if options.voltage and options.cs and apix and options.ampcont:
			ctfs = sptctffit( options, apix, imagefilenames, angles, icethickness )
			
	
	else:
		if options.voltage and options.cs and options.defocus and apix and options.bfactor and options.ampcont:
			print """\nExplicit parameters --defocus,--cs,--apix,--bfactor,--voltage and 
				--ampcont will be used."""
			
			print '\ndefocus to set is', options.defocus	
			line = 'defocus=' + str( options.defocus ) + 'ampcont=' + str( options.ampcont ) + 'voltage=' + str( options.voltage ) + 'cs=' + str( options.cs ) + 'apix=' + str( apix ) + 'bfactor=' + str( options.bfactor )
			ctf = ctfparamparser( line )
			print '\nafter parsing and making ctf object, set defocus is', ctf.defocus
			for kk in range( len(angles) ):
				#if kk not in indxstoexclude:
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
	
	finalangles = list( angles )
	if options.exclude:
		excludeindxs = options.exclude.split( ',' )
		excludeindxsint = [ int(i) for i in excludeindxs ]
	
		finalangles = [ angle for i, angle in enumerate( angles ) if i not in excludeindxsint ]
		excludedangles = [ angle for i, angle in enumerate( angles ) if i in excludeindxsint ]
		
		print "\n %d angles were excluded" %( len(excludeindxsint) )
		print excludedangles
		print "\ntherefore, final included angles are", finalangles
	
	return finalangles


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
	
	print "In phase flipper, PRJ min, max, sigma, mean", prj['minimum'],prj['maximum'],prj['sigma'],prj['mean']
	
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
	print "ctf to apply is", ctf
	print "to prj_fft min, max, sigma, mean", prj_fft['minimum'],prj_fft['maximum'],prj_fft['sigma'],prj_fft['mean']
	ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)

	prj_fft.mult(flipim)

	intermediate = prj_fft.copy()

	prj_flipped=prj_fft.do_ift()
	
	
	#print "prj_flipped to return, min, max, sigma, mean", prj_flipped['minimum'],prj_flipped['maximum'],prj_flipped['sigma'],prj_flipped['mean']

	return prj_flipped, intermediate
	
	
def fullcorrection(options,img,ctf):
	
	prj=img.copy() 
	
	prj_fft = prj.do_fft()

	flipim = prj_fft.copy()	
	ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_TOTAL)

	prj_fft.mult(flipim)

	intermediate = prj_fft.copy()

	prj_flipped=prj_fft.do_ift()
	
	return prj_flipped, intermediate
	
	
	
def tilerfft(options, angle, imgt, currentstrip, nstrips, start, end, step, savestriptiles, saveffts):	
	kk = 0
	fftcumulative=None	
	nbx=0
	
	ny = imgt['ny']
	
	
	signtag = 'p'
	if angle < 0.0:
		signtag ='m'
	step = int(step)
	if options.tileoverlap:
		step = int(options.tilesize/2)
		print "tiles are overlapping, therefore step is", step
	#for x in range( micrographstarts[m], micrographstarts[m] + micrographwidth - options.tilesize + 1, options.stripstep ):
	for x in range( start, end, step ): 
		#print "\nsumming over tiles along y"
		#for y in range(0, ny - options.tilesize+1, options.tilesize):
		for y in range(0, ny - options.tilesize+1, step):
			#print "tile at y", y
			clipr = imgt.get_clip(Region(x,y, options.tilesize, options.tilesize))
			
			if clipr['sigma']:
				
				allgood = 1
				
				if float(options.prunetest) >= 0.0:
					allgood = checkcorners( clipr, options )
				
				if allgood:			
				
					#fftcumulativeimgfile = options.path  + '/angle_' + signtag + str( int(math.fabs( round(angle) ) + '_strip' + str(currentstrip).zfill(len(str(nstrips))) + '_fft.hdf'

								
					clipr.process_inplace("normalize.edgemean")
					
					if savestriptiles:
						clipout = options.path + '/angle_' + signtag + str( int(math.fabs( round(angle) ))) + '_strip' + str(currentstrip).zfill(len(str(nstrips))) + '.hdf'
						clipr.write_image( clipout, kk )
		
					kk+=1
			
					fft = clipr.do_fft()
					fft.ri2inten()
					if fftcumulative==None: 
						fftcumulative=fft
					else: 
						fftcumulative+=fft
					nbx+=1
				else:
					"WARNING: tile excluded because it has at least one bad corner!"
				
			else:
				print "WARNING: tile excluded because sigma is zero!"

	if nbx > options.mintiles:
	
		if fftcumulative:
			fftcumulative.mult(1.0/(nbx*options.tilesize**2))
			fftcumulative.process_inplace("math.sqrt")
			fftcumulative["is_intensity"]=0				# These 2 steps are done so the 2-D display of the FFT looks better. Things would still work properly in 1-D without it
	
			signtag = 'p'
			if angle < 0.0:
				signtag ='m'
	
			#plotname = options.path + '/fit_' + str( imgindx ).zfill( len(str( nangles ))) + '_' + signtag + str( int(math.fabs( round(angle) ) )).zfill(3) +'.png'
	
			if saveffts:
				fftcumulativeimgfile = options.path  + '/angle_' + signtag + str( int(math.fabs( round(angle) ))) + '_strip' + str(currentstrip).zfill(len(str(nstrips))) + '_fft.hdf'
				fftcumulative.write_image(fftcumulativeimgfile,0)	
	
			return fftcumulative
		else:
			print "\nWARNING: bad strip!"
			return None
	else:
		print "\nWARNING: strip has too few good tiles and thus is being skipped"
		return None

def checkcorners( img, options ):
	
	nx = img['nx']
	#ny = img['ny']
	
	cornernx = round( nx * options.prunetest)
	
	if not cornernx:
		cornernx = 4
	#cornerny
	#clipr = imgt.get_clip(Region(x,y, options.tilesize, options.tilesize))
	
	corner1 = img.get_clip( Region(0,0, cornernx, cornernx))
	corner2 = img.get_clip( Region(nx-cornernx,0, cornernx, cornernx))
	corner3 = img.get_clip( Region(nx-cornernx,nx-cornernx, cornernx, cornernx))
	corner4 = img.get_clip( Region(0,nx-cornernx, cornernx, cornernx))
	
	if not corner1['sigma']:
		return 0
	elif not corner2['sigma']:
		return 0
	elif not corner3['sigma']:
		return 0
	elif not corner4['sigma']:
		return 0
	else:
		return 1
	

def fitdefocus( ffta, angle, apix, options, nsubmicros, currentsubmicro, defocusmin, defocusmax, defocusstep, x=0 ):


	
	fftbg = ffta.process("math.nonconvex")
	fft1d = ffta.calc_radial_dist(ffta.get_ysize()/2,0.0,1.0,1)	# note that this handles the ri2inten averages properly
	
	#print "fft1d is", fft1d, type(fft1d)
	
	signtag = 'p'
	if angle < 0.0:
		signtag ='m'
	
	
	if 0:
		fft1doutfile = options.path + '/angle_' + signtag + str( int(math.fabs( round(angle) ))) + '_strip' + str(currentsubmicro).zfill(len(str(nsubmicros))) + '_fft1d.txt'
		g=open( fft1doutfile, 'w' )	
		fft1dstr = [ str(ii)+ ' ' + str(fft1d[ii]) + '\n' for ii in range(len(fft1d)) ]
		fft1dstrf = fft1dstr[1:-1]
	
		#print "type fft1dstrf", type(fft1dstrf), fft1dstrf
		g.writelines( fft1dstrf )
		g.close()

	
	ctf = None
	try:
		# Compute 1-D curve and background
		ds = 1.0/( options.apix * options.tilesize )
		bg_1d = e2ctf.low_bg_curve(fft1d,ds)
		
		#initial fit, background adjustment, refine fit, final background adjustment
		#ctf = e2ctf.ctf_fit(fft1d,bg_1d,bg_1d,ffta,fftbg,options.voltage,options.cs,options.ampcont,apix,1,dfhint=( defocusmingradient, defocusmaxgradient, step ))
	

		ctf = e2ctf.ctf_fit( fft1d, bg_1d, bg_1d, ffta, fftbg, options.voltage, options.cs, options.ampcont, apix, 1,dfhint=( defocusmin, defocusmax, defocusstep ) )
		bgAdj(ctf,fft1d)
	except:
		print "ctf fit failed! first try"
		print "len fft1d is", len(fft1d)
		print "ffta is", ffta
		ctf = None
			
	try:
		#ctf = e2ctf.ctf_fit(fft1d,ctf.background,ctf.background,ffta,fftbg,options.voltage,options.cs,options.ampcont,apix,1,dfhint=( defocusmingradient, defocusmaxgradient, step))
		ctf = e2ctf.ctf_fit(fft1d,ctf.background,ctf.background,ffta,fftbg,options.voltage,options.cs,options.ampcont,apix,1,dfhint=( defocusmin, defocusmax, defocusstep ))	
		bgAdj(ctf,fft1d)

		#if options.astigmatism: 
		#	e2ctf.ctf_fit_stig(ffta,fftbg,ctf)
	except:
		print "ctf fit failed! second adjustment try"
		print "len fft1d is", len(fft1d)
		print "ffta is", ffta
		ctf = None
		

	if 0:
		'''
		Plot background subtracted curve?
		'''
		fz=int(ctf.zero(0)/(ds*2)) 	#jesus
		bs =[ str(i) + ' ' + str(fft1d[i]-ctf.background[i]) for i in xrange(fz)] #jesus
	
		bsoutfile = options.path + '/angle_' + signtag + str( int(math.fabs( round(angle) ))) + '_strip' + str(currentsubmicro).zfill(len(str(nsubmicros))) + '_bs1d.txt'
	
		gg=open( bsoutfile, 'w' )
	
		bsstr = [ str(ii)+ ' ' + str( bs[ii] ) + '\n' for ii in range(len( bs )) ]
		#bsstrf = fft1dstr[1:-1]
	
		#print "type fft1dstrf", type(fft1dstrf), fft1dstrf
		gg.writelines( bsstr )
		gg.close()
	

		'''
		'''
	
	stripdefocus = None
	if ctf:
		ctf.background = bg_1d
		ctf.dsbg = ds
		
		ctf1d = ctf.compute_1d
		#ctfbs = ctf1d - bg_1d
		
		#print '\n\nctf1d is', ctf1d
		#print '\n\nctfbs is', ctfbs
		
		#c = open( options.path + '/ctf.txt', 'w' )
		#c.writelines( ctf1d )
		#c.close()
		
		#print "\nds is", ds

		#db=js_open_dict(info_name(arg,nodir=not options.usefoldername))
		#db["ctf_frame"]=[512,ctf,(256,256),set(),5,1]

		#print info_name(arg,nodir=not options.usefoldername),ctf
	
		#print "\nctf and type", ctf, type(ctf)
		#print "ctf atrributes", ctf.get_attr_dir()
		#print "ctf dir", dir(ctf)
	
		stripdefocus = ctf.defocus
		
		 
		
	
	#print "background is", ctf.background

	return stripdefocus

	
	
def sptctffit( options, apix, imagefilenames, angles, icethickness ):
	
	print "e2spt_ctf (sptctffit)"
	
	#if options.invert: cmd += " --mult=-1"
	#if options.edgenorm: cmd += " --process=normalize.edgemean"
	#if options.xraypixel: cmd += " --process=threshold.clampminmax.nsigma:nsigma=4"
		
	# We estimate the defocus and B-factor (no astigmatism) from the micrograph and store it in info and the header
	
	#d=EMData(output,0)
	
	#targetdefocus = 3
	#defocusmin = 3
	#defocusmax = 3
	
	#if options.defocus:
	#targetdefocus = options.defocus
	
	defocusmin = 0.5
	defocusmax = 15
	defocusstep = 0.1
	defocuswiggle = None	
	
	ctfs={}
	print "\n\nimagefilenames", imagefilenames, len(imagefilenames)
	
	#angles.sort()	#use to get image index
	
	
	allglobaldefocuses = {}
	angerrors = {}
	#angerrorlines = []
	imgnum = 0
	
	maxangle = max( [ math.fabs( a ) for a in angles] )
	skip = []
	if options.skipstripping:
		skip = [ int(i) for i in options.skipstripping.split(',') ]
	
	
	anglestoexclude=[]

	'''
	calculate the maximum allowable depth of focus or variation in defocus across particles to achieve 2/3 Nyquist resolution
	'''
	nyquist = apix*2.0
	twothirdsnyquist = nyquist*3.0/2.0

	#set the wavelength lambd of the electrons depending on acceleration voltage
	lambd = 0.0197
	if int(options.voltage) == 200:
		lambd = 0.0251
	elif int(options.voltage) == 100:
		lambd = 0.037

	allowabledz = twothirdsnyquist*twothirdsnyquist/(2.0*lambd)

	print "\n\n\n\n\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nthe theoretical allowable depth of focus or defocus variation limit in angstroms is %.2f to reach 2/3 nyquist resolution %.2f" %( allowabledz, twothirdsnyquist )

	if options.defocusvariationlimit:
		allowabledz = options.defocusvariationlimit
		print "however, the allowable depth of focus in micrometers has been set via --defocusvariationlimit and is", allowabledz
		print "which in angstroms is", allowabledz*10000

	
	progresscount=1
	for label in imagefilenames:
		angle = imagefilenames[label][1]
		imgindx = angles.index( angle )
		
		print "\n\n\n\n\n\n\nautofitting ctf for stem %s, image %s, imgindx %d, angle %f, apix %f, progress = %d/%d" %( label, imagefilenames[label][0],imgindx,angle,apix,progresscount,len(angles))
		img = EMData( imagefilenames[label][0], 0 )
		progresscount+=1

		nx = int( img['nx'] )
		ny = int( img['ny'] )
		
		
		globaldefocus = None
		globalmiddle = None
		centerdefocus = None
		defocuserror=0
		if options.firstfitglobal:
			'''
			First fit by tiling the ENTIRE micrograph to get the average defocus without gradient compensation,
			just to get a first 'aproximate' estimate for the range of defocuses to look for
			'''
			#options, angle, imgt, currentstrip, nstrips, start, end, step, savestriptiles, saveffts
			fftg = tilerfft( options, angle, img, 0, 1, 0, nx - options.tilesize + 1, options.tilesize, 0, 0 )
			
			defocusmin = 0.5
			defocusmax = 15.0

			if options.targetdefocus:
				defocusmin = options.targetdefocus -1.5
				defocusmax = options.targetdefocus + 1.5

			#ffta, angle, apix, options, nsubmicros, currentsubmicro, defocusmin, defocusmax, defocusstep, x=0 ):
			globaldefocus = fitdefocus( fftg, angle, apix, options, 1, 0, defocusmin, defocusmax, defocusstep, 0)
			
			if globaldefocus:
				defocusmin = globaldefocus -1.5
				defocusmax = globaldefocus + 1.5


			globalmiddle = ( nx/2.0 ) * apix / 10000	#convert to micrometers, so x and y axis are in the same units
			
			
			
			
			#get central region of 2 tile width
			#r = Region( img['nx']/2 - options.tilesize - int(ceil(options.tilesize/2.0)), 0, 3*options.tilesize, img['ny'])
			
			r = Region( img['nx']/2 - options.tilesize, 0, 2*options.tilesize, img['ny'])

			
			fixcenter = 0
			
			if fixcenter:
				imgcenterstrips = img.get_clip( r )
			
				centerstripsnx = imgcenterstrips['nx']
				print "centerstripsnx", centerstripsnx
				imgcenterstrips_nx = imgcenterstrips['nx']
			
				fftcenter = tilerfft( options, angle, imgcenterstrips, 0, 1, 0, imgcenterstrips_nx - options.tilesize + 1, options.tilesize, 0, 0 )
				
				centerdefocus = fitdefocus( fftcenter, angle, apix, options, 1, 0, defocusmin, defocusmax, defocusstep, 0)
			
				print "centerdefocus %f , globaldefocus %f, for angle %f" % ( centerdefocus, globaldefocus, angle )
			
			
			if globaldefocus and not centerdefocus:
				print "\nglobal defocus for this image is", globaldefocus
				defocusmin = globaldefocus - 1.5
				defocusmax = globaldefocus + 1.5	
				
				allglobaldefocuses.update({ angle:globaldefocus} )
			
			elif globaldefocus and centerdefocus:
				defocuserror = globaldefocus-centerdefocus
			
				if defocuserror:
					averagedefocus = (globaldefocus+centerdefocus)/2.0
					globaldefocus = averagedefocus
					
					#defocusmin = averagedefocus - defocuswiggle/2.0
					#defocusmax = averagedefocus + defocuswiggle/2.0	
					
					print "averagedefocus",averagedefocus
					#print "defocusmin after wiggle",defocusmin
					#print "defocusmax after wiggle",defocusmax
					
					allglobaldefocuses.update({ angle:averagedefocus} )
				else:
					#defocusmin = globaldefocus - defocuswiggle/2.0
					#defocusmax = globaldefocus + defocuswiggle/2.0	
					
					#print "defocusmin after wiggle",defocusmin
					#print "defocusmax after wiggle",defocusmax
					
	
					allglobaldefocuses.update({ angle:globaldefocus} )
					
			else:
				print "\nWARNING! global defocus fitting failed! for image", label
		
		
		print "\n\nfound this globaldefocus %.2f for label %s" %(globaldefocus,label)

		xs = []
		imgdefocuses = []
		micrographmiddle = img['nx']/2
		
		#do not strip-fit if an image is bad
		if imgindx not in skip:
			
		
			imgdefocuses = []
		
			faileddefs = []
			failedmids = []
		
			'''#
			#Position of the tilt axis in micrometers
			#'''
			pxta = ( nx/2.0 ) * apix/10000
		
			'''#
			#Find the distance dx100dz away from the tilt axis, across which the vertical distance 
			dzx changes by 100 nm, i.e., 0.1 micrometers, given the specified ice thickness
			#'''
		
			#deltata = pxta
			
			#dx100dz = nx/2.0						#at zero degrees dzx100 is half the micrograph size, since there should be no defocus variation due to tilt
			
			stripwidth = nx 						#at zero degrees dzx100 is the micrograph size, since there should be no defocus variation due to tilt


			if math.fabs( angle ) > 0.5:
		
				#defocusvariationlimit is in micrometers. it is the variation in defocus to tolerate across a strip due to tilt and still consider the defocus "constant"
			
				#dx100dz = math.fabs( options.defocusvariationlimit / numpy.sin( math.radians( angle ) ) * 10000/apix )	#calculate strip half-width in pixels
				
				stripwidthold = math.fabs( 0.1 / numpy.sin( math.radians( angle ) ) * 10000/apix )

				if icethickness:
					depthdefocus = math.fabs( icethickness / numpy.cos( math.radians( angle ) ) )	#icethickness, in pixels, will contribute to defocus variation with tilt, in addition to the tilting itself
				
					stripwidthold += depthdefocus				
				
				stripwidthold*=2 #this factor of two accounts for the fact that the old method only calculated things for half of the micrograph (from the tilt axis)

				print "old strip width in pixels was", stripwidthold

				#cos/sin is cot, but the function does not exist in numpy
				stripwidth = int( math.fabs( allowabledz * (numpy.cos( math.radians (angle) ) / numpy.sin( math.radians (angle) ) ) / apix - icethickness / numpy.sin( math.radians( angle )) ) )#this is in pixels

				print "however, new strip width is", stripwidth

				if stripwidth > nx:
					stripwidth = nx
		
			'''#
			#submicrograph width (the region at "equal/constant defocus" to tile) will be twice dx100dz (stripwidthold) or just stripwidth 
			#'''
			micrographwidth = int( stripwidth )
		
			#micrographwidth = 3*options.tilesize
		
		
			adjuststart = 0								#Flag parameter used (below) for highly tilted images in which the region of "constant defocus" is smaller than the tilesize
			if micrographwidth < options.tilesize:
				options.stripstep = micrographwidth
				micrographwidth = options.tilesize
				adjuststart = 1
		
			print "for angle %f, micrographwidth is %d" %(angle, micrographwidth)
		
			#micrographsboundaries = []
		
			'''#
			#Find out how many submicrographs of the determined width fit in the whole image.
			Since during tiling these regions need to be centered at the tilt axis, any even
			number of micrographs needs to be rounded down to the closest smallest odd number.
			For example, if submicrograph width is 10 and the whole image has length 20, you 
			can fit two submicrographs in the whole image, yes; but, since the first submicrograph
			is to be centered at pixel 10, it will go from pixel 5 to 15; thus, the remaning
			5 pixels at either side cannot fit a whole submicrograph. Therefore, you find the 
			largest odd number of micrographs that fit in the whole image and treat the rest 
			as excedent. 
			#'''
			nmicros = nx/micrographwidth
			nmicrosint = int( math.floor( nmicros ) )	
			if not nmicrosint % 2:
				nmicrosint -= 1
		
			print "nmicrosint is %d since (nmicrosint+2)*micrographwidth is %d while nx is %d" %( nmicrosint, (nmicrosint+2)*micrographwidth, nx )
		
			excedent = nx - nmicrosint * micrographwidth
		
			print "excendent is", excedent
		
			#excedent = nmicros - nmicrosint
		
			#nmicrosroundedup = math.ceil( nmicros )
		
			micrographstarts = []
		
			aux = 0
			#if excedent >= options.tilesize:	
			if excedent:
				if excedent >= options.tilesize*2 and not options.excludeedges:
					aux = 2
		
			print "\nadjuststart is", adjuststart

			include = 1
			if not adjuststart:
				for iz in range( nmicrosint + aux ):		#You'll usually have two extra micrographs to include the excedent at either side of the central region where you can fit an exact odd multiple of submicrograph widths
					print "iz is", iz
				
					if aux:
				
						if iz == 0:							#first submicrograph starting point
							start = 0
						elif iz == nmicrosint + 1:				#last submicrograph starting point; range goes from 0 to microsint + 2 at most, without including the upper bound
							start = nx - micrographwidth
						else:
							start = (iz-1)*micrographwidth + excedent/2
						
					else:
						start = (iz)*micrographwidth + excedent/2		#Don't bother with edge micrographs if aux is 0
					
					print "withOUT adjuststart, start to append is", start
				
				
					#if include:
				
					micrographstarts.append( int(start) )
				
			else:
				nmicrosint = int( math.floor( ( nx - options.tilesize ) / options.stripstep ) )
			
				excedentnew = nx - options.tilesize * nmicrosint
			
				init = 0
				if options.excludeedges:
					init = 1
			
				for ii in xrange(init,nmicrosint):
				
					#if int(excedent)/2 < int(options.stripstep):
					#	start = i*options.stripstep
					#else:
				
					start = ii*options.stripstep
				
					print "with adjuststart, start to append is", start
				
					micrographstarts.append( int(start) )
			
				if excedentnew and not options.excludeedges: 	#plus add the final one if there's excedentnew
					start = nx - options.tilesize
					micrographstarts.append( start )
				
			micrographstarts.sort()
			print "\nfor img %d micrographstarts are" % (imgindx) 
			print micrographstarts	
		
			#micromids = [h+micrographwidth/2 for h in micrographstarts]
			micromids = []
		
		
		
			for m in range( len (micrographstarts)):
			

			
				'''#
				#The defocus is fitted per submicrograph (regions of pseudo-constant defocus); therefore, things are reset here
				#'''
			
				#print "m and type", m, type(m)
				#print "submicrographwidth, corresponding to a region of pseudo constant defocus is", micrographwidth
				print "options.tilesize and type", options.tilesize, type( options.tilesize )
			
			
			
				fftc = tilerfft( options, angle, img, m, len(micrographstarts), micrographstarts[m], micrographstarts[m] + micrographwidth - options.tilesize + 1, options.stripstep, options.savestriptiles, options.saveffts )
				#tilerfft(options, img, start, end, step)
				#fft( micrographstarts[m], micrographstarts[m] + micrographwidth - options.tilesize + 1, options.stripstep ):
	
				
				micrographmiddle =  ( list(micrographstarts)[m] + micrographwidth/2 ) * apix / 10000.00 #xaxis in micrometers too

				stripdefocus = None
				if fftc:
					
					if icethickness and options.coords:
						icethicknessm = icethickness * img['apix_x']/10000
					
						defocuswiggley = math.fabs( icethicknessm / math.cos( math.radians(angle) ))
						if defocuserror:
							defocuswiggley += defocuserror
						
						defocuswigglex = 2*math.fabs( -1*(micrographmiddle - img['nx']/2.0)*math.sin( math.radians(angle) ) * img['apix_x']/10000)
						
						defocuswiggle = defocuswigglex + defocuswiggley
						
						print "\nicethicknessm is", icethicknessm
						print "\ndefocuswiggle is", defocuswiggle
						
						if defocuswiggle:
							defocusmin = globaldefocus - 2*defocuswiggle
							defocusmax = globaldefocus + 2*defocuswiggle
							
							print "after wiggle, defocusmin is", defocusmin
							print "after wiggle, defocusmax is", defocusmax
										
					
					stripdefocus = fitdefocus( fftc, angle, apix, options, len(micrographstarts), m, defocusmin, defocusmax, defocusstep, micrographstarts[m])
				
					if stripdefocus:
						print "defocus for strip at x %d is %.6f" %( micrographstarts[m] , stripdefocus )
					else:
						print "WARNING! bad strip; defocus for strip at x %d is None" %( micrographstarts[m] )

					
			
			
				#print "micrographmiddlein micrometers is", micrographmiddle
				print "stripdefocus", stripdefocus
				if stripdefocus:
					#imgdefocuses.append( stripdefocus*10000/apix )		#defocus in pixels
					imgdefocuses.append( stripdefocus )					#defocus in micrometers
				
					micromids.append( micrographmiddle )	
					print "\nappended (good) micrographmiddle is", micrographmiddle	
				else:
					print "\nappending to failed results"
					faileddefs.append( (defocusmin+defocusmax)/2 )
					failedmids.append( micrographmiddle )
			
			#xs = numpy.array( [i*options.stripstep + options.tilesize/2.0 for i in range(len(imgdefocuses))] )
		
			#print "micromids are", micromids
			xs =numpy.array( micromids )
		
			imgdefocuses = numpy.array( imgdefocuses )
		
			#print 'xs are', xs, type(xs)
			print "\n\n\n\nimgdefocuses", imgdefocuses
			#print "\bPLOTTING\n\n"
		
			m=0
			b=0
			defocuscalc = 0	
			nxMicrometers = nx*apix/10000.00
		else:
			print "\n\n\n\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!img being skipped for strip-based fitting"
			xs.append(micrographmiddle)
			imgdefocuses.append(globaldefocus)
			anglestoexclude.append(angle)
		
		xs = numpy.array( xs )
		imgdefocuses = numpy.array( imgdefocuses )
		
		print "\nfor label",label
		print "imgdefocuses are",imgdefocuses
		print "xs are",xs
		
		if xs.any() and imgdefocuses.any():
		
			for y in range( len(imgdefocuses)):
				#if math.fabs(imgdefocuses[y]-globaldefocus)
				imgdefocuses[y] = imgdefocuses[y] *-1
			


			
			#if m < 0:
			#	if angle > 0:
			#		m *= -1
			#if m > 0:
			#	if angle < 0:
			#		m *= -1
			
			
			
			
			if len(xs) < 2 and len (imgdefocuses) < 2:
				anglecalc = angle
				print "defocus calc is derived from middle strip only! it would have been", defocuscalc
				defocuscalc = imgdefocuses[0]
				print "but is", defocuscalc
				if defocuscalc < 0.0:
					defocuscalc*=-1.0
			else:
				m, b = numpy.polyfit(xs, imgdefocuses, 1)
				
				defocuscalc = m * nxMicrometers/2.0 + b
				print "\n\n\n\nSSSSSSSS\nbased on angle, m, b, xs, imgdefocuses",angle,m,b,xs,imgdefocuses
				print "defocuscalc is",defocuscalc
				#print "the mirrored defocus is", (-1*m) * (-1*nxMicrometers)/2.0 + b
				print "while globaldefocus is",globaldefocus


				#positive angle --> positive slope --> more defocused left side, less defocused right side
				anglecalc = math.degrees(numpy.arctan( m ))
				
			if angle not in anglestoexclude:
				angerror = math.fabs( angle - anglecalc )
				
				print "angle, anglecalc, and angerror are", angle, anglecalc, angerror
				angerrors.update( { angle:angerror } )
				
				if math.fabs(angerror) > math.fabs(options.maxangerror):
					print "\nusing globaldefocus for label %s due to large angular error" %(label)
					defocuscalc = globaldefocus
					middef = imgdefocuses[0]
					if defocuscalc < 0.0:
						defocuscalc*=-1.0
					
					if len(imgdefocuses) > 2:
						middef = imgdefocuses[len(imgdefocuses)/2]
						defocuscalc = (globaldefocus+middef)/2.0
					
					
										
			
			

			sptctfplotter( options, nxMicrometers, xs, imgdefocuses, maxangle, angle, angles.index( angle ), len(angles), imgindx, m, b, globaldefocus, globalmiddle, faileddefs, failedmids )		

			
		else:
			print "\nWarning: All defocuses failed for this submicrograph. Nothing to plot."
			defocuscalc = globaldefocus
			if defocuscalc < 0.0:
					defocuscalc*=-1
		
		#if xs.any() and imgdefocuses.any():
			#pass nx in micrometers
		
		print "\n\n\n\n\n\n\n\n@@@@@@@@@@@@@@@@@@@@@\nbefore setting params, defocuscalc %f and globaldefocus %f" %(defocuscalc,globaldefocus)

		params = {'ampcont':options.ampcont,'apix':apix,'bfactor':options.bfactor,'cs':options.cs,'defocus':math.fabs(defocuscalc),'voltage':options.voltage}
	
		ctf = EMAN2Ctf()
		#ctf.from_dict({'defocus':params['defocus'],'bfactor':params['bfactor'],'ampcont':params['ampcont'],'apix':params['apix'],'voltage':params['voltage'],'cs':params['cs']})	
		ctf.from_dict(params)
		
		ctfs.update( {angle:ctf} )
				
		imgnum += 1
		
	
	angles.sort()
	
	angerrors = collections.OrderedDict(sorted(angerrors.items()))
	if angerrors:
	
		avgangerror = sum( [  math.sqrt(angerrors[a]*angerrors[a]) for a in angerrors.keys() ] ) /len( angerrors )
		
		a=open(options.path + '/angular_error_avg.txt','w')
		a.writelines([str(avgangerror)+'\n'])
		a.close()
		
		lines = []
		finalangles = []
		finalangerrors = []
		
		
		
		for angle in angerrors:
			line = str(angle) + ' ' + str(angerrors[angle]) + '\n'
			lines.append( line )
			finalangles.append( angle )
			finalangerrors.append( angerrors[angle] )
		
		fa = open( options.path + '/tiltangle_vs_angular_error.txt',  'w' )
		fa.writelines( lines )
		fa.close()
	
		xlabel = 'Tilt angle (degrees)'
		ylabel = 'Angular error (degrees)'
	
		plotname = options.path + '/tiltangle_vs_angular_error.png'	
		title = 'Tiltangle vs angular error'
		
		generalplotter( options, finalangles, finalangerrors, xlabel, ylabel, plotname, title, False )

	
	if allglobaldefocuses:
		
		allglobaldefocusesvals = []
		finalangles = []
		lines=[]
	
		for angle in angles:
			if angle in allglobaldefocuses:
				line = str( angle ) + ' ' + str( allglobaldefocuses[angle] ) + '\n'
				lines.append(line)
				#print "Line to write is", line
				allglobaldefocusesvals.append(allglobaldefocuses[angle])
				finalangles.append(angle)
			else:
				print "Global defocus fit failed for image at angle", angle
		
		f = open( options.path + '/tiltangle_vs_globaldefocus.txt', 'w' )
		f.writelines(lines)
		f.close()
		
		f = open( options.path + '/eucentricity_variation.txt', 'w' )
		f.writelines([ str( max(allglobaldefocusesvals) - min(allglobaldefocusesvals) ) +'\n' ] )
		f.close()
		
		ylabel = 'Global defocus (micrometers)'
		xlabel = 'Tilt angle (degrees)'
		plotname = options.path + '/tiltangle_vs_globaldefocus.png'
		title = 'Tiltangle vs global defocus'
		
		generalplotter( options, finalangles, allglobaldefocusesvals, xlabel, ylabel, plotname, title )
	else:
		print "WARNING! All global defocuses estimations failed!"
	

	return ctfs
	



def sptctfplotter( options, nx, xdata, ydata, maxangle, angle, angleindx, nangles, imgindx, m=0, b=0, gdefocus=None, gmid=None, failedys=[], failedxs=[] ):
	import matplotlib
	matplotlib.use("TkAgg")
	import matplotlib.pyplot as plt
	import pylab
	
	
	#change defocus values to negative so that the slope of the plot makes sense; they are passed in negative already
	#for y in range( len(ydata)):
	#	ydata[y] = ydata[y] *-1
	
	if gdefocus:
		gdefocus *= -1
	else:
		print "WARNING! No global defocus gdefocus for angle", angle
	
	print "failed data in plotter is", failedys, failedxs
	
	#fig = plt.figure()

	proportionx = 6
	proportiony = 6
	
	
	'''	
	xdatarange = max(xdata)-min(xdata)
	ydatarange = max(ydata)-min(ydata)
	
	if xdatarange and ydatarange:
		
		print "xdatarange", xdatarange
		print "ydatarange", ydatarange
		
		if xdatarange > ydatarange:
			scalefactor = xdatarange/ydatarange
			proportionx *= scalefactor
		elif ydatarange > xdatarange:
			scalefactor = ydatarange/xdatarange
			proportiony *= scalefactor
	
		print "scalefactor",scalefactor
	

	proportiony = math.fabs( 10* sin( math.degrees( maxangle) ) )
	proportionx = math.fabs( 10* cos( math.degrees( maxangle) ) )
	
	print "proportionx", proportionx
	print "proportiony", proportiony
	'''
	
	fig = plt.figure( figsize=(10, 10) )

	ax = fig.add_subplot(111)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.tick_params(axis='both',reset=False,which='both',length=8,width=3)
	
	matplotlib.rc('xtick', labelsize=16) 
	matplotlib.rc('ytick', labelsize=16) 
	
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)
	
	#pylab.plot(azs, values[ele], color=RGB_tuples[kont], linewidth=2)
	
	pylab.rc("axes", linewidth=3.0)
	
	pylab.xlabel('X Axis', fontsize=16, fontweight='bold')
	pylab.ylabel('Y Axis', fontsize=16, fontweight='bold')

	#print "max y is", max(yaxis)
	
	print "ydata is", ydata
	
	'''
	maxy = max(ydata)
	miny = min(ydata)
	
	'''
	
	
	yavg = sum( ydata )/len(ydata)
	miny = yavg
	maxy = yavg
	
	if gdefocus:
		if float(gdefocus) > float( maxy ):
			maxy = gdefocus
		if float(gdefocus) < float( miny ):
			miny = gdefocus
		
	if failedys:
		minyf = min( failedys )
		if float(minyf) < float(miny):
			miny = minyf * -1
		maxyf = max( failedys )
		if float(maxyf) > float( maxy ):
			maxy = maxyf * -1
	
	
	pylab.ylim([ miny-1, maxy+1 ])

	#print "max x is", max(xaxis)

	pylab.xlim([ 0, nx + 0.5 ])

	ax.set_xlabel('X coordinate (micrometers)', fontsize=18, fontweight='bold')
	ax.set_ylabel('Defocus (micrometers)', fontsize=18, fontweight='bold')
	
	title ="Tilt image " + str( angleindx ).zfill( len(str( nangles ))) + ", angle=" + str(angle)
	pylab.title( title, fontweight='bold', fontsize=18 )
	
	#plt.plot(xdata, ydata, '.', markersize=10)
	plt.scatter(xdata, ydata, marker='.', s=200,alpha=0.9, color='k')
	
	if failedxs and failedys:
		print "plotting failed data"
		
		for y in range( len(failedys)):
			failedys[y] = failedys[y] *-1
		
		plt.scatter( numpy.array(failedxs), numpy.array(failedys), marker='.', s=100,alpha=0.9, color='r')

	if m and b:
		plt.plot(xdata, m*xdata + b, '-', linewidth=3, alpha=0.75)
	
	print "\nglobal defocus and mid inside plotter are", gdefocus, gmid

	if gdefocus and gmid:	
		print "\nplotting global fit at", gdefocus, gmid
		gdefs= numpy.array([gdefocus])
		gmids= numpy.array([gmid])
		print "\nplotting global fit at", gdefs, gmids
		plt.scatter(gmids,gdefs,alpha=0.6,zorder=1,s=400,marker='o',facecolors='none', edgecolors='g',linewidth=3)
		#plt.plot( gdefs, gmids )
		#plt.show()
		
		imodm = math.tan( math.radians( angle ) )
		print "\nIMOD slope is",imodm
		
		#y = mx + b ; therefore, b = y - mx
		b = gdefocus - imodm * gmid
		
		x1 = xdata[0]
		y1 = imodm * x1 + b
		
		x2 = xdata[-1]
		y2 = imodm * x2 + b
		
		imodxs = [x1,x2]
		imodys = [y1,y2]
		
		plt.plot( imodxs, imodys, color='g', linewidth=3)
		 
	
	signtag = 'p'
	if angle < 0.0:
		signtag ='m'
	
	plotname = options.path + '/fit_' + str( imgindx ).zfill( len(str( nangles ))) + '_' + signtag + str( int(math.fabs( round(angle) ) )).zfill(3) +'.png'

	plt.savefig( plotname )
	#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nSSSSSSSSSSSSSSSSSSSS\nSaved plot"
	#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"

	plt.clf()
	plt.close(fig)
	plt.close("all")

	return


def generalplotter( options, xaxis, yaxis, xlabel, ylabel, plotname, title, flipyaxis=True, fit=False, ptclnx=0 ):
	import matplotlib
	matplotlib.use("TkAgg")
	import matplotlib.pyplot as plt
	import pylab
	
	sizerangex = math.fabs( max(xaxis)-min(xaxis) )
	sizerangey = math.fabs( max(yaxis)-min(yaxis) )
	
	proportionfactor = 1.0
	
	sizeplotx=15.0
	sizeploty=15.0
	
	if sizerangex > sizerangey:
		proportionfactor = sizerangey/sizerangex
		sizeploty =  int( round( sizeploty*proportionfactor ) )
		
	elif sizerangey > sizerangex:
		proportionfactor = sizerangex/sizerangey
		sizeplotx = int( round( sizeplotx*proportionfactor ) )
	print "\nsizerangex=%f, sizerangey=%f, proportionfactor=%f therefore sizeplotx=%d, sizeploty=%d" %(sizerangex,sizerangey,proportionfactor,sizeplotx,sizeploty)
		
	fig = plt.figure(figsize=(30, 3))
	
	'''
	FORMAT AXES
	'''
	
	plt.axis('equal')

	#change defocus values to negative so that the slope of the plot makes sense
	if flipyaxis:
		for y in range( len(yaxis)):
			yaxis[y] = yaxis[y] *-1
				
	matplotlib.rc('xtick', labelsize=16) 
	matplotlib.rc('ytick', labelsize=16) 
	
	font = {'weight':'bold','size':16}
	matplotlib.rc('font', **font)
	
	fig2 = plt.figure()
	
	#fig = plt.figure(figsize=(15, 6))
		
	ax = fig2.add_subplot(111)
	
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.tick_params(axis='both',reset=False,which='both',length=8,width=3)
	
	maxy = max(yaxis)
	miny = min(yaxis)
	yrange = maxy - miny
	extray = yrange/20
	
	#if options.yrange:
	#	ylim1 = int( options.yrange.split(',')[0] )
	#	ylim2 = int( options.yrange.split(',')[1] )
	#	yrange = options.yrange
	#	extray = 1
		

	maxx = max(xaxis)
	minx = min(xaxis)
	xrange = maxx - minx
	extrax = xrange/20
	
	
	if 'pixels' in xlabel or 'pixels' in ylabel:
		if options.radius:
			extray = extrax = options.radius * 2
		elif ptclnx:
			extray = extrax = int(ptclnx/2)
		
			
	#print "\nmax y is", maxy
	ylim1 = miny - extray
	ylim2 = maxy + extray
	
	xlim1 = minx - extrax
	xlim2 = maxx + extrax
	
	#if options.xrange:
	#	xlim1 = int( options.xrange.split(',')[0] )
	#	xlim2 = int( options.xrange.split(',')[1] )
	#	xrange = options.yrange
	#	extrax=1
		

	pylab.ylim([ylim1, ylim2])
	
	#print 'yrange', ylim1,ylim2
	
	print "\nmax x is", max(xaxis)
	
	pylab.xlim([xlim1, xlim2])
	
	#print 'xrange', xlim1,xlim2
	
	
	
	
	
	
	
	
	
	
	
	
	ax.set_xlabel(xlabel, fontsize=18, fontweight='bold')
	ax.set_ylabel(ylabel, fontsize=18, fontweight='bold')
	
	#title ="Tilt image " + str( angleindx ).zfill( len(str( nangles ))) + ", angle=" + str(angle)
	#pylab.title( title )
	
	pylab.title( title )
	
	plt.scatter(xaxis,yaxis,alpha=0.70,zorder=1,s=100,facecolors='b', edgecolors='b')
	
	if fit:
		m, b = numpy.polyfit(xaxis, yaxis, 1)
	
		if m and b:
			xarray = numpy.array( xaxis )
			plt.plot(xaxis, m*xarray + b, '-', linewidth=3, alpha=0.75,color='k',linestyle='--')
			
	plt.savefig( plotname )

	#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nSSSSSSSSSSSSSSSSSSSS\nSaved plot"
	#print "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
	
	plt.clf()
	plt.close(fig)
	plt.close("all")

	
	return	

	

def bgAdj(ctf,fg_1d):
	"""Smooths the background based on the values of the foreground near the CTF zeroes and puts the
	smoothed background into the CTF object"""
	ds=ctf.dsbg
	#print "\n (bgAdj) ds is", ds
	ctf=ctf
	#print "\n (bgAdj) ctf is", ctf
	bg_1d=list(ctf.background)
	#print "\n (bgAdj) bg_1d is", bg_1d

	xyd=XYData()

	# Find the minimum value near the origin, which we'll use as a zero (though it likely should not be)
	mv=(fg_1d[1],1)
	fz=int(ctf.zero(0)/(ds*2))
	
	

	
	for lz in xrange(1,fz):
		mv=min(mv,(fg_1d[lz],lz))

	xyd.insort(mv[1],mv[0])

	# now we add all of the zero locations to our XYData object
	for i in xrange(100):
		z=int(ctf.zero(i)/ds)
		if z>=len(bg_1d)-1: break
		if fg_1d[z-1]<fg_1d[z] and fg_1d[z-1]<fg_1d[z+1]: mv=(z-1,fg_1d[z-1])
		elif fg_1d[z]<fg_1d[z+1] : mv=(z,fg_1d[z])
		else : mv=(z+1,fg_1d[z+1])
		xyd.insort(mv[0],mv[1])

	# new background is interpolated XYData
	ctf.background=[xyd.get_yatx_smooth(i,1) for i in xrange(len(bg_1d))]

	# if our first point (between the origin and the first 0) is too high, we readjust it once
	bs=[fg_1d[i]-ctf.background[i] for i in xrange(fz)]
	
	#print "bs first is", bs
	
	if min(bs)<0 :
		mv=(bs[0],fg_1d[0],0)
		for i in xrange(1,fz): mv=min(mv,(bs[i],fg_1d[i],i))
		xyd.set_x(0,mv[2])
		xyd.set_y(0,mv[1])
		
		ctf.background=[xyd.get_yatx_smooth(i,1) for i in xrange(len(bg_1d))]
		
		#bs2=[fg_1d[i]-ctf.background[i] for i in xrange(fz)]
		#print "bs second is", bs2
		
			 

if '__main__' == __name__:
	main()
