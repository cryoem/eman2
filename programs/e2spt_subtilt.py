#!/usr/bin/env python

# Author: Jesus Galaz-Montoya, 02/Feb/2013, last update 12/Feb/2015
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


import os
from EMAN2 import *
import sys
import numpy
import math


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """UNDER DEVELOPOMENT. Extracts particles from each image in an aligned tilt 
		series based on A) their position in the reconstructed tomogram
		or B) their position in the 0 degrees tilt image of a tilt series."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	'''
	Parameters for adding ctf and noise
	'''
	parser.add_argument('--tiltseries',type=str,default='',help="""File in .ali, .mrc or .hdf format of the aligned tiltseries.""")
		
	parser.add_argument('--tiltangles',type=str,default='',help="""File in .tlt or .txt format containing the tilt angle of each tilt image in the tiltseries.""")
		
	parser.add_argument('--coords3d',type=str,default='',help="""File in .txt format containing the coordinates of particles determined from the reconstructed tomogram of the supplied tiltseries.""")
	
	parser.add_argument('--coords2d',type=str,default='',help="""File in .txt format containing the coordinates of particles determined from the aligned 0 tilt image in the supplied tiltseries.""")
	
	parser.add_argument('--zerotiltindx',type=int,default=-1,help="""The default is the image at the middle of the stack. Since the stack might have more images to the left or the right of the actual 0-tilt (or lowest tilt) image, you can explicitly provide the index of the lowest tilt image here. This is used for tracking images.""")
	
	parser.add_argument('--centerzerotilt',action='store_true',default=False,help="""Default=False. If specified, this option will center the zerotilt (or least tilted image) for each particle by using as a reference a sharp-circle of radius=box/2 or the value specified through --radius.""")
	
	parser.add_argument('--excludeedge',type=float,default=0.0,help="""Integer number of pixels away from the edge of each image in the tilt series to not extract particles from. For example, if you specify 100, and the images are 4096x4096 pixels, any particle with its center lying between 0 and 200 or 3896 and 4096 will nto be extracted.""") 
	
	parser.add_argument('--cshrink', type=int, default=1, help="""Specifies the factor by which to multiply the coordinates in the coordinates file, so that they can be at the same scale as the tomogram. For example, provide 2 if the coordinates are on a 2K x 2K scale, but you want to extract the particles' subtiltseries from the UN-shrunk 4K x 4Ktiltseries.""")
	
	parser.add_argument('--tomogram',type=str,default='',help="""Path to raw, unbinned tomogram.""")
	
	parser.add_argument('--saveanglestacks',type=str,default='',help="""Default=None. Comma separated values of tilt angle indexes for which you want to save all particles as a stack. For example, if you want all particles from the 0 tilt image, you would provide the index for that image in the tilt series. In a tilt series with 61 images (1-61), the 0 tilt image is probably image number 31, so you would say --saveanglestakcs=31, and all the particles from the 0 tilt image would be put into a single HDF stack.""")

	parser.add_argument('--tiltaxislocation',type=int,default=-1,help="""By default, the tilt axis will be assumed to run through the middle of the tomogram in X, parallel to the Y axis. For example, if the dimensions of the tomogram are 4096x3000x500, the tilt axis will be assumed to be at X=2048. Provide a different integer number to change the location of the tilt axis (it will still be assumed to be parallel to Y though).""") 
	
	parser.add_argument('--tiltaxisptcls',type=int,default=-1,help="""Specifies the distance from the tilt axis to consider particles for extraction. By default, all particles will be extracted. However, if you provide, for example, --tiltaxisptls=10, only particles with centers -10 to 10 pixels away from the tilt axis will be extracted.""")
	
	parser.add_argument('--ntiltslow',type=int,default=0,help="""Default=0 (not used). If you supply an even number 1 will be added to it (for example, 4 will be turned into 5). If --ntiltslow>0, it specifies the number of tiltimages to keep in each subtiltseries, starting from the zero-tilt image and incorporating particles from right and left, one at a time. For example, in a tiltseries from -60 to 60 degress with a step size of 2 degrees, --ntiltslow=5 would keep tiltimages at angles 0,2,-2,-4,-4.""")
	
	parser.add_argument('--ntiltslowneg',type=int,default=0,help="""Default=0 (not used). If --ntiltslowneg>0, it specifies the number of tiltimages to keep in each subtiltseries, starting from the zero-tilt image and progressively incorporating particles from negatively tilted images only. For example, in a tiltseries from -60 to 60 degress with a step size of 2 degrees, --ntiltslowneg=5 would keep tiltimages at angles 0,-2,-4,-6,-8.""")

	parser.add_argument('--ntiltslowpos',type=int,default=0,help="""Default=0 (not used). If --ntiltslowpos>0, it specifies the number of tiltimages to keep in each subtiltseries, starting from the zero-tilt image and progressively incorporating particles from positively tilted images only. For example, in a tiltseries from -60 to 60 degress with a step size of 2 degrees, --ntiltslowpos=5 would keep tiltimages at angles 0,+2,+4,+6,+8.""")
	
	#parser.add_argument('--ntiltshighneg',type=int,default=0,help="""Default=0 (not used). If --ntiltshighneg>0, it specifies the number of tiltimages to keep in each subtiltseries, starting from the zero-tilt image and progressively incorporating particles from negatively tilted images only. For example, in a tiltseries from -60 to 60 degress with a step size of 2 degrees, --ntiltshighneg=5 would keep tiltimages at angles -60,-58,-56,-54,-52.""")

	#parser.add_argument('--ntiltshighpos',type=int,default=0,help="""Default=0 (not used). If --ntiltshighpos>0, it specifies the number of tiltimages to keep in each subtiltseries, starting from the zero-tilt image and progressively incorporating particles from positvely tilted images only. For example, in a tiltseries from -60 to 60 degress with a step size of 2 degrees, --ntiltshighpos=5 would keep tiltimages at angles  +60,+58,+56,+54,+52.""")
	
	parser.add_argument('--tomosides',type=str,default='',help="""Comma separated values for the tomogram dimensions. Alternatively, provide the path to the tomogram itself through --tomogram.""")
		
	parser.add_argument("--icethicknessauto",action='store_true',default=False,help="""Default=False. If supplied, the thickness of the tomogram in Z will be calculated by computing the difference between the largest and the smallest Z coordinate found in the --coords3d coordinates file.""")
	
	parser.add_argument("--zshift",type=str,default='half',help="""By default, the tomogram will be shifted -half the ice thickness so that the middle of the tomogram is at z=0. Provide a positive or negative integer to shift the z position by a different amount""")
	
	parser.add_argument("--radius",type=int,default=0,help="""Default=0 (not used). Radius of the particle in pixels. 2*radius will be added to the icethickness if --radius AND --icethicknessauto are supplied.""") 
	
	parser.add_argument("--invertangles",action='store_true',default=False,help="""Default=False. If True, this will multiple all angles by -1, in case the directionality is messed up.""")
	
	parser.add_argument('--path',type=str,default='spt_subtilt',help="""Directory to save the results.""")
	
	parser.add_argument('--boxsize',type=int,default=128,help="""Size of the 2D "tiles" or images for each particle from each image in the tiltseries.""")
	
	parser.add_argument('--apix',type=float,default=0.0,help="""If provided, this value will be used for apix instead of the one read from the header of --tiltseries""")
		
	parser.add_argument("--ppid", type=int, help="""Set the PID of the parent process, used for cross platform PPID""",default=-1)
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	parser.add_argument('--subset', type=int, default=0, help='''Specify how many sub-tiltseries (or particles) from the coordinates file you want to extract; e.g, if you specify 10, the first 10 particles will be boxed.\n0 means "box them all" because it makes no sense to box none''')
	
	#parser.add_argument('--tomogramthickness',type=int,default=None,help='Z dimension of the reconstructed tomogram.')
		
	#parser.add_argument("--everyother", type=int, help="""Pick every other tilt. For example, --tilt=3 would pick every third tilt only.""",default=-1)

	parser.add_argument("--subtractbackground",action='store_true',default=False,help="""(Experimental. Not working yet). This will extract a box from the tomogram much larger than the subtomogram. Projections will be generated. You MUST provide --tomogram for this.""")

	parser.add_argument("--normproc",type=str,default='',help="""WARNING: Not used anywhere yet. Default=None""")
	
	#parser.add_argument("--yshort",action='store_true',help="""Not used anywhere yet. Default=False""", default=False)
	
	parser.add_argument("--shrink", type=int,default=1,help="""Default=1 (no shrinking). Integer shrinking factor, part of preprocessing to facilitate particle tracking.""")
	
	parser.add_argument("--lowpass", type=str,default='',help="""Default=None. Requires --track. Low pass filtering processor (see e2help.py processors at the command line), part of preprocessing to facilitate particle tracking.""")

	parser.add_argument("--highpass", type=str,default='',help="""Default=None (no highpass). High pass filtering processor (see e2help.py processors at the command line), part of preprocessing to facilitate particle tracking.""")
	
	parser.add_argument("--mask", type=str,default='',help="""Default=None. Requires --track. Masking processor (see e2help.py processors at the command line), part of preprocessing to facilitate particle tracking.""")

	parser.add_argument("--threshold", type=str,default='',help="""Default=None (no threshold). Requires --track. Thresholding processor (see e2help.py processors at the command line), part of preprocessing to facilitate particle tracking.""")

	parser.add_argument("--preprocess", type=str,default='',help="""Default=None (no additional preprocessing). Requires --track. Any additional preprocessing processor (see e2help.py processors at the command line), part of preprocessing to facilitate particle tracking.""")
	
	parser.add_argument("--track",action='store_true',default=False,help="""Default=False (not used). If supplied, this option will track particles from one tilt image to another.""")
	
	(options, args) = parser.parse_args()
	
	if options.ntiltslow:
		if options.ntiltslowneg:
			print "ERROR: Cannot specify --ntiltslow and --ntiltslowneg at the same time"
			sys.exit()
		if options.ntiltslowpos:
			print "ERROR: Cannot specify --ntiltslow and --ntiltslowpos at the same time"
			sys.exit()
		
	if options.ntiltslowneg:
		if options.ntiltslow:
			print "ERROR: Cannot specify --ntiltslowneg and --ntiltslow at the same time"
			sys.exit()
		if options.ntiltslowpos:
			print "ERROR: Cannot specify --ntiltslowneg and --ntiltslowpos at the same time"
			sys.exit()
	
	if options.ntiltslowpos:
		if options.ntiltslow:
			print "ERROR: Cannot specify --ntiltslowpos and --ntiltslow at the same time"
			sys.exit()
		if options.ntiltslowneg:
			print "ERROR: Cannot specify --ntiltslowpos and --ntiltslowneg at the same time"
			sys.exit()
		
		
	print "\nI've read the options"	
	
	'''
	Check that all needed parameters are properly supplied
	'''
	if not options.tiltseries:  
		print "ERROR: You must provide --tiltseries."
		sys.exit()
	
	if not options.tiltangles:
		print "ERROR: You must provide --tiltangles."
		sys.exit()
	
	if not options.coords2d and not options.coords3d:
		print "ERROR: You must provide EITHER --coords2d OR --coords3d." 
		sys.exit()
	
	if options.coords2d and options.coords3d:
		print "ERROR: You must provide EITHER --coords2d OR --coords3d, not both." 
		sys.exit()
	
	if not options.tomogram and not options.tomosides:
		print "ERROR: You must provide EITHER --tomogram OR --tomosides." 
		sys.exit()
		
	if options.tomogram and options.tomosides:
		print "ERROR: You must provide EITHER --tomogram OR --tomosides, not both." 
		sys.exit()
	
	'''
	Parse tilt angles from file supplied via --tiltangles
	'''
	anglesfile = open(options.tiltangles,'r')				#Open tilt angles file
	alines = anglesfile.readlines()							#Read its lines
	anglesfile.close()										#Close the file
	
	#tiltangles = [ alines[i].replace('\n','') for i in range(len(alines)) ]	#Eliminate trailing return character, '\n', for each line in the tiltangles file
	
	tiltanglesfloatabs = []
	for line in alines:
		print "line is", line
		ang = math.fabs( float( line.replace('\n','') ) )
		print "ang is", ang
		tiltanglesfloatabs.append( ang )
	#tiltanglesfloatabs = [ math.fabs( float( alines[i].replace('\n','') ) ) for i in range(len(alines)) ]

	tiltanglesfloat = [ float( alines[i].replace('\n','') ) for i in range(len(alines)) ]
	
	#if options.invertangles:
	#	tiltanglesfloat = [ -1*float( alines[i].replace('\n','') ) for i in range(len(alines)) ]
	#	print "INVERTED tiltanglesfloat", tiltanglesfloat
	#else:
	print "tiltanglesfloat", tiltanglesfloat

		
	
	#tiltanglesfloat.sort()
	
	#print "sorted tiltangles",tiltanglesfloat
	
	ntiltangles = len( tiltanglesfloat )
	
	'''
	Check that there are as many images in --tiltseries as angles in --tiltangles
	'''
	serieshdr = EMData(options.tiltseries,0,True)
	nslices = serieshdr['nz']
	
	if int( nslices ) != int( ntiltangles ):
		print """ERROR: The tiltangles file doesn't seem to correspond to the tiltseries provided.
				The number of images in --tiltseries (z dimension of MRC stack) must be equal to the number
				of lines in --tiltangles."""
		sys.exit()
	
	'''
	Get apix
	'''
	nx = serieshdr['nx']
	ny = serieshdr['ny']
	apix = serieshdr['apix_x']
	if options.apix:
		apix=options.apix
	
	'''
	If error free up to this point, make path to store results
	'''
	from e2spt_classaverage import sptmakepath
	options = sptmakepath(options,'sptSubtilt')
	
	'''
	Start logging this run of the program at this point
	'''
	logger = E2init(sys.argv, options.ppid)
	
	
	#"""
	#(CRAZY)
	#You do not need to keep track of the mathematics of tilting and rotating and find the correspondence between
	#the tomogram and each image in the tilt series.
	#Instead, let's make a simple 3D model, containing a bright dot at the position of each particle.
	#Then, rotate that using the known tilt angles, generate a projection for each, and find the dots (maxima) in each.
	#Use the location of these maxima, to extract each particle from each image in the tilt series.
	#Would fail if they overlap though
	#"""
	
	tomox=tomoy=tomoz=0
	if options.tomosides:							#Read tomogram dimensions.
		sides=options.tomosides.split(',')
		tomox = int(sides[0])
		tomoy = int(sides[1])
		tomoz = int(sides[2])
	
	elif options.tomogram:
		tomohdr = EMData(options.tomogram,0,True)	#Read tomogram dimensions from tomogram header, if the tomogram is provided.
		tomox = int(tomohdr['nx'])
		tomoy = int(tomohdr['ny'])
		tomoz = int(tomohdr['nz'])

	icethickness = tomoz
	
	#if float( options.shrink ) > 1.0:								#The 'MODEL' to build for the coordinates need not use the full size of the tomogram.
	#	tomox = int(tomox)/options.shrink
	#	tomoy = int(tomoy)/options.shrink
	#	tomoz = int(tomoz)/options.shrink
	
	#tomovol = EMData(tomox,tomoy,tomoz)			#Create empty volume for the MODEL to build
	#tomovol.to_zero()								#Make sure it's empty
	
	clines=[]
	if options.coords2d:
		cfile = open(options.coords2d,'r')				
		clines = cfile.readlines()						
		cfile.close()									
	
	elif options.coords3d:
		cfile = open(options.coords3d,'r')				
		clines = cfile.readlines()						
		cfile.close()									
		
	
	'''
	Clean the coordinate file lines (clines) if there's garbage in them.
	Some people might manually make ABERRANT coordinates files with commas, tabs, or more 
	than one space in between coordinates. Then, parse each clean line.
	'''
	ppp = 0	
	cleanlines=[]
	zs = []
	
	for line in clines:
		
		if options.subset:
			if int( ppp ) >= (options.subset):
				break
			
		line =line.replace(", ",' ')	
		line = line.replace(",",' ')
		line = line.replace("x",'')
		line = line.replace("y",'')
		line = line.replace("z",'')
		line = line.replace("=",'')
		line = line.replace("_",' ')
		line = line.replace("\n",'')
		line = line.replace("\t",' ')
		line = line.replace("  ",' ')
		
		finallineelements=line.split(' ')
		
		if options.coords3d:
			if line and len(finallineelements) == 3:
				cleanlines.append(line)
				zs.append( int( line.split()[-1] ))		
				ppp += 1
			else:
				print "\nBad line removed", line
	
		elif options.coords2d:
			if line and len(finallineelements) == 2:
				cleanlines.append(line)
				ppp += 1
			else:
				print "\nBad line removed", line
	
	print "icethicknessauto before if is", options.icethicknessauto
	if zs and options.icethicknessauto and not options.tomogram and not options.tomosides:
		print "icethicknessauto after if is", options.icethicknessauto
	
		icethicknessfile = options.path + '/icethickness_estimate.txt'
	
		itf = open( icethicknessfile, 'w' )
	
		autoIcethickness = max( zs ) -  min( zs )
	
		autoIceLine = 'max(z) - min(z) = ' + str( autoIcethickness ) + ' pixels '
	
		autoIcethicknessUnshrunk = autoIcethickness
		if options.cshrink and options.cshrink > 1:
			autoIcethicknessUnshrunk = autoIcethickness * options.cshrink
			autoIceLine += ' , unshrunk = ' + str( autoIcethicknessUnshrunk ) + ' pixels'
	
		icethickness = autoIcethicknessUnshrunk
		print "\nIce thickness has been estimated from z coordinates to be", icethickness
	
		autoIcethicknessInApix = autoIcethicknessUnshrunk*apix
		autoIceLine += ' , ' + str( autoIcethicknessInApix ) + ' angstroms'
		
		itf.write( autoIceLine )
	
		itf.close()		
	
	'''
	Iterate over the correct number of viable lines from the coordinates file.
	'''
	
	nptcls = len(cleanlines)
	if int( options.subset ) > 0:
		if int( options.subset ) > len(cleanlines):
			print """WARNING: The total amount of lines in the coordinates files is LESS 
				than the --subset of particles to box you specified; therefore, ALL particles 
				will be extracted."""
		else:
			nptcls = int(options.subset)
			print "\nThe SUBSET of particles to work with is", nptcls
	else:
		print """\nBased on the number of coordinates, the size of the ENTIRE SET of 
			subtiltseries to extract is""", nptcls
	
	#print "There are these many clean lines", len(cleanlines)
	#print "Clean lines are", cleanlines
	
	#everyotherfactor = 1
	#if options.everyother > 1:
	#	everyotherfactor = options.everyother
	
	
	maxtilt=0
	
	
	if options.subtractbackground:
		maxtilt = max( tiltanglesfloat )
		print "\n(e2spt_subtilt.py) maxtilt is", maxtilt
		
		from e2spt_boxer import unbinned_extractor
		
		bgboxsize = (2 * options.boxsize / math.cos( math.radians(maxtilt+5)  )) + 10
		invert=0
		center=0
	
	
	
	
	tiltaxisloc = tomox/2.0 
	if options.tiltaxislocation > -1:
		tiltaxisloc = options.tiltaxislocation/2.0 
	
	xclowerthresh = 0
	xcupperthresh = tomox
	
	if options.radius and options.tiltaxisptcls == -1:
		xclowerthresh = options.radius
		xcupperthresh = tomox - options.radius
	
	if options.tiltaxisptcls > -1:
		xclowerthresh = tiltaxisloc - options.tiltaxisptcls
		xcupperthresh = tiltaxisloc + options.tiltaxisptcls
	
	
	'''
	Divide the tiltangles into negative and positive ranges since iterations will progress
	from the 0 tilt angle and then gradually 1 to the left, 1 to the right, etc, etc.
	Figure out which side has the most angles too.
	'''
	nimgsOriginal = EMData( options.tiltseries, 0, True )['nz']
	middleIndx = nimgsOriginal/2
	
	absangles = []	#absolute value of all angles
	
	for ang in tiltanglesfloat:
		
		absangles.append( math.fabs( ang ) )
	
	zerotiltangleabs = min( absangles )
	zerotiltindx = absangles.index( zerotiltangleabs )
	zerotiltangle = tiltanglesfloat[ zerotiltindx ]
	
	print "\nBy middle tilt series index, and by smallest tilt angle index, zerotilt is at", middleIndx, zerotiltindx
	print "And zerotiltangle is", zerotiltangle
	
	
	anglesBelowZeroTilt = []
	anglesAboveZeroTilt = []

	for ang in tiltanglesfloat:
		#print "\nAngle is", ang
		if float(ang) < float(zerotiltangle):
			#if not options.invertangles:
			anglesBelowZeroTilt.append( float(ang) )
			#else:
			#	anglesAboveZeroTilt.append( -1*float(ang) )
				
			#print "Therefore it was added to the below group"
		elif float(ang) > float(zerotiltangle):
			#if not options.invertangles:
			anglesAboveZeroTilt.append( float(ang) )
			#else:
			#	anglesBelowZeroTilt.append( -1*float(ang) )
			#print "Therefore it was added to the above group"
	
	nLangles = len( anglesBelowZeroTilt )
	nUangles = len( anglesAboveZeroTilt )
	mostangles = max( nLangles, nUangles )
	
	#if not options.invertangles:
	print "anglesBelowZeroTilt", anglesBelowZeroTilt
	print "anglesAboveZeroTilt", anglesAboveZeroTilt
	#else:
	#print "INVERTED anglesBelowZeroTilt", anglesBelowZeroTilt
	#print "INVERTED anglesAboveZeroTilt", anglesAboveZeroTilt
	#print "nLangles and nUangles are", nLangles, nUangles
	#sys.exit()
	anglestackslist = options.saveanglestacks.split(',')


	ptclNum3D = ptclNum = 0
	
	if options.ntiltslow:				#To pick a symmetric number of tilt images, left and right of the middle tilt?
		if not options.ntiltslow %2:
			options.ntiltslow += 1
	
	for line in cleanlines:
		line = line.split()	
		
		print "\n\n\n\n\n+=================+\nAnalyzing particle number\n+=================+\n", ptclNum3D
		xc = 0
		yc = 0
		zc = 0
		if len(line) == 3:
			xc = float(line[0])				#Determine x y z coordinates for each line
			yc = float(line[1])
			zc = float(line[2])	
		elif len(line) == 2:
			xc = float(line[0])				#Determine x y coordinates for each line, and default zc to 0 if --coords2d supplied instead of coords3d
			yc = float(line[1])
			zc = 0
		else:	
			print "\nThere's still an aberrant line in your coordinates file, see", line
			sys.exit()
		
		if options.verbose:
			print "\nRead these coordinates", xc,yc,zc
	
		if options.cshrink:
			xc*=options.cshrink
			yc*=options.cshrink
			zc*=options.cshrink
			
			if options.verbose:
				print "\nThe real coordinates after multiplying cshrink are", xc, yc, zc
		
		print "Before entering algorithm, xc, yc, zc are", xc,yc,zc
		print "Both conditions xc > lowerth and xc < upperth together are", int(xc) > int( xclowerthresh ) and int(xc) < int( xcupperthresh )
		
		if int(xc) > int( xclowerthresh ) and int(xc) < int( xcupperthresh ):
		
			
		
			#outIndx=0
			#ret=0
			wholebox=0
			
			
			'''
			The middle or zero-tilt or lowest-tilt image will serve to align all the rest if --track is on.
			If --centerzerotilt is on, find it, extract it, recenter it using autocentering based on mirror images, and reextract.
			'''
			
			retm = extract2D( options, zerotiltangle, icethickness, tomox, xc, yc, zc, 0, 0, zerotiltindx )
			middleslice = retm[0]
			
			if options.coords3d:
				sptcoords = ( xc, yc, zc )
				middleslice['ptcl_source_coord']=sptcoords

			cumulativeLdx = cumulativeLdy = cumulativeUdx = cumulativeUdy = 0	
			
			ptclfile = options.path + '/subtiltPtcl_' + str(ptclNum).zfill( len( str( len (cleanlines)))) + '.hdf'
				
			print "Found least tilted image at index and angle ", zerotiltindx, zerotiltangle
			
			'''
			print "Autocentering it"
			
			midslicemirrorX = middleslice.process('xform.mirror',{'axis':'x'})
			midslicemirrorY = middleslice.process('xform.mirror',{'axis':'y'})
			
			
			
			ccfpmx = middleslice.calc_ccf( midslicemirrorX )
			ccfpmy = middleslice.calc_ccf( midslicemirrorY )

			ccfpmxC = ccfpmx.process('xform.phaseorigin.tocorner')
			ccfpmyC = ccfpmy.process('xform.phaseorigin.tocorner')

			maxccfpmxC = ccfpmxC.calc_max_location()
			maxccfpmyC = ccfpmyC.calc_max_location()
			
			
			
			midsxt = -1 * ( middleslice['nx'] /2.0 - maxccfpmxC[0])/ 2.0
			midsyt = -1 * ( middleslice['ny'] /2.0 - maxccfpmyC[1])/ 2.0
			
			print "Autocentering translations for tilt 0 are", midsxt, midsyt
			retm2 = extract2D( options, zerotiltangle, icethickness, tomox, xc, yc, zc, midsxt, midsyt, zerotiltindx )
			
			middleslice = retm2[0]
			'''
			
			midscoordx = xc
			midscoordy = yc
			
			if options.centerzerotilt:					
				box = middleslice['nx']
				radius = box / 4.0
				
				if options.radius:
					radius = options.radius
				else:
					print """WARNING: --centerzerotilt requires --radius. Since the latter wasn't provided, it will be assumed to be 1/4 of the --boxsize"""
					
				template = EMData( box, box )
				template.to_one()
				template.process_inplace('mask.sharp',{'outer_radius':radius})
				template.mult(-1)
			
				if options.lowpass or options.highpass or options.mask or options.shrink or options.preprocess or options.threshold:
					middleslice = preprocImg( middleslice, options )
			
				ccft = middleslice.calc_ccf( template )
				ccftC = ccft.process('xform.phaseorigin.tocorner')
				maxccftC = ccftC.calc_max_location()
			
				midsxt = -1 * ( middleslice['nx'] /2.0 - maxccftC[0])			
				midsyt = -1 * ( middleslice['nx'] /2.0 - maxccftC[1])
				print "Autocentering translations for tilt 0 are", midsxt, midsyt
				#cumulativeLdx += midsxt
				#cumulativeLdy += midsyt
			
				retm2 = extract2D( options, zerotiltangle, icethickness, tomox, xc, yc, zc, midsxt, midsyt, zerotiltindx )
			
				middleslice = retm2[0]
				midscoordx = retm2[1]
				midscoordy = retm2[2]
				
				if options.coords3d:
					sptcoords = ( midscoordx, midscoordy, zc )
					middleslice['ptcl_source_coord']=sptcoords
			
			
			#middleslice.process_inplace('normalize')
			
			middleslice['spt_tiltangle'] = zerotiltangle
			middleslice['spt_tiltaxis'] = 'y'
			middleslice['spt_subtilt_x'] = midscoordx
			middleslice['spt_subtilt_y'] = midscoordy
			middleslice['origin_x'] = middleslice['nx']/2.0
			middleslice['origin_y'] = middleslice['ny']/2.0
			middleslice['origin_z'] = 0

			middleslice['apix_x'] = apix
			middleslice['apix_y'] = apix
			middleslice['apix_z'] = apix
			
			
				
			middleslice.write_image( ptclfile, 0 )
			
			if str(zerotiltindx) in anglestackslist and options.saveanglestacks:	
				zeroanglestackname = options.path + '/anglestack_' + str(zerotiltindx) + '_angle' + str(int(zerotiltangle)) + '.hdf' 
				middleslice.write_image( zeroanglestackname, -1 )
			
			
			print "Iteration over all tilt angles (other than the least tilted one) starting now."
			
			refL = middleslice.copy()
			refU = middleslice.copy()
			
			print "mostangles is", mostangles
			print "nLangles is", nLangles
			print "nUangles is", nUangles
			print "tiltanglesfloat are", tiltanglesfloat
			
			
			leftanglesN = len( tiltanglesfloat[ 0: zerotiltindx ] )
			
			rightanglesN = len( tiltanglesfloat[ zerotiltindx: -1 ] )
			
			for k in range( mostangles ):
				
				if options.ntiltslow:
					if k == options.ntiltslow:
						break
				
				elif options.ntiltslowneg:
					if k == options.ntiltslowneg:
						break
				
				elif options.ntiltslowpos:
					if k == options.ntiltslowpos:
						break
				
				if k < leftanglesN and not options.ntiltslowpos: # and zerotiltindx - (k+1) > -1:
					#print "\n k and nLangles are", k, nLangles
					lowerindx = zerotiltindx - (k+1)
					
					#if not options.negativetiltseries:
					#	lowerindx = zerotiltindx + (k+1)
						
					
					#if options.ntiltslow:
					#	lowerindx = zerotiltindx - (k+1) - zerotiltiltidx
					
					lowerangle = tiltanglesfloat[ lowerindx ]
					
					#if options.invertangles:
					#	lowerangle *= -1
					#	print "stacking INVERTED lowerangle", lowerangle
					
					print "stacking angle %f from the LEFT of lowest tiltangle in the angles list" %( lowerangle )
					print "found at lowerindx", lowerindx
					
					retL = write2D( options, lowerangle, icethickness, tomox, tomoy, xc, yc, zc, cumulativeLdx, cumulativeLdy, refL, apix, 'lower', ptclfile, maxtilt, lowerindx )
					
					if retL:
						refL = retL[0]
						cumulativeLdx = retL[1]
						cumulativeLdy = retL[2]
				
						if str(lowerindx) in anglestackslist and options.saveanglestacks:	
							anglestackname = options.path + '/anglestack_' + str(lowerindx) + '_angle' + str(int(lowerangle)) + '.hdf' 
							
							
							
							
							#middleslice['spt_tiltangle'] = zerotiltangle
							#middleslice['spt_tiltaxis'] = 'y'
							#middleslice['spt_subtilt_x'] = midscoordx
							#middleslice['spt_subtilt_y'] = midscoordy
							#middleslice['origin_x'] = middleslice['nx']/2.0
							#middleslice['origin_y'] = middleslice['ny']/2.0
							#middleslice['origin_z'] = 0

							#middleslice['apix_x'] = apix
							#middleslice['apix_y'] = apix
							#middleslice['apix_z'] = apix
							
							
							
							refL.write_image( anglestackname, -1 )
				
				if k < rightanglesN and not options.ntiltslowneg:
					#print "\n k and nUangles are", k, nUangles
					upperindx = zerotiltindx + (k+1)
					#if not options.negativetiltseries:
					#	upperindx = zerotiltindx - (k+1)
					
					upperangle = tiltanglesfloat[ upperindx ]
					
					#if options.invertangles:
					#	upperangle *= -1
					#	print "stacking INVERTED upperangle", upperangle
						
					print "stacking angle %f from the RIGHT lowest tiltangle in the angles list" %( upperangle )
					print "found at upperindx", upperindx
					retU = write2D( options, upperangle, icethickness, tomox, tomoy, xc, yc, zc, cumulativeUdx, cumulativeUdy, refU, apix, 'upper', ptclfile, maxtilt, upperindx )
					if retU:
						refU = retU[0]
						cumulativeUdx = retU[1]
						cumulativeUdy = retU[2]
				
						if str(upperindx) in anglestackslist and options.saveanglestacks:	
							anglestackname = options.path + '/anglestack_' + str(upperindx) + '_angle' + str(int(upperangle)) + '.hdf' 
							refU.write_image( anglestackname, -1 )
				
			
			ptclNum += 1	
			
		else:
			print "Particle skipped because it's center in X, xc=%d, is outside the lower and upper boundaries to consider [%d,%d]; i.e., too far away from the tilt axis" % ( xc, xclowerthresh, xcupperthresh ) 
			print "First condition, xc > lowerthresh", int(xc) > int( xclowerthresh )
			print "First condition, xc < upperthresh", int(xc) < int( xcupperthresh )
		
			print "Both conditions", int(xc) > int( xclowerthresh ) and int(xc) < int( xcupperthresh )
			print "IF all true you shouldn't be reading this message!"
		
		
		
		ptclNum3D += 1
	
	E2end(logger)
	
	return




def checkcorners( img, options ):
	
	nx = img['nx']
	#ny = img['ny']
	
	#cornernx = round( nx * 0.05 ) #cornernx = round( nx * options.prunetest ) 
	
	#if not cornernx:
	
	cornernx = 8
	
	#cornerny
	#clipr = imgt.get_clip(Region(x,y, options.tilesize, options.tilesize))
	
	r1 = Region(0,0, cornernx, cornernx)
	corner1 = img.get_clip( r1 )
	
	r2 = Region(nx-cornernx,0, cornernx, cornernx)
	corner2 = img.get_clip( r2 )
	
	r3 = Region(nx-cornernx,nx-cornernx, cornernx, cornernx)
	corner3 = img.get_clip( r3 )
	
	r4 = Region(0,nx-cornernx, cornernx, cornernx)
	corner4 = img.get_clip( r4 )
	
	if not corner1['sigma']:
		print "\nimgsize is %d and sigma is %f" % (img['nx'], img['sigma'])
		print "region r1 is", r1
		print "empty sigma for corner1", corner1['sigma']
		return 0
	elif not corner2['sigma']:
		print "\nimgsize is %d and sigma is %f" % (img['nx'], img['sigma'])
		print "region r2 is", r2
		print "empty sigma for corner2", corner2['sigma']
		return 0
	elif not corner3['sigma']:
		print "\nimgsize is %d and sigma is %f" % (img['nx'], img['sigma'])
		print "region r3 is", r3
		print "empty sigma for corner3", corner3['sigma']
		return 0
	elif not corner4['sigma']:
		print "\nimgsize is %d and sigma is %f" % (img['nx'], img['sigma'])
		print "region r4 is", r4
		print "empty sigma for corner4", corner4['sigma']
		return 0
	else:
		return 1









def write2D( options, angle, icethickness, tomox, tomoy, xc, yc, zc, cumulativedx, cumulativedy, ref, apix, tag, ptclfile, maxtilt, sliceindx ):

	'''
	Extract an image immediately to the left and another immediately to the right of the zero tilt image for k=0; 
	then the next image to the left and the next to the right (in the tilt series) for k=1; so on and so forth...
	'''

	ret1 = extract2D( options, angle, icethickness, tomox, xc, yc, zc, cumulativedx, cumulativedy, sliceindx )
	img = ret1[0]
	fx = ret1[1]
	fy = ret1[2]
	finalimg = img.copy()
	
	if not finalimg['sigma']:
		print "\n(e2spt_subtilt)(write2D)ERROR: the extracted image is completely empty. Mean and sigma are", finalimg['mean'],finalimg['sigma']
		#sys.exit()
		return None
			
	if options.track:
		'''
		Preprocess the extracted images inside align2D, since preprocessing is needed only if alignment is performed.
		Start from k=0, then align them to the immediate previous image (k-1); for k=0, the k-1 image is the zero tilt image.
		(The zero tilt image is probably handled independently before calling this function, write2D).
		'''
		ret2 = align2D( options, ref, img )
		rdx = ret2[0]
		#rdy = ret2[1]
		rdy=0 	#particles shouldn't move at all in y
		
		kurtosis = ret2[2]
		
		#if float(kurtosis) > 1.0:
		#	print "kurtosis was > 1.0", kurtosis
		#	return None

		'''
		Reextract better centered images taking into account accumulated x and y shifts
		'''

		cumulativedx += rdx
		cumulativedy += rdy
	
		print "cumulativedx", cumulativedx
		print "cumulativedy", cumulativedy

		retf = extract2D( options, angle, icethickness, tomox, xc, yc, zc, cumulativedx, cumulativedy, sliceindx )
		
		
		#print "Slice index used for extraction was", sliceindx
	
		e = retf[0]
		finalimg = e.copy()
		
		fx = retf[1]
		fy = retf[2]
	
	threshy1 = float( options.excludeedge )
	threshx1 = float( options.excludeedge )
	
	threshy2 = float( tomoy ) - options.excludeedge
	threshx2 = float( tomox ) - options.excludeedge
	
	
	
	allgood = checkcorners( finalimg, options )	
	
	#if float( fx ) > threshx1 and float(xc) > threshx1 and float( fx ) < threshx2 and float (xc) < threshx2 and float( fy ) > threshy1 and float(yc) > threshy1 and float( fy ) < threshy2 and float(yc) < threshy2:
	if allgood and float( fx ) > threshx1 and float( fx ) < threshx2 and float( fy ) > threshy1 and float( fy ) < threshy2:
		
		finalimg['spt_tiltangle'] = angle
		finalimg['spt_tiltaxis'] = 'y'
		finalimg['spt_subtilt_x'] = fx
		finalimg['spt_subtilt_y'] = fy
		finalimg['origin_x'] = finalimg['nx']/2.0
		finalimg['origin_y'] = finalimg['ny']/2.0
		finalimg['origin_z'] = 0

		finalimg['apix_x'] = apix
		finalimg['apix_y'] = apix
		finalimg['apix_z'] = apix

		if options.coords3d:
			sptcoords = ( fx, fy, zc )
			finalimg['ptcl_source_coord']=sptcoords

		finalimg.process_inplace('normalize.edgemean')
	
	
		if options.subtractbackground and maxtilt:
			print "WARNING: \nBackground subtraction not working yet!"			
			#subtractBackground()	
	
		tmpimgfile = options.path + '/tmp.hdf'
	 
		finalimg.write_image( tmpimgfile, 0 )
	
		cmd = 'e2proc2d.py ' + tmpimgfile + ' ' + ptclfile + ' && rm ' + tmpimgfile
	
		if tag == 'lower':
			cmd = 'e2proc2d.py ' + ptclfile + ' ' + tmpimgfile + ' && mv ' + tmpimgfile + ' ' + ptclfile
	
		os.popen(cmd)
	
		return [finalimg, cumulativedx, cumulativedy]
	else:
		print "\nWARNING! Particle excluded from angle view %.2f since its center %.2f, %.2f is outside the --excludeedge limits x=[%.2f,%.2f] and y[%.2f,%.2f]" %(angle,fx,fy,threshx1,threshx2,threshy1,threshy2) 
		
		print "float( fx ) > threshx1", float( fx ) > threshx1
		print "float(xc) > threshx1", float(xc) > threshx1  
		print "float( fx ) < threshx2", float( fx ) < threshx2  
		print "float (xc) < threshx2", float (xc) < threshx2
		print "float( fy ) > threshy1", float( fy ) > threshy1
		print "float(yc) > threshy1", float(yc) > threshy1 
		print "float( fy ) < threshy2", float( fy ) < threshy2 
		print "float(yc) < threshy2", float(yc) < threshy2
		
		print "but 'allgood' is", allgood
		
		return None
		

def extract2D( options, angle, icethickness, tomox, xc, yc, zc, cumulativedx, cumulativedy, sliceindx ):
	
	
	print "(e2spt_subtilt)(extract2D) tomox is %d icethickness is %d" %( tomox, icethickness )
	
	if options.invertangles:
		angle *= -1
	
	tAxisShift = tomox/2.0
	xcToAxis = xc - tAxisShift
	oldx = xcToAxis
	
	
	#zSectionShift = tomoz/2.0
	zSectionShift = -1 * icethickness/2.0

	if options.zshift != 'half':
		#zSectionShift += icethickness/2.0 
		zSectionShift = int( options.zshift )
	
	zcToMidSection = zc + zSectionShift
	oldz = zcToMidSection

	
	'''
	Particles to the left of the tilt axis experience a POSITIVE shift in X when tilted, regardless of what the tilt angle is.
	Particles to the right of the tilt axis experience a NEGATIVE shift in X, regardless of what the tilt angle is.
	To account for this, we multiply times -1.
	'''
	cosTerm = xcToAxis * math.cos( math.radians(angle)  )
	
	'''
	Particles undergo additional movement in X upon tilting, due to their distance from the midZ plane of the tomogram.
	
	Particles ABOVE the middle Z plane, tilted NEGatively, are displaced NEGatively, to the left, compared to particles in the middle Z plane.
	Particles ABOVE the middle Z plane, tilted POSitively, are displaced POSitively, to the right, compared to particles in the middle Z plane.

	Particles BELOW the middle Z plane, tilted NEGatively, are displaced POSTitively, to the right, compared to particles in the middle Z plane.
	Particles BELOW the middle Z plane, tilted POSTitively, are displaced NEGatively, to the left, compared to particles in the middle Z plane.
	'''
	
	sinTerm = zcToMidSection * math.sin( math.radians(angle)  )
	
	xtToAxis = zcToMidSection * math.sin( math.radians(angle)  ) + 1*xcToAxis * math.cos( math.radians(angle)  )

	yt = yc

	xt = xtToAxis + tAxisShift
		
	
	
	'''
	alternatively, model the xc, yc, zc coordinates as a vector and just rotate in 3D,
	then retrieve the new xc2d and yc2d coordinates
	'''
	
	vect = numpy.array([xc - tAxisShift,0,zc-icethickness/2.0,0])
	newvector = yrotate( vect,angle )
	
	xc2d = int(newvector[0]) + tAxisShift
	#yc2d = int(newvector[1])
	yc2d = yt
	zc2d = int(newvector[2]) + icethickness/2.0
	
	
	
		
	newx = oldx * math.cos( math.radians(angle) ) + oldz * math.sin( math.radians( angle ) )
	
	newuncompensatedx = newx + tAxisShift
	
	
	daxis = cosTerm - oldx
	dmidz = sinTerm
	
	dx = xt -xc
	dx2 = daxis + dmidz
	
	
	#Translations in y should only occur when --track is on
	print "\nFor original x=%d, compensated x=%d, original z=%d, compensated z=%d, angle=%f, newx=%f, newUNcompx=%f, xt=%f, dx=%f, dx2=%f, daxis=%f, costerm=%f, dmidz=%f" %(xc,oldx, zc,oldz, angle, newx,newuncompensatedx,xt,dx,dx2,daxis,cosTerm,dmidz)
	
	print "yt %d, yc %d" %(yt,yc)

	if float(xt) < 0.0:
		print "Something went awfully wrong; you have a negative X coordinate",xt
		print "tomox and tomox/2.0 are", tomox, tomox/2.0

	#if float(yt) < 0.0:
	#	print "Something went awfully wrong; you have a negative Y coordinate",yt
	#	print "yc is", yc

	#if float(xt) < 0.0 or float(yt) < 0.0:
	#	print "Either X or Y are negative, see", xt, yt
	#	sys.exit()

	if float(xt) < float(options.boxsize)/2.0: #or float(yt) < float(options.boxsize)/2.0:
		print "Pick a smaller boxsize; otherwise, some of your particles will contain empty regions outside the image"
		print "Particle is centered at", xt, yc
		print "And boxsize/2 is", options.boxsize/2.0
		sys.exit()

	xt += cumulativedx
	
	#yt += cumulativedy		#particles shouldn't move at all in y
	
	##NEW
	#r = Region( (2*xc2d-options.boxsize)/2, (2*yc2d-options.boxsize)/2, sliceindx, options.boxsize, options.boxsize, 1)
	
	##OLD
	r = Region( (2*xt-options.boxsize)/2, (2*yt-options.boxsize)/2, sliceindx, options.boxsize, options.boxsize, 1)
	
	#print "\n\n\nRRRRRRRRRR\nThe region to extract is", r


	#print "\n(e2spt_subtilt.py) Extracting image for tilt angle %f and particle %d" %( angle, ptclNum )
	e = EMData()
	e.read_image(options.tiltseries,0,False,r)
	e.process_inplace('normalize.edgemean')
	
	print "\nold way", xt, yt
	print "new way", xc2d, yc2d
	
	##NEW
	#return [ e, xc2d, yc2d ]

	##OLD
	return [ e, xt, yt ]


def yrotate(M,theta):
	t = math.pi*theta/180
	cosT = math.cos( t )
	sinT = math.sin( t )
	#print 'M is', M
	R = numpy.array(
	[[ cosT,  0.0, sinT, 0.0 ],
	[ 0.0,   1.0,  0.0, 0.0 ],
	[-sinT,  0.0, cosT, 0.0 ],
	[ 0.0,  0.0,  0.0, 1.0 ]], dtype=numpy.float32)
	Mrot = numpy.dot(M ,R)

	return Mrot



def align2D( options, ref, img ):
	
	refp = ref.copy()
	imgp = img.copy()
	
	kurtosis = imgp.get_attr('kurtosis')
	#if kurtosis < 1:
	
	if options.lowpass or options.highpass or options.preprocess or options.threshold or options.mask or options.shrink > 1:
		refp = preprocImg( refp, options )
		imgp = preprocImg( imgp, options )
	
	#kurtosis = imgp.get_attr('kurtosis')
	
	imgp.process_inplace( 'filter.matchto',{'to':refp})
	
	ccf = refp.calc_ccf( imgp )
	ccf.process_inplace("xform.phaseorigin.tocorner") 
	
	locmax = ccf.calc_max_location()
							
	drecenterX = refp['nx']/2.0 - locmax[0]
	drecenterY = refp['ny']/2.0 - locmax[1]
	
	if options.shrink:
		print "(e2spt_subtilt)(align2D) on shrunk images, translation to recenter are", drecenterX, drecenterY

	if options.shrink:
		drecenterX *= options.shrink
		drecenterY *= options.shrink
	if options.shrink:
		print "(e2spt_subtilt)(align2D) on the actual images, translation to recenter are", drecenterX, drecenterY
	
	return [drecenterX,drecenterY,kurtosis]



def preprocImg( iimg, options ):
	
	img = iimg.copy()
	
	img.process_inplace('normalize.edgemean')
	
	if options.threshold and options.threshold != 'None' and options.threshold != 'none': 
		threshold=''
		try:
			threshold=parsemodopt(options.threshold)
		except:
			print "Failed to parse threshold"
		print "Parsed threshold is", threshold
		img.process_inplace( threshold[0], threshold[1] )

	
	if options.mask and options.mask != 'None' and options.mask != 'none':
		mask=''
		try:
			mask=parsemodopt(options.mask)
		except:
			pass
		img.process_inplace( mask[0], mask[1] )
	
	
	print "Raw image size", img['nx'],img['ny']
	if options.shrink and int(options.shrink) > 1:
		img.process_inplace('math.meanshrink',{'n': options.shrink })
	
	print "Shrunk image size", img['nx'],img['ny']
	
	if options.highpass and options.highpass != 'None' and options.highpass != 'none': 
		highpass=''
		try:
			highpass=parsemodopt(options.highpass)
		except:
			pass
		img.process_inplace( highpass[0], highpass[1] )
	
	
	if options.preprocess and options.preprocess != 'None' and options.preprocess != 'none': 
		preprocess=''
		try:
			preprocess=parsemodopt(options.preprocess)
		except:
			pass
		img.process_inplace( preprocess[0], preprocess[1] )
	
	
	if options.lowpass and options.lowpass != 'None' and options.lowpass != 'none': 
		lowpass=''
		try:
			lowpass=parsemodopt(options.lowpass)
		except:
			pass
		img.process_inplace( lowpass[0], lowpass[1] )
	
	#img.write_image( options.path + '/preproc.hdf',-1 )
	return img


def subtractBackground():
	
	'''
	Extract a large volume around each particle (only for k==0), to distinguish ptcl from background
	and generate background-substracted re-projections (for each tilt angle)
	'''
	if k == 0:			
		#ret = unbinned_extractor(options,bgboxsize,xc,yc,zc,options.cshrink,invert,center,options.tomogram)
		rw =  Region( (2*xc-bgboxsize)/2, (2*yc-bgboxsize)/2, (2*zc-bgboxsize)/2, bgboxsize, bgboxsize, bgboxsize)
	
		wholebox = EMData()
		wholebox.to_zero()
		wholebox.read_image(options.tomogram,0,False,rw)
	
		if int( options.shrink ) > 1:
			wholebox.process_inplace('math.meanshrink',{'n':options.shrink})
	
		#wholebox.process_inplace('normalize.edgemean')
	
		nsz = wholebox['sigma_nonzero']
	
		img.process_inplace("threshold.belowtozero",{"minval":snz*-1.5})
	
		#img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.5})
		#img.process_inplace("threshold.belowtozero",{"minval":snz/100.0})
	
		print "(e2spt_subtilt.py) Extracted whole 3D box " + str(angle) + " and mean " + str(wholebox['mean']) + " for particle " + str(ptclNum)
		wholebox.write_image(options.path + '/subtiltPtcl_' + str(ptclNum) + '_whole3D.hdf',0)

	if wholebox:
	
		#e.process_inplace('normalize')
	
		'''
		Rotate the larger box extracted and then clip into prism
		to avoid density gradient (a rotated cube rotates the data inside)
		causing higher density in projections along the diagonal of the cube
		'''
		#angle = angle *-1
		t= Transform({'type':'eman','az':-90,'alt':angle,'phi':90})
		#t= Transform({'type':'eman','az':90,'alt':angle})
		#t= Transform({'type':'eman','alt':angle})


		print "\nTransform to ROTATE volume is", t
	
		wholeboxRot = wholebox.copy()
		wholeboxRot.transform(t)
	
		finalbox = options.boxsize
		finalbgbox = bgboxsize
		if int( options.shrink ) > 1:
			finalbox = options.boxsize / options.shrink
			finalbgbox = bgboxsize / options.shrink
						
		rbgprism =  Region( (wholebox['nx'] - finalbox)/2, (wholebox['ny'] - finalbox)/2, (wholebox['nz'] - finalbgbox)/2, finalbox, finalbox, finalbgbox)
		wholeboxRot.clip_inplace( rbgprism )
		print "\nSizes of prism are", wholeboxRot['nx'],wholeboxRot['ny'],wholeboxRot['nz']
	
		if k == 0 :
			wholeboxRot.write_image(options.path + '/subtiltPtcl_' + str(ptclNum) + '_whole3DROT.hdf',0)
	
		ptclreprj = wholeboxRot.project("standard",Transform())
		print "\nGenerated ptclreprj with mean and XY sizes", ptclreprj['mean'],type(ptclreprj), ptclreprj['nx'],ptclreprj['ny'],ptclreprj['nz']
	
		ptclreprj.mult( math.cos( math.radians(angle) ) )
	
		mask = EMData(ptclreprj['nx'],ptclreprj['ny'],ptclreprj['nz'])
		mask.to_one()
		mask.process_inplace('mask.sharp',{'outer_radius':-2})
	
		ptclreprj.process_inplace('normalize.mask',{'mask':mask})
	
		if k==0:
			ptclreprj.write_image(options.path + 'reprj0nomatch.hdf',0)
		ptclreprj.process_inplace('filter.matchto',{'to':e})
	
		if k==0:
			ptclreprj.write_image(options.path + 'reprj0yesmatch.hdf',0)
		#ptclreprj.rotate(90,0,0)
		ptclreprj.write_image(options.path + '/subtiltPtcl_' + str(ptclNum) + '_reprj.hdf',outIndx)

	
		'''
		Generate projections of the background density by masking out the particle
		'''
		maskrad = finalbox/3.0
	
		bgbox = wholeboxRot.process('mask.sharp',{'inner_radius':maskrad})
							
		print "\nMasked bgbox with inner_radius, and ptclbox with outer_radius", maskrad
	
		bgprj = bgbox.project("standard",Transform())
	
		bgprj.mult( math.cos( math.radians(angle) ) )
	
		print "\nGenerated ptclreprj with mean and XY sizes", bgprj['mean'],type(bgprj), bgprj['nx'],bgprj['ny'],bgprj['nz']
	
		bgprj.process_inplace('normalize')
	
	
		bgprj.process_inplace('filter.matchto',{'to':e})
		bgprj.rotate(90,0,0)
		bgprj.write_image(options.path + '/subtiltPtcl_' + str(ptclNum) + '_bgprj.hdf',outIndx)
	
		clean = e - bgprj
		clean.write_image(options.path + '/subtiltPtcl_' + str(ptclNum) + '_clean.hdf',outIndx)
		print "\nComputed clean. Max e and Max bgprj are",e['maximum'], bgprj['maximum']

		cleanreprj = ptclreprj - bgprj
		print "\nComputed cleanprj"
	
		cleanreprj.write_image(options.path + '/subtiltPtcl_' + str(ptclNum) + '_cleanreprj.hdf',outIndx)

	outIndx+=1
	

	print "\n\n\n"	
	
	return


if __name__ == '__main__':
	
	main()




"""
PROGRAM SRCRAPS

		#		(cos q  0  -sin q   0)
		#Ry(q) = (0      1    0      0)
		#        (sin q  0  cos q    0)
		#        (0      0    0     1) 

		
		
		
	#	if float( options.shrink) > 1.0:					#Shrink them if the model will be smaller than the original tomgoram
	#		xc = xc/options.shrink
	#		yc = yc/options.shrink
	#		zc = zc/options.shrink
			
		#tomovol.set_value_at(xc,yc,zc,1)	#Set the value of the center pixel where any particles were picked to 1.
		
	
	#for tilt in tiltangles:
	#	t=Transform({'type':'eman','az':0,'phi':0,'alt':tilt})
	#	prj = tomovo.project('standard',t)
	#	
	#	for i in range(len(clines)):
	
	'''
	
	Old mathematical approach
	for j in range(nslices):
		
		2Dregion = Region(0,0,j,nx,ny,j+1)
		#a = EMData(s)
		slice=EMData()
		slice.read_image(options.tiltseries,0,False,2Dregion)
		tiltaxisx = int(nx)/2
		
		alpha = tiltangles[j]
		#k = 1
		for i in range(nptcls):
			#Some people might manually make ABERRANT coordinates files with commas, tabs, or more than once space in between coordinates
			clines[i] = clines[i].replace(", ",' ')	
			clines[i] = clines[i].replace(",",' ')
			clines[i] = clines[i].replace("x",'')
			clines[i] = clines[i].replace("y",'')
			clines[i] = clines[i].replace("z",'')
			clines[i] = clines[i].replace("=",'')
			clines[i] = clines[i].replace("_",' ')
			clines[i] = clines[i].replace("\n",' ')
			clines[i] = clines[i].replace("\t",' ')
			clines[i] = clines[i].replace("  ",' ')
			clines[i] = clines[i].split()		
		
			xc = int(clines[i][0])
			yc = int(clines[i][1])
			zc = int(clines[i][2])
		
			y2 = y
			
			xp2ta = tiltaxisx - xc
			zp = zc - options.thickness/2
			
			if alpha > 0:
				zp = -1 * zp
			
			dxa = (xp2ta + zp) * cos(alpha)
			
			x2 = tiltaxis - la
			
			r = Region((x2 - boxsize)/2,(y2 - boxsize)/2, boxsize, boxsize)
        		e = EMData()
			e.read_image(s,0,False,r)
			
			name = 'particle#' + str(k).zfill(len(pcoords)) + '_slice' + str(j).zfill(len(pcoords)) + '.mrc'
	'''


"""