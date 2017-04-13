#!/usr/bin/env python
#
# Author: Jesus Galaz, 11/01/2012; last update 31/Aug/2016
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
	
	'''
	parameters to be passed on to e2tomo_ctf.py
	'''
		
	parser.add_argument("--ampcont",type=float,default=0.05,help="""Default=0.05. Amplitude contrast to use for CTF correction phase flipping. Supply it to replace the value in ctfparamsfile(s), or if ctfparamsfile(s) are lacking altogether.""")
	
	parser.add_argument("--angles",type=str,default='',help="""Default=None. Overwrites --anglesfile. Comma separated list of tilt angles.""")
		
	parser.add_argument("--anglesfile",default='',type=str,help="""Default=None. File containing a list of tilt angles corresponding to the tilt angles of images 0 to n of an aligned tilt series""")
	
	parser.add_argument("--apix",type=float,default=0.0,help="""Default=whatever is on the header of the images. Sampling of the images in angstroms/pixel. Supply --apix here to replace the value in ctfparamsfile(s), or if ctfparamsfile(s) are lacking altogether.""")	
	
	parser.add_argument("--bfactor",type=int,default=1000,help="""Default=1000. Bfactor or temperature factor to use. Supply it to replace the value in ctfparamsfile(s), or if ctfparamsfile(s) are lacking altogether.""")
	
	parser.add_argument("--calcglobaldefocus", action='store_true', default=False, help="""Default=False. Calculate the global defocus. If --fitgradient is nor provided, these will be the final defocus values to use.""")
	
	parser.add_argument("--correctionwidth", type=int, default=0, help="""Default=tile size. Width of the strip to be phase-flipped. 1 would mean that the correction moves pixel by pixel.""")

	parser.add_argument("--cs", type=float,default=2.1,help="""Default=2.1. Cs of the microscope with which the images were collected. Supply it to replace the value in ctfparamsfile(s), or if ctfparamsfile(s) are lacking altogether.""")

	parser.add_argument("--ctfparamsfile",type=str,default='',help="""This should be a text file with ctf parameters in the following format; defocus=value voltage=value cs=value apix=value bfactor=value ampcont=value. A single space should separate each parameter from the next. Do not write any unit symbols for the values; just the numerical value. Defocus should be in microns, voltage in kV, apix in angstroms per pixel, and ampcont (amplitude contrast) should be a decimal; for example, 0.1 for 10 percent amplitude contrast. IF you want to use DIFFERENT PARAMETERS PER IMAGE, then the file must contain multiple rows with the different values. The first row will be used to phase flip the first image, the second row to phase flip the second, etc.""")
	
	parser.add_argument("--defocilist",type=str,default='',help='''Text file containing a single column of defocus values in microns. The file should have as many defocus values as images in the tiltseries or subtiltseries supplied.''')
	
	parser.add_argument("--depthfocus",type=float,default=0.0,help="""Default=0. Not used. Total variation in defocus (in Angstroms) tolerated within a strip to still consider it a region of 'constant defocus'.""")

	parser.add_argument("--defocusmin",type=float,default=0.0,help=""" If --autofit, minimum autofit defocus. Default=0.0, not used. A value will be estimated based on tilt angle and distance from the tilt axis.""")
	
	parser.add_argument("--defocusmax",type=float,default=0.0,help="""Default=0.0, not used. If --autofit, maximum autofit defocus. A value will be estimated based on tilt angle and distance from the tilt axis.""")
	
	parser.add_argument("--exclude",type=str,default='',help="""Default=None. Only works if --tiltseries is provided. Comma-separated list of image indexes in the --tiltseries to exclude from CTF fitting. For example, --exclude 0,3,4,6,7.""")
	
	parser.add_argument("--excludeedges",action='store_true',default=False,help='''Ignore 'excedent' (smaller than the width of a strip) at the edge of micrographs after dividing them into strips.''')
	
	parser.add_argument("--fitgradient", action='store_true', default=False,help="""Runs automated CTF fitting on the input images, based on tiling.""")
	
	parser.add_argument("--fixctfhighpass", action='store_true', default=False,help="""Applies the 'filter.ctfcorr.simple' processor. Type 'e2help.py processors' at the commnadline for an explanation.""")
	
	parser.add_argument("--imagestem",type=str,default='',help="""Default=None. If the images to apply ctf correction on are already unstacked and are individual mrc files, supply a common string to all of them.""")
	
	parser.add_argument("--infodir",type=str,default='',help="""Folder typically produced by e2evalimage.py or previous runs of this program containing info.json files, one per tilt image in a tilt series. Each .json file should contain the fitted ctf and all associated parameters for each tilt image.""")

	parser.add_argument("--invert",action='store_true',default=False,help='''Invert the contrast of the output data, compared to the input data.''')
	
	parser.add_argument("--mintiles", type=int, default=0, help="""Minimum number of 'good tiles' in strip to consider it.""")

	parser.add_argument("--overlaptiles",action='store_true',default=False,help="""Default=False. If provided, it will cause tiles to overlap by 50 percent in x and y for power spectrum computation (periodogram averaging).""")

	parser.add_argument("--path",type=str,default='sptctf',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'sptctf'; for example, sptctf_02 will be the directory by default if 'sptctf_01' already exists.""")
					
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	parser.add_argument("--phaseflipwhole", action='store_true',default=False,help="""This will perform phase flipping on the entire image for each image in an aligned tiltseries using the CTF parameters supplied.""")
	
	parser.add_argument("--phaseflipstrips",action='store_true',default=False,help="""This will perform phase flipping on images of an aligned tilt series on a strip-by-strip basis, assuming the supplied ctf parameters correspond to the proper values at the tilt axis, either the same values for all images (--defocus,--ampcont,--cs,--apix,--voltage,--bfactor) or a different set for each (--ctfparamsfile), taking into account the tilt angle for  each image (--anglesfile), which should be supplied through an IMOD-like .tlt file.""")
		
	parser.add_argument("--predictdefocus", action='store_true', default=False, help="""Default=False. Works only when --autofitting is supplied and --calcglobaldefocus is NOT supplied. This constrains defocus search during strip-based defocus gradient fitting around the predicted value based on tilt angle and global defocus (the average defocus of the entire micrograph).""")
	
	parser.add_argument("--prunetest",type=float,default=0.1,help="""Default=0.1. Decimal number that indicates the percentage of --tilesize (in terms of side length) to tolerate of 'bad' values (i.e., empty regions of constant density) at the corners, and still include the tile for CTF fitting. For example, if --tilesize=256, and --prunetest=0.1, a box of ~25-26 pixels each corner of every tile will be analyzed and if the standard deviation of any of the corners is 0, the tile will be excluded. To turn off this option supply --prunetest=-1.0. The program automatically adjusts things so that the minimum size of regions at the corners to check will be 4x4 pixels.""")
	
	parser.add_argument("--radius",type=int,default=0,help="""Radius of the particle in pixels.""")
	
	parser.add_argument("--savestriptiles",action='store_true',default=False,help="""Saves all tiles for all strips, for all images, in one stack per strip.""")

	parser.add_argument("--saveffts",action='store_true',default=False,help="""Saves ffts of each average of tiles per strip, for all images.""")
				
	parser.add_argument("--stripflipstep",type=int,default=0,help="""Default=0. Not used. This is automatically calculated based on the allowable depth of focus, in turn calculated based on apix, voltage, angle, thickness, etc. If this parameter is on, it will determine the strips in each image to flip using the same defocus. To get flipping with no "seams" (flips pixel by pixel) you would provide a value of 1 for this parameter.""")
	
	parser.add_argument("--skipimgstrips", type=str, default='', help="""Default=None. Comma-separated list of image indexes to exclude from strip-based fitting (in this case, only global defocus tiling the entire image wil be measured for those images).""")
	
	parser.add_argument('--subset', type=int, default=0, help='''Requires --subtiltsdir. Specify how many subtiltseries (or particles) to ctf correct. If you specify 10, the first 10 subtiltseires in --subtiltsdir will be corrected. 0 means "process all" because it makes no sense to process none''')
	
	parser.add_argument("--subtiltsdir",type=str,default='',help="""Provide a directory containing individual stacks, where each stack is a 'mini tilt series' or a 'subtilt series' for single particles. Then, each image for each particle in the dir will be phase-phlipped using the ctf parameters you provide. If each image in the subtilt series is at a different defocus, then the parameters should be provided through --ctfparamsfile, whith a different defocus value per row. (There should be as many rows as images in each subtiltseries).""")
	
	parser.add_argument("--thickness", type=int,default=0,help="""This corresponds to the spread of the specimen in Z in the tomogram, at the same binning (sampling) as the provided tiltseries, images or subtiltseries.""")	
	
	#parser.add_argument("--thicknessauto",action='store_true',default=False,help="""Default=False. Requires --radius. The thickness of the specimen will be --radius*2. If --coords is provided (and if --thickness is *not* provided), the thickness of the specimen in Z will be calculated by computing the difference between the largest and the smallest Z coordinate in --coords, plus the size of the specimen as --radius*2.""")
	
	parser.add_argument("--tilesize",type=int,default=512,help="""Tile size to use for strips when --autofit is provided.""")
		
	parser.add_argument("--tiltseries", type=str, default='', help="""Aligned tilt series. File format must be MRC and must have .mrc or .st or .ali extension.""")
					
	parser.add_argument("--reconstructor", type=str,default="fourier:mode=gauss_2",help="""Default=fourier:mode=gauss_2. The reconstructor to use to reconstruct the tilt series into a tomogram. Type 'e2help.py reconstructors' at the command line to see all options and parameters available. To specify the interpolation scheme for the fourier reconstruction, specify 'mode'. Options are 'nearest_neighbor', 'gauss_2', 'gauss_3', 'gauss_5', 'gauss_5_slow', 'gypergeom_5', 'experimental'. For example --reconstructor=fourier:mode=gauss_5 """)
	
	parser.add_argument("--targetdefocus", type=float,default=0.0,help="""Default=0 (not used). Target defocus at the tilt axis. In the absence of ctfparamsfile(s) this value will be assumed to be the defocus at the tilt axis for all tilt images and will be used to constrain calculation of the global defocus of each image.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")
	
	parser.add_argument("--voltage", type=int,default=200,help="""Default=200. Voltage of the microscope with which the images where collected. Supply it to replace the value in ctfparamsfile(s), or if ctfparamsfile(s) are lacking altogether.""")

	
	'''
	parameters exclusive to e2stp_ctf.py
	'''
	#parser.add_argument("--coords",default='',type=str,help="""Default=None. Text file containing x y z (or just z) coordinates for the particles, used to calculate thickness if --thicknessauto is specified for ctf fitting.""")
	
	#parser.add_argument("--defocusbottom",action='store_true',default=False,help="""Requires --coords. Assumes the signal for defocus measurement (e.g., carbon film) is at the top layer of the tomogram. By default, the average signal is assumed to correspond to the middle Z plane of the tomogram.""")
	
	#parser.add_argument("--defocustop",action='store_true',default=False,help="""Requires --coords. Assumes the signal for defocus measurement (e.g., carbon film) is at the top layer of the tomogram. By default, the average signal is assumed to correspond to the middle Z plane of the tomogram.""")
	
	#read this from the header?	
	#parser.add_argument("--framexsize",type=int,default=0,help="""This correspond to the X size in pixes of the images/frames in the raw tilt series; that is, the size of the entire frame along the X axis (perpendicular to the direction of the tilt axis in the aligned tilt series). It is used to calculate the distance of each particle (subtiltseries) to the tilt axis, since this will induce different shifts in defocus in 3-D for the actual particles. Particles right at the tilt axis don't move "up" or "down" as they are tilted.""")

	#parser.add_argument("--pad2d", type=float,default=0.0,help="""Default=0. Padding factor to zero-pad the 2d images in the tilt series prior to reconstruction. (The final reconstructed subvolumes will be cropped to the original size).""")

	#parser.add_argument("--pad3d", type=float,default=0.0,help="""Padding factor to zero-pad the reconstruction volume. (The final reconstructed subvolumes will be cropped to the original size).""")	
	
	#parser.add_argument("--save2d",action='store_true',default=False,help="""If on, the CTF corrected subtiltseries will be saved as 2-D image stacks [one per particle].""")
	
	#parser.add_argument("--nozcorrection",action='store_true',default=False,help="""If you turn on this option and --subtiltsdir is provided, the position in Z of each subtomogram will not be considered for CTF correction""")
		
	'''
	deprecated parameters
	'''
	#parser.add_argument("--save3d",action='store_true',default=False,help="""If on, the CTF corrected subtiltseries will be reconstrcuted into subvolumes and saved into a stack. Options --reconstructor, --pad2d, --pad3d are used if --save3d is on.""")
	#parser.add_argument("--outputstem",type=str,default='',help="""Stem common to all output image stacks. For example, if --outputstem=myvirus and --save2d is provided, the phase-flipped images for each subtiltseries will be saved to myvirus_subtiltptclXXXX.hdf. The stack of reconstructed subvolumes will be saved to myvirus_stack3d.hdf""")
	#parser.add_argument("--stripstep",type=int,default=0,help="""Default=0. Not used. This is automatically calculated based on strip width, which is automatically calculated based on the allowable depth of focus, in turn calculated based on apix, voltage, angle, thickness, etc. If this parameter is on, it will determine the number of strips and the overlap between them for defocus gradient fitting. For example, for a 4000x4000 pixels image, a strip step of 400 would yield 10 strips, by default. If --stripstep=1 were provided, the image would be devided into 4000-400=3600 strips. The first strip would go from pixel 0 to pixel 400, the second strip from pixel 1 to pixel 401, the third from pixel 2 to 402, etc... up to the las strip going from pixel 3600 to 4000.""")

	(options, args) = parser.parse_args()
	
	'''
	c:detect errors to prevent crashes, e.g. from having defective files
	'''
	errordetector(options)

	'''
	c:if no errors detected thus far, log current run of the program
	'''
	logger = E2init(sys.argv, options.ppid)
		
	'''
	c:if no crashes so far, make the directory where results will be stored
	'''
	from e2spt_classaverage import sptmakepath
	options = sptmakepath (options, 'sptctf')
	
	'''
	c:store parameters defined for this run in a text file for easy review (kind of redundant with logging, 
	c:but not everyone is aware of the invisible log file; plus, the logfile can grow with runs of OTHER programs)
	'''
	from e2spt_classaverage import writeParameters
	cmdwp = writeParameters(options,'e2spt_ctf.py', 'sptctf')
	
	'''
	c:if any images are to be excluded, prune them out from the start and redefine a new, clean --tiltseries and --anglesfile 
	'''
	if options.exclude:
		options.tiltseries,options.anglesfile = pruneexcluded(options)
	
	'''
	c:get the tiltangles into an indexed and ordered dictionary
	'''
	#angles = {}
	#if options.anglesfile:
	angles = getangles( options )
	print "returned angles are", angles, len(angles)
	
	'''
	c:get the images to process into an ordered list
	'''
	imagefilenameslist = getimages(options)

	'''
	c:make sure the number of angles and number of images match
	'''
	nangles = len( angles )
	nimgs = len( imagefilenameslist )
	if nangles != nimgs:
		print "\nERROR: The number of angles %d does not coincide with number of images %d" % ( nangles, nimgs)
		sys.exit(1)	

	'''
	c:pair images and angles into and ordered dictionary
	'''
	imagefilenamesdict = imgangpairs(options,angles,imagefilenameslist)

	'''
	c:downsize to a subset if so specified
	'''
	if options.subset:
		imagefilenamesdict = { key:value for key, value in imagefilenamesdict.items() if key < options.subset }
		imagefilenamesdict = collections.OrderedDict(sorted(imagefilenamesdict.items()))
		
		angles = { key:value for key, value in angles.items() if key < options.subset }
		angles = collections.OrderedDict(sorted(angles.items()))
	
		print "subset of imagefilenamesdict",imagefilenamesdict
		print "subset of angles",angles
	
	'''
	c:figure out apix to use, the x size of the images, and the number of images nimgs to process.
	c:if not explicitly provided through the command line
	'''
	apix = 0.0
	framexsize = 0
	nimgs = len(imagefilenamesdict)
					
	if not options.apix:# or not options.framexsize:

		count=0
		for indx in imagefilenamesdict:
			if count < 1:
				apix = EMData( imagefilenamesdict[indx][0], 0, True )['apix_x']
				#framexsize = EMData( imagefilenamesdict[indx][0], 0, True )['nx']
			count+=1

	if options.apix:
		apix = options.apix
		
	#if options.framexsize:
	#	framexsize = options.framexsize
	#'''
	#c:calculate the specimen thickness value to use
	#'''
	#thicknesscalc(options)

	'''
	c:read or generate the CTF parameters to use
	'''
	ctfs = genctfparamlines( options, apix, nimgs, angles, imagefilenamesdict )
	
	'''
	c:if a file with defocuses to use was not provided, and the defocus was not fitted,
	c:then a defocus file will be generated from the defocus values in the ctfs, which 
	c:should be the --targetdefocus value provided
	'''
	if not options.defocilist and not options.fitgradient and not options.calcglobaldefocus:
		defocusesfile = options.path + '/defocuses.txt'
		os.system( 'touch ' + defocusesfile )
		fd = open(defocusesfile,'w')
		linesd = []
		for indx in angles:
			angle = angles[indx]
			defocus = ctfs[indx].defocus
				
			lined = str(angle)+ '\t' +str(defocus) + '\n'
			linesd.append( lined )
			
		fd.writelines(linesd)
		fd.close()
	
	#print "ctfs len is", len(ctfs)
	#print "and ctfs are" 
	#for ctf in ctfs:
	#	print ctf, ctfs[ctf]
	

	'''
	#c:verify that you have CTF parameters for each image, returned from genctfparamlines
	'''
	if ctfs:
		if len( ctfs ) != len( imagefilenamesdict ):
			print """\n(e2spt_tomo)(main) ERROR: it seems like you have fewer parameter
				lines in the --ctfparamsfile than images in the tilt series.
				You need one line of parameters per image.
				To apply the same correction to all images, enter the parameters directly,
				through --defocus, --apix, --ampcont, --bfactor, --cs and --voltage"""
			sys.exit(1)
	else:
		print """\n(e2spt_tomo)(main) ERROR: there is no CTF information for any of the images. If you were
			using --autofit, this means autofitting failed for all images."""
		sys.exit(1)
			
	'''
	#c:once the images are organized, proceed to CTF correct them, depending on the strategy selected
	'''
	if options.phaseflipwhole:
		'''
		#c:apply the same correction to the entire image. Only works on low tilt or high mag, 
		#c:when the defocus gradient is irrelevant.
		'''
		flipstack(options,imagefilenamesdict,ctfs)

	elif options.phaseflipstrips:
		
		flipstripbystrip( options, imagefilenamesdict, ctfs, angles )
	
	E2end(logger)
	
	print "ABCDEFG"
	
	return


def getimages(options):
	'''
	ABERRANT Linux does NOT order lists by default, whereas OSX seemingly does. Find filenames FIRST, then explicitly ORDER them; 
	eventually ordered filenames will be paired with ordered tilt angles
	'''

	imagefilenameslist = [] 
	
	'''
	#If input consists of individual image files, find them and put them into an imagefilenames dictionary
	'''
	if options.imagestem:
		print "\n(e2spt_ctf)(getimages)processing all images in the current directory containing the following string", options.imagestem
		findir=os.listdir(os.getcwd())
		
		for f in findir:
			if options.imagestem in f:
				imagestem = f.replace('.mrc','')
				imagefilenameslist.append(f)
			
		print "\n(e2spt_ctf)(getimages) found these many tilt images",len(imagefilenameslist)

	'''
	#If input is a tiltseries, unstack it, then put the individual images into an imagefilenames dictionary
	'''
	if options.tiltseries:
		#nz=EMData(options.tiltseries,0,True)['nz']
		if '.st' in options.tiltseries or '.mrc' in options.tiltseries or '.mrcs' in options.tiltseries or '.ali' in options.tiltseries or '.hdf' in options.tiltseries:
			print "\n(e2spt_ctf)(getimages) unstacking tiltseries", options.tiltseries
			
			nimgs = EMData(options.tiltseries,0,True)['nz']
			print "\n(e2spt_ctf.py)(main) containing %d tilt images" %(nimgs)
			
			#if nimgs > 1:
			
			#c:make sure the unstacked tiltseries images will have .hdf extension				
			outname = options.tiltseries

			outname = 'tomoctftmp.hdf'
			
			'''
			if '.mrc' in options.tiltseries[-4:]:
				outname = options.tiltseries.replace('.mrc','.hdf')
			if '.mrcs' in options.tiltseries[-5:]:
				outname = options.tiltseries.replace('.mrcs','.hdf')
			if '.st' in options.tiltseries[-3:]:
				outname = options.tiltseries.replace('.st','.hdf')	
			if '.ali' in options.tiltseries[-4:]:
				outname = options.tiltseries.replace('.ali','.hdf')
			if '.MRC' in options.tiltseries[-4:]:
				outname = options.tiltseries.replace('.MRC','.hdf')
			if '.MRCS' in options.tiltseries[-5:]:
				outname = options.tiltseries.replace('.MRCS','.hdf')
			if '.ST' in options.tiltseries[-3:]:
				outname = options.tiltseries.replace('.ST','.hdf')	
			if '.ALI' in options.tiltseries[-4:]:
				outname =  options.tiltseries.replace('.ALI','.hdf')	
			

			outname = outname.replace('.hdf','_UNSTACKED_sptctftmp.hdf')
			'''
			
			if options.path not in outname:
				outname = options.path + '/' + outname
			
			#c:define the command to unstack the tiltseries and execute it
			cmdun = 'e2proc2d.py ' + options.tiltseries + ' ' + outname

			if nimgs > 1:
				cmdun += ' --unstacking'
			
			#print "\n\n.st in tiltseries", ".st" in options.tiltseries[-3:]
			
			print "\nunstacking cmdun is", cmdun

			
			runcmd( options, cmdun )
			
			#c:images are unstacked directly into the output directory; and add them to imagefilenameslist as the images to work with
			outnamestem = os.path.basename(outname).replace('.hdf','')
			
			c = os.getcwd()
			findir = os.listdir( options.path )
			#print "\nstem to look for in unstacked images is",outnamestem
			for f in findir:
				if outnamestem in f:
					#print "found image", f
					finalimagefile=options.path + '/' + f
					imagefilenameslist.append( finalimagefile )
					#print "and appended it to imagefilenameslist"				
			
			#elif nimgs == 1:
			#	imagefilenameslist.append( options.tiltseries )

			nfiles = len(imagefilenameslist)
			if nimgs != nfiles:
				print """\n(e2spt_ctf.py)(getimages) ERROR: it seems like not all the images
				in the tilt series were properly unstacked. nimages and nfiles are""",nimgs,nfiles
				sys.exit(1)
			
		else:
			print "\n(e2spt_ctf)(getimages) --tiltseries must be in .st or .mrc or .mrcs or .ali or .hdf extension"
			sys.exit(1)	
		
		
	'''
	sort the list of images to avoid cross platform issues stemming from differences in default ordering;
	then, fill in an indexed dictionary
	'''
	imagefilenameslist.sort()

	return imagefilenameslist


def imgangpairs(options,angles,imagefilenameslist):

	imagefilenamesdict = {}
	kki=0
	#print "\nimagefilenameslist is", imagefilenameslist
	for f in imagefilenameslist:
		#imagestem = ''
		#if '.mrc' in f[-4:]:
		#	imagestem = f.replace('.mrc','')
		#elif '.MRC' in f[-4:]:
		#	imagestem = f.replace('.MRC','')
		#elif '.st' in f[-3:]:
		#	imagestem = f.replace('.st','')
		#elif '.ST' in f[-3:]:
		#	imagestem = f.replace('.ST','')
		#elif '.ali' in f[-4:]:
		#	imagestem = f.replace('.ali','')
		#elif '.ALI' in f[-4:]:
		#	imagestem = f.replace('.ALI','')
		#elif '.hdf' in f[-4:]:
		#	imagestem = f.replace('.hdf','')
		#elif '.HDF' in f[-4:]:
		#	imagestem = f.replace('.HDF','')

		#if imagestem:
		
		imagefilenamesdict.update({kki:[ f,angles[kki] ]})
		kki+=1

	imagefilenamesdict = collections.OrderedDict(sorted(imagefilenamesdict.items()))
	
	return imagefilenamesdict

'''
def thicknesscalc(options):
	print "\n(e2tomo_ctf)(thicknesscalc) determining thickness"
	if not options.thickness:
		#print "--thickness not provided"
		if options.thicknessauto:
			#print "--thicknessauto is on"
			if options.coords:
				f=open(options.coords,'r')
				lines=f.readlines()
				zs = [ int( line.replace('\t',' ').split(' ')[-1].replace('\n','') ) for line in lines ]
				f.close()
				
				maxz = max( zs ) 
				minz =  min( zs ) 
				
				autoicethickness = maxz - minz
				
				if options.radius:
					autoicethickness += 2*options.radius
				else:
					print "\nWARNING: --radius was not provided and therefore --thicknessauto will be calculated solely based on the particle centers distribution."

				icefile = options.path+'/autoicethickness.txt'
				f = open( icefile, 'w')
				
				line = [str(autoicethickness)+'\n']
				f.writelines(line)
				f.close()

				options.thickness = autoicethickness
						
			elif not options.coords:
				#print "--coords not provided"
				if options.radius:
					#print "--radius is", options.radius 
					options.thickness = options.radius*2
					
				else:
					print "\nERROR: --thicknessauto requires --coords and/or --radius (ideally both)."

		elif not options.thicknessauto:
			if options.radius:
				options.thickness == options.radius*2
			elif not options.radius:
				print "\nWARNING: no --thickness or --radius to calculate it with provided, and --thicknessauto is also turned off."

	elif options.thickness:
		if options.radius:
			if options.thickness < options.radius*2:
				options.thickness == options.radius*2
				print "\nWARNING: --thickness was smaller than the diameter of the specimen (--radius *2), which makes no sense. Resetting --thickness to be --radius*2."
	
	print "returning thickness value", options.thickness
	#sys.exit(1)
	return
'''

def pruneexcluded(options):
	#cmdex = 'e2proc2d.py ' + options.tiltseries + ' ' + newtiltseries + ' --outmode int16'
	unstackedimgsbase = os.path.splitext(options.tiltseries)[0] + '_sptctftmp.hdf'
	cmdun = 'e2proc2d.py ' + options.tiltseries + ' ' + unstackedimgsbase + ' --unstacking'

	runcmd(options,cmdun)
	
	excludefile = options.path + '/excluded.lst'
	f = open( excludefile,'w')
	excludeindxs = options.exclude.split(',')
	lines = [ line + '\n' for line in excludeindxs ]
	f.writelines( lines )
	f.close()
	
	#cmdex +=  ' --exclude ' + excludefile
	
	stem=os.path.splitext(options.tiltseries)[0] + '_sptctftmp' 
		
	findir=os.listdir( os.getcwd() )
	imgslist=[]
	imgsdict={}

	i=0
	for fi in findir:
		if stem in fi:
			#imgsdict.update( { i:fi } )
			imgslist.append(fi)
	
	imgslist.sort()
	 
	imgskeep=[]
	
	for i in range(len(imgslist)):
		if str(i) not in excludeindxs:
			imgskeep.append( imgslist[i] )
		elif str(i) in excludeindxs:
			os.remove( imgslist[i] )
	
	newtiltseries = options.path + '/' + os.path.splitext(options.tiltseries)[0] + '_sptctftmp' + os.path.splitext(options.tiltseries)[1]
	
	cmdstack = 'e2buildstacks.py --stackname tmpstack.hdf ' + stem + '* && e2proc2d.py tmpstack.hdf ' + newtiltseries + ' --twod2threed --outmode int16'
	print "command to rebuild tiltseries", cmdstack
	runcmd(options,cmdstack)
	
	#os.rename(newtiltseries, options.path + '/' + newtiltseries)
	
	anglestmp=getangles(options)
	tltlines=[]
	for i in anglestmp:
		if str(i) not in excludeindxs:
			lin=str(anglestmp[i])+'\n'
			tltlines.append(lin)
	
	#tltlines.sort()
	
	newtltfile = options.path + '/' + os.path.splitext(options.anglesfile)[0] + '_sptctftmp' + os.path.splitext(options.anglesfile)[1]
	g=open(newtltfile,'w')
	g.writelines(tltlines)
	g.close()

	findir = os.listdir(os.getcwd())
	for fi in findir:
		if '_sptctftmp' in fi:
			os.remove(fi)
	
	os.remove('tmpstack.hdf')
	
	return newtiltseries,newtltfile


'''
c:function to run commands and the command line
'''
def errordetector( options ):
	
	if options.tiltseries and options.imagestem:
		print """ERROR: You either 1) supply a tiltseries with .mrc, .st, or .ali extension (in MRC format),
		or with .hdf extension (in HDF format), for an entire tomogram, 2) a stem (a 'string') common 
		to all individual .mrc or .hdf images corresponding to a tiltseries, OR 3) a directory 
		with subtiltseries in .hdf format for individual subtomograms. You cannot supply both
		--tiltseries, --subtiltsdir and --imagestem at the same time. Pick one."""
		sys.exit()
	
	nimgs=0
	if options.tiltseries:
		try:
			hdr=EMData(options.tiltseries,0,True)
			nimgs=hdr['nz']
		except:
			print "\n(e2tomo_ctf)(errordetector) ERROR: bad image", options.tiltseries
			sys.exit()	
	
	if options.phaseflipstrips:
		if not options.fitgradient:
			print "\n(e2tomo_ctf)(errordetector) ERROR: --phaseflipstrips requires --fitgradient"
	
	return


'''
c:function to run commands and the command line
'''
def runcmd( options, cmd ):
	if options.verbose > 9:
		print "\n(e2spt_ctf)(runcmd) running command", cmd
	
	p=subprocess.Popen( cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	text=p.communicate()	
	p.stdout.close()
	
	if options.verbose > 9:
		print "\n(e2spt_ctf)(runcmd) done"
	return


'''
c:function to determine the ctf to use for each image at the tilt axis, depending on input parameters
'''
def genctfparamlines( options, apix, nimgs, angles, imagefilenames ):
	
	#indxstoexclude = [int(i) for i in options.exclude.split(',')]
	
	print "e2spt_ctf (genctfparamlines)"
	#print "received imagefilenames", imagefilenames
	
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
				ctfs.update( { kk:ctf } )
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
				ctfs.update( { kk:ctf } )
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
			print "\nread defocuses are",defoci
			for d in defoci:
				#if kk not in indxstoexclude and kk < len(angles):
				#de = d.replace('\n','').replace(' ','').replace('\t','')
				de=d.split()[-1]
				line = 'defocus=' + str( de ) + ' ampcont=' + str( options.ampcont ) + ' voltage=' + str( options.voltage ) + ' cs=' + str( options.cs ) + ' apix=' + str( apix ) + ' bfactor=' + str( options.bfactor )
				ctf = ctfparamparser( line ) 
				angle = angles[ kk ]
				ctfs.update( { kk:ctf } )
				kk+=1
		else:
			print """\nERROR: There's nothing to do. If you supply --defocilist, you also 
			have to provide the following 4 parameters:
			--ampcont,--cs,--voltage,--bfactor (and, optionally, --apix, if the images 
			don't have the correct apix in their headers.)"""
			sys.exit(1)
		
	else:
		print "(e2spt_ctf)(main) autofitting using voltage=%.2f, cs=%.2f, apix=%.2f, ampcont=%.2f" %( float(options.voltage), float( options.cs), float( apix), float( options.ampcont))
		#print "imagefilenames", imagefilenames
		if options.voltage and options.cs and apix and options.ampcont:
			defocusmean=0 
			defocusstd=0
			finaldefoci={}
			if options.calcglobaldefocus:
				ret = calcglobaldefocus( options, apix, imagefilenames, angles )
				finaldefoci = ret[0]
				defocusmean = ret[1]
				defocusstd = ret[2]
			elif options.targetdefocus:
				for i in range(nimgs):
					finaldefoci.update({ i:options.targetdefocus })
			
			if options.fitgradient:
				finaldefoci = fitdefocusgradient( options, apix, imagefilenames, angles, finaldefoci, defocusmean, defocusstd )
			
			if finaldefoci:
				for kk in angles:
					#print '\nkk is', kk
					angle=angles[kk]
					#print 'angle is', angle
					#print "len angles", len(angles)
					#print 'len finaldefoci is', len(finaldefoci)
					
					print 'defocus, finaldefoci[kk]',finaldefoci[kk]
					params = {'ampcont':options.ampcont,'apix':apix,'bfactor':options.bfactor,'cs':options.cs,'defocus':finaldefoci[kk],'voltage':options.voltage}	
					ctf = EMAN2Ctf()
					ctf.from_dict(params)
					ctfs.update( {kk:ctf} )
				
		else:
			if not options.voltage:
				print "\n(e2spt_ctf)(main) ERROR: --voltage required with --fitgradient."
			if not options.cs:
				print "\n(e2spt_ctf)(main) ERROR: --cs and required with --fitgradient."
			if not options.ampcont:
				print "\n(e2spt_ctf)(main) ERROR: --ampcont required with --fitgradient."
			if not apix:
				print "\n(e2spt_ctf)(main) ERROR: apix required with --fitgradient."
			sys.exit(1)
			
	return ctfs


'''#
#c:function to read angles and return them in an indexed, ordered dictionary
#'''
def getangles( options ):
	
	anglesdict = {}
	
	lines=[]

	if options.anglesfile:
		if options.verbose > 9:
			print "\n(e2spt_ctf)(getangles) reading angles from --anglesfile %s" %(options.anglesfile)
		f = open( options.anglesfile, 'r' )
		lines = f.readlines()
		f.close()
		
		print "lines are len(lines)", lines, len(lines)
		
		ai=0
		for line in lines:
			line = line.replace('\t','').replace('\n','')
	
			if line:
				anglesdict.update( {ai:float(line) } )
			ai+=1

	elif options.angles:
		if options.verbose > 9:
			print "\n(e2spt_ctf)(getangles) parsing angles from --angles %s" %(options.anglesfile)
		lines=[ float(a) for a in options.angles.split(',')]
		
		for i in range(len(lines)):
			anglesdict.update( {i:lines[i]} )
	
	if options.verbose > 9:
		print "\n(e2spt_ctf.py)(getangles) angles are", anglesdict
	
	anglesdict = collections.OrderedDict(sorted(anglesdict.items()))
	
	return anglesdict
	
	
def ctfparamparser( pline ):
	
	defocus = pline.replace('\n',' ').split("defocus=")[-1].split(' ')[0]
	ampcont = pline.replace('\n',' ').split("ampcont=")[-1].split(' ')[0]
	voltage = pline.replace('\n',' ').split("voltage=")[-1].split(' ')[0]
	cs = pline.replace('\n',' ').split("cs=")[-1].split(' ')[0]
	apix = pline.replace('\n',' ').split("apix=")[-1].split(' ')[0]
	bfactor = pline.replace('\n',' ').split("bfactor=")[-1].split(' ')[0]

	params = {'ampcont':ampcont,'apix':apix,'bfactor':bfactor,'cs':cs,'defocus':defocus,'voltage':voltage}
	print "\nparameters are",params
	print "\n(e2spt_ctf.py)(ctfparamparser) The parsed parameters are:\n"
	for key in params.keys():
		print key + '=' + params[key] +'\n'
	
	ctf = EMAN2Ctf()
	#ctf.from_dict({'defocus':params['defocus'],'bfactor':params['bfactor'],'ampcont':params['ampcont'],'apix':params['apix'],'voltage':params['voltage'],'cs':params['cs']})	
	ctf.from_dict(params)
	
	return ctf


def flipstack(options,imgsdict,ctfs,apix):
	
	for indx in imgsdict:
		print "\n(e2tomo_ctf)(flipstack) working on image",indx,imgsdict[indx]
	
		print "\nI will phaseflip these images ignoring the defocus gradient",len(imgsdict)
		
		ctf = ctfs[indx]
		
		if ctf:
			imgfile = imgsdict[indx][0]
			imgfilebase = os.path.basename(imgsdict[indx][0])
			imgfilextension = os.path.splitext(imgfilebase)[1]
			print "\nLoading image file", imgfile
			img = EMData( imgfile )
		
			originalx = img['nx']
			originaly = img['ny']

			imgf = flipimage(options,img,ctf)	

			imgffileout = options.path + '/' + os.path.splitext(imgfilebase)[0]+'_flipped' + imgfilextension

			imgf['spt_phaseflipped']='entire image'
			imgf['spt_defocus_mean'] = ctf.defocus
			imgf['apix_x'] = ctf.apix
			imgf['apix_y'] = ctf.apix
			
			if options.verbose:
				print "\nwriting image out to %s, with size %d,%d, sigma %f, mean %f" %(imgffileout,imgf['nx'],imgf['ny'],imgf['sigma'],imgf['mean'])
			
			imgf.write_image(imgffileout,0)
			print "\nWrote flipped image to", imgffileout
		else:
			print "\nERROR: no ctf for image at index %d" %( indx )


	outflippedstack = ''
	extension = '.mrc'
	images2stack = options.path + '/*_flipped_w*' 
	if options.tiltseries:
		base=os.path.basename(options.tiltseries)
		extension = os.path.splitext(base)[1]
		outflippedstack = options.path + '/' + os.path.splitext(base)[0]+ '_flipped_w.hdf'
	elif options.imagestem:
		outflippedstack = options.path + '/' + options.imagestem + '_flipped_w.hdf'
		
	cmd = 'e2buildstacks.py --stackname ' + outflippedstack + ' ' + images2stack
	cmd += ' && e2fixheader.py --input ' + outflippedstack + ' --stem apix --stemval ' + str(apix) + ' --valtype float'
	
	outflippedstack3d = outflippedstack.replace('.hdf','.mrc')
	cmd += ' && e2proc2d.py ' + outflippedstack + ' ' + outflippedstack3d + ' --outmode int16 --twod2threed'
	#cmd += ' && e2fixheader.py --input ' + outflippedstack3d + ' --stem apix --stemval ' + str(apix) + ' --valtype float'
	if extension != '.mrc':
		finalname = outflippedstack3d.replace('.mrc',extension)
		os.rename( outflippedstack3d, finalname )
	
	findir = os.listdir( options.path )
	for f in findir:
		if 'sptctftmp' in f:
			os.remove( options.path + '/' + f)
	
	runcmd( options, cmd )
	return

	
def flipimage(options,img,ctf):	

	#prj=img.copy() 
	#print "In phase flipper, PRJ min, max, sigma, mean", prj['minimum'],prj['maximum'],prj['sigma'],prj['mean']
	
	nx=img['nx']
	ny=img['ny']
	#apix=img['apix_x']
	
	#print "\n(e2tomo_ctf)(flipimage) apix",ctf.apix
	#ctf.ac = 0.05
	#print "ctf.ac", ctf.ac
	
	if nx != ny:
		#img.process_inplace('filter.ctf',{"ctf":ctf, "type":1}) 
		print "image is NOT square nx, ny are ",nx,ny
		
		maxsize = max(nx,ny)
		img = clip2d( img, maxsize )
		#img_corrected = img.process('filter.ctfcorr.simple',{"ac":ctf.ac, "voltage":ctf.voltage,"apix":apix,"defocus":ctf.defocus}) 
	
	if options.fixctfhighpass:
		img_corrected = img.process('filter.ctfcorr.simple',{"voltage":ctf.voltage,"apix":ctf.apix,"defocus":ctf.defocus}) 

	img_fft = img.do_fft()
	flipim = img_fft.copy()	
	ctf.compute_2d_complex( flipim,Ctf.CtfType.CTF_SIGN )

	img_fft.mult( flipim )
	img_corrected = img_fft.do_ift()
	img_corrected['apix_x'] = ctf.apix
	img_corrected['apix_y'] = ctf.apix

	print "\nflipped image with this ctf", ctf
	print "ctf.defocus", ctf.defocus

	return img_corrected
	

def flipstripbystrip(options,imgsdict,ctfs,angles):
	
	for indx in imgsdict:
		print "\n(e2tomo_ctf)(flipstripbystrip) working on image",indx,imgsdict[indx]
	
		print "\nI will phaseflip these images computing the defocus gradient strip by strip",len(imgsdict)
		
		ctf = ctfs[indx]
		angle = angles[indx]
		
		if ctf:
			imgfile = imgsdict[indx][0]
			imgfilebase = os.path.basename(imgsdict[indx][0])
			imgfilextension = os.path.splitext(imgfilebase)[1]
			print "\nLoading image file", imgfile
			
			img = EMData( imgfile, 0, True )
			
			#hdr = EMData( imgfile, 0, True)
			nx = img['nx']
			ny = img['ny']
			apix = img['apix_x']
			
			maxsize = max(nx,ny)
			#padsize = maxsize + options.tilesize
			#img_padded = clip2d( img, padsize, padsize )
			
			nstripsx = int(nx/options.tilesize)
			nstripsy = int(ny/options.tilesize)

			if options.correctionwidth:
				nstripsx = int(nx/options.correctionwidth)

			
			finalimage = EMData( nx, ny )
			finalimage.to_zero()
			
			for i in range(nstripsx):

				startx = i*options.tilesize
				#startx = i*options.correctionwidth
				
				if options.correctionwidth:
					startx = i*options.correctionwidth

				xpixels = options.tilesize
				if options.correctionwidth:
					xpixels = options.correctionwidth

				
				#on the last strip, you have to go all the way to the end of the image; so this last strip might be larger than the others	
				if i == nstripsx-1:
					xpixels = nx-startx
				
				#stripcenterXpixels = int( ( startx + (startx + xpixels) )/ 2 )
				stripcenterXpixels = int( ( 2*startx + xpixels )/ 2 )

				stripcenterXmicrons = stripcenterXpixels * apix/10000 #this is in micrometers
				
				'''
				For positive tilt angles (counter clockwise) the defocus decreases (the particle is more overfocus, less defocused) for positions px right of the tilt axis
				while defocus increases for particles left of the tilt axis (they are more defocused).
				For negative tilt angles (clockwise) the defocuses increases (the particle is more defocused)for px right of the tilt axis while
				defocus decreases (more overfocused, less defocused) for particles left of the tilt axis.
				'''
				px = ( stripcenterXpixels - nx/2.0 ) * apix/10000	#this is the position of the strip center relative to the tilt axis
				if px < 0:
					print "\npx (in microns) is left of the tilt axis", px
				elif px > 0:
					print "\npx (in microns) is right of the tilt axis", px
				elif px==0:
					print "\npx (in microns) is on the tilt axis", px
			
				dzx = -1 * px * numpy.sin( math.radians( angle ) )		#the -1 accounts for the fact that positive tilt angles are clockwise, negative counter clockwise
		
				if angle < 0.0:
					print "\ngiven a negative, CLOCKWISE tilt angle=%f, and coordx=%f pixels, px=%f microns, THEN dzx=%f microns" %( angle,stripcenterXmicrons,px,dzx) 
				if angle > 0.0:
					print "\ngiven a positive, COUNTER CLOCKWISE tilt angle=%f, and coordx=%f pixels, px=%f microns, THEN dzx=%f microns" %( angle,stripcenterXmicrons,px,dzx) 
					
				newdefocus = ctf.defocus + dzx 
				print "\ntherefore, for angle=%f, and defocus=%f, the first corrected defocus is NEWdefocus1=%f" % ( angle, ctf.defocus, newdefocus )
				
				stripctf = EMAN2Ctf()
				stripctf.from_dict({ 'defocus':newdefocus, 'bfactor':ctf.bfactor, 'ampcont':ctf.ampcont, 'apix':ctf.apix, 'voltage':ctf.voltage, 'cs':ctf.cs })			
				
				pad = options.tilesize/2
				
				
				for j in range(nstripsy):
					starty = j*options.tilesize	
					ypixels = options.tilesize

				
					#on the last strip, you have to go all the way to the end of the image; so this last strip might be larger than the tilesize	
					
					if j == nstripsy-1:
						ypixels = ny-starty
				
					#r = Region(startx - pad, starty - pad, xpixels + options.tilesize, ypixels + options.tilesize)
					#r = Region(startx - pad, starty - pad, options.tilesize + pad, options.tilesize + pad)
					#r = Region(startx - pad, starty - pad, options.tilesize + pad, options.tilesize + pad)
					
					stripcenterYpixels = int( (2*starty + ypixels) / 2 )

					r = Region( stripcenterXpixels - options.tilesize, stripcenterYpixels - options.tilesize, 2*options.tilesize, 2*options.tilesize )
					

					print "\nclipping region", r
				
					clipr_padded = EMData()
					clipr_padded.read_image( imgfile, 0, False, r )
					
					#clipr_padded.write_image(options.path + '/stripspadded_' + str(indx).zfill(len( str(nstripsx) )) +'.hdf',i) 
				
					'''
					try:
						clipr['nx']
					except:
						print "\nERROR: clipr is of type", type(clipr)
						print "used region r=",r
						print "imgfile", imgfile
						try:
							hdr=EMData(imgfile, 0, True)
							print "with hdr", hdr
						except:
							print "invalid imgfile", imgfile
							sys.exit()
						sys.exit()
					'''
								
					#clipr_padded = clip2d( clipr, options.tilesize + 40, ny + 40 )
				
					clipr_flipped = flipimage( options, clipr_padded, stripctf )
					#clipr_flipped.write_image(options.path + '/stripsflipped_' + str(indx).zfill(len( str(nstripsx) )) +'.hdf',i) 

					#stripxsize = options.tilesize
					#if i == nstrips-1:
					#	stripxsize = xpixels
					
					clipr_flipped_originalsize = clip2d( clipr_flipped, xpixels, ypixels ) 
					#r = Region(startx - pad, starty - pad, options.tilesize + pad, options.tilesize + pad)

					#clipr_flipped_originalsize.write_image(options.path + '/stripsflippedcropped_' + str(indx).zfill(len( str(nstripsx) )) +'.hdf',i) 

					#if j != 0:
					#	clipr_flipped_originalsize.to_zero()
					#	print "\nMMMMMMADE clip 0"

					print "\nfor strip, column=%d, row=%d, clipping corrected strip to size nx=%d, nx=%d, ny=%d, ny=%d" %(i,j,xpixels,ypixels,clipr_flipped_originalsize['nx'],clipr_flipped_originalsize['ny'])
					print "and inserting it at x=%d, y=%d" %(startx,starty)
					print "mean=%.4f, sigma=%.4f" %( clipr_flipped_originalsize['mean_nonzero'], clipr_flipped_originalsize['sigma'] )


					#strip_final = clip2d( img_flipped_originalsize, options.tilesize, ny )
					#finalimage.insert_clip( strip_final, r )
				
				
					#SIMPLE ADDITION OF CORRECTED IMAGE			
					#edge2 = nx - startx - options.tilesize
					#strip_final = img_flipped_originalsize.process('mask.zeroedge2d',{'x0':startx,'x1':edge2})
					#finalimage += strip_final 

					finalimage.insert_clip( clipr_flipped_originalsize, (startx,starty) )
				
					#WEN JIANG
					#clip = getClip(image=img, xc=i, yc=ny/2, nx=boxsize, ny=ny, angle=0)
					#clip.process_inplace("mask.fillflat", {"random": 1})
					#clip.process_inplace("filter.ctf", {"ctf":ctf, "type":1})
					#clipcenter = clip.get_clip(EMAN2.Region(boxsize/2, 0, 1, ny))	# take the center line
					#fixed.insert_clip(clipcenter, (i, 0))
				

			imgffileout = options.path + '/' + os.path.splitext(imgfilebase)[0]+'_flipped_s' + imgfilextension
			#imgffileout = options.path + '/' + os.path.splitext(imgfilebase)[0]+'_flipped_s.hdf'

			finalimage['spt_phaseflipped']='strips'
			finalimage['spt_defocus_mean'] = ctf.defocus
			print "\nbefore normalization, finalimage mean=%.4f, sigma=%.4f" %(finalimage['mean_nonzero'],finalimage['sigma'])
			
			finalimage.process_inplace('normalize')

			
			print "\nwriting image out to %s, with size %d, %d, sigma %f, mean %f" %(imgffileout,finalimage['nx'],finalimage['ny'],finalimage['sigma'],finalimage['mean'])
			
			finalimage.write_image(imgffileout,0)
			print "\nwrote flipped image to", imgffileout


		else:
			print "\nERROR: no CTF for image at indx %d" %(indx)
			
			
	outflippedstack = ''
	extension = '.mrcs'
	images2stack = options.path + '/*_flipped_s*' 
	if options.tiltseries:
		base = os.path.basename(options.tiltseries)
		extension = os.path.splitext(base)[1]
		outflippedstack = options.path + '/' + os.path.splitext(base)[0]+ '_flipped_s.mrcs' #+ extension
	elif options.imagestem:
		outflippedstack = options.path + '/' + options.imagestem + '_flipped_s.mrcs'
		
	cmdstack = 'e2buildstacks.py --stackname ' + outflippedstack + ' ' + images2stack
	runcmd(options,cmdstack)

	cmdfixheader = 'e2fixheader.py --input ' + outflippedstack + ' --stem apix --stemval ' + str(apix) + ' --valtype float'
	runcmd(options,cmdfixheader)
	
	if extension != '.mrcs':
		os.rename(outflippedstack,outflippedstack.replace('.mrcs',extension))
		outflippedstack = outflippedstack.replace('.mrcs',extension)
	#outflippedstack3d = outflippedstack.replace('.hdf','.mrc')
	#cmdtype = 'e2proc2d.py ' + outflippedstack + ' ' + outflippedstack3d + ' --outmode int16 --twod2threed'	
	#runcmd( options, cmdtype )

	findir = os.listdir( options.path )
	for f in findir:
		if 'tomoctftmp' in f:
			os.remove( options.path + '/' + f )
		if 'sptctftmp' in f:
			os.remove( options.path + '/' + f)
	
	#testdisplay = EMData(outflippedstack,0)
	#display(testdisplay)
	
	return

	
def fullcorrection(options,img,ctf):
	
	prj=img.copy() 
	
	prj_fft = prj.do_fft()

	flipim = prj_fft.copy()	
	ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_TOTAL)

	prj_fft.mult(flipim)

	intermediate = prj_fft.copy()

	prj_flipped=prj_fft.do_ift()
	
	return prj_flipped, intermediate
		
	
def tilerfft(options, angle, imgf, currentstrip, nstrips, start, end):	
	kk = 0
	fftcumulative=None	
	nbx=0
	
	imgfhdr = EMData( imgf,0, True)
	ny = imgfhdr['ny']

	signtag = 'p'
	if angle < 0.0:
		signtag ='m'
		
	step=options.tilesize
	if options.overlaptiles:
		step = int(options.tilesize/2)
		print "\n(e2spt_ctf)(tilerfft) --overlaptiles is on, therefore step is", step
	
	#for x in range( micrographstarts[m], micrographstarts[m] + micrographwidth - options.tilesize + 1, options.stripfitstep ):
	
	tilesgood=0
	tilesbad=0
	
	for x in range( start, end, step ): 
		#print "\nsumming over tiles along y"
		for y in range(0, ny - options.tilesize+1, options.tilesize):
			#print "tile at y", y
			#clipr = imgt.get_clip(Region(x,y, options.tilesize, options.tilesize))

			r=Region(x,y, options.tilesize, options.tilesize)
			clipr = EMData()
			clipr.read_image(imgf,0,False,r)

			
			
			if clipr['sigma']:
				
				allgood = 1
				
				if float(options.prunetest) >= 0.0:
					allgood = checkcorners( clipr, options )
				
				if allgood:			
				
					#fftcumulativeimgfile = options.path  + '/angle_' + signtag + str( int(math.fabs( round(angle) ) + '_strip' + str(currentstrip).zfill(len(str(nstrips))) + '_fft.hdf'

								
					clipr.process_inplace("normalize.edgemean")
					
					if options.savestriptiles:
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
					tilesgood+=1
				else:
					tilesbad+=1
					pass
					#"WARNING: tile excluded because it has at least one bad corner!"
				
			else:
				tilesbad+=1
				pass
				#print "WARNING: tile excluded because sigma is zero!"

	if nbx > options.mintiles:
	
		if fftcumulative:
			fftcumulative.mult(1.0/(nbx*options.tilesize**2))
			fftcumulative.process_inplace("math.sqrt")
			fftcumulative["is_intensity"]=0				# These 2 steps are done so the 2-D display of the FFT looks better. Things would still work properly in 1-D without it
	
			signtag = 'p'
			if angle < 0.0:
				signtag ='m'
	
			#plotname = options.path + '/fit_' + str( imgindx ).zfill( len(str( nangles ))) + '_' + signtag + str( int(math.fabs( round(angle) ) )).zfill(3) +'.png'
	
			if options.saveffts:
				fftcumulativeimgfile = options.path  + '/angle_' + signtag + str( int(math.fabs( round(angle) ))) + '_strip' + str(currentstrip).zfill(len(str(nstrips))) + '_fft.hdf'
				fftcumulative.write_image(fftcumulativeimgfile,0)	
	
			return fftcumulative, tilesgood
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
	

def fitdefocus( ffta, angle, apix, options, nsubmicros, currentsubmicro, defocusmin, defocusmax, defocusstep, x=0, imgindx=0 ):

	ctfi = EMAN2Ctf()
	#ds=0.0
	ds = 1.0/( options.apix * options.tilesize )

	bg_1d=[]
	stripdefocus = 0.0
	
	if imgindx:
		f=open(options.path+'/ds_new.txt','a')
		line=str(imgindx).zfill(2)+'\t'+str(ds)+'\n'
		f.write(line)
		f.close()
	
	fftbg = ffta.process("math.nonconvex")
	fft1d = ffta.calc_radial_dist(ffta.get_ysize()/2,0.0,1.0,1)	# note that this handles the ri2inten averages properly
		
	signtag = 'p'
	if angle < 0.0:
		signtag ='m'
	
	try:
		# Compute 1-D curve and background
		
		
		
		bg_1d = e2ctf.low_bg_curve(fft1d,ds)
		
		#initial fit, background adjustment, refine fit, final background adjustment
		#ctf = e2ctf.ctf_fit(fft1d,bg_1d,bg_1d,ffta,fftbg,options.voltage,options.cs,options.ampcont,apix,1,dfhint=( defocusmingradient, defocusmaxgradient, step ))
	

		ctfi = e2ctf.ctf_fit( fft1d, bg_1d, bg_1d, ffta, fftbg, options.voltage, options.cs, options.ampcont, apix, 1,dfhint=( defocusmin, defocusmax, defocusstep ) )
		bgAdj(ctfi,fft1d)
	except:
		print "\nctf fit failed! first try"
		print "len fft1d is", len(fft1d)
		print "ffta is", ffta
		#ctfi = None
			
	try:
		#ctf = e2ctf.ctf_fit(fft1d,ctf.background,ctf.background,ffta,fftbg,options.voltage,options.cs,options.ampcont,apix,1,dfhint=( defocusmingradient, defocusmaxgradient, step))
		ctfi = e2ctf.ctf_fit(fft1d,ctfi.background,ctfi.background,ffta,fftbg,options.voltage,options.cs,options.ampcont,apix,1,dfhint=( defocusmin, defocusmax, defocusstep ))	
		bgAdj(ctfi,fft1d)

		#if options.astigmatism: 
		#	e2ctf.ctf_fit_stig(ffta,fftbg,ctf)
	except:
		print "ctf fit failed! second adjustment try"
		print "len fft1d is", len(fft1d)
		print "ffta is", ffta
		#ctfi = None
		

	
	ctfi.background = bg_1d
	ctfi.dsbg = ds
	ctf1d = ctfi.compute_1d
	stripdefocus = ctfi.defocus
	
	return stripdefocus


'''
function to calculate the global defocus of an image (tiling the entire image and computing 
the incoherent average and FFT using all the tiles)
'''
def calcglobaldefocus(options, apix, imagefilenames, angles):
	
	print "\n(e2spt_tomo)(calcglobaldefocus)"
	
	defocuswiggle = None
	globaldefocus = None
	globalmiddle = None
	centerdefocus = None
	defocuserror=0
	
	defocusmean = 0.0
	
	defocusmin = 0.75
	defocusmax = 10
	defocusstep = 0.1
	
	'''
	read image size from the header of the first image
	'''
	nx=EMData(imagefilenames[0][0],0,True)['nx']
	endstart = nx - options.tilesize + 1
	
	'''
	constrain defocus search around --targetdefocus if provided; otherwise, calculate
	the global defocus of the lowest tilt image first, as this is the most likely to 
	succeed, and use this to constrain defocus fitting
	'''
	if options.targetdefocus:
		defocusmin = options.targetdefocus - 3.0
		defocusmax = options.targetdefocus + 3.0
	else:
		
		lowesttiltangleindx = len(angles)/2.0 	#c:by default, the lowest tilt angle should be in the middle of the tiltseries
		
		lowesttiltangle = min([ math.fabs(ang) for ang in angles.values() ])	#c:find the actual smallest tiltangle, then its index
		try:
			lowesttiltangleindx = angles.keys()[angles.values().index(lowesttiltangle)]	#c:if the lowest tiltangle is actually positive
		except:
			lowesttiltangle=-1*lowesttiltangle
			try:
				lowesttiltangleindx = angles.keys()[angles.values().index(lowesttiltangle)]	#c:if the lowest tiltangle is negative
			except:
				print "(e2tomo_ctf)(calcglobaldefocus) ERROR: determined lowest tilt angle %f somehow is not in --angles or --anglesfile. Using the middle index by default" %(lowesttiltangle)
				#sys.exit()
		lowesttiltimgfile = imagefilenames[lowesttiltangleindx][0]
		#lowesttiltimg = EMData(lowesttiltimgfile,0)
		
		fftlowesttilt, ntiles = tilerfft( options, lowesttiltangle, lowesttiltimgfile, 0, 1, 0, endstart )
		gdeflowesttilt = fitdefocus( fftlowesttilt, lowesttiltangle, apix, options, 1, 0, defocusmin, defocusmax, defocusstep, 0, str(lowesttiltangleindx))
		
		if gdeflowesttilt > 0.75 and gdeflowesttilt < 10:
			defocusmin = gdeflowesttilt - 2.0
			defocusmax = gdeflowesttilt + 2.0
		
	if defocusmin < 0.0:
		defocusmin = 0.75
	
	'''
	iterate trough images to calculate the global defocus of each
	'''
	allglobaldefocuses = {}
	for imgindx in imagefilenames:
		'''
		First fit by tiling the ENTIRE micrograph to get the average defocus without gradient compensation,
		just to get a first 'aproximate' estimate for the range of defocuses to look for
		'''
		angle = imagefilenames[imgindx][1]
		
		imgfile = imagefilenames[imgindx][0]
		#img = EMData(imgfile,0)
		
		imgfile = imagefilenames[imgindx][0]
		fftg, ntiles = tilerfft( options, angle, imgfile, 0, 1, 0, endstart )
		fftg.write_image(options.path +'/fftg_new' + str(imgindx).zfill(2) +'_'+ str(angle)+ '.hdf',0)
		globaldefocus = fitdefocus( fftg, angle, apix, options, 1, 0, defocusmin, defocusmax, defocusstep, 0, str(imgindx))
		
		if globaldefocus:
			allglobaldefocuses.update({ imgindx:globaldefocus })
			
			print "\n(e2spt_ctf)(calcglobaldefocus) global defocus is", globaldefocus
			
		else:
			print "\n(e2spt_ctf)(calcglobaldefocus) WARNING: global defocus measurement failed for image index %d, image file %s" %(imgindx, imgfile)
			if options.targetdefocus:
				allglobaldefocuses.update({ imgindx:options.targetdefocus })
					
		
		
		#print "bad version: global defocus middle is", globalmiddle
		#sys.exit()
		
		"""
		#get central region of 3 tile width
		r = Region( img['nx']/2 - options.tilesize - int(ceil(options.tilesize/2.0)), 0, 3*options.tilesize, img['ny'])
		globalmiddle = ( nx/2.0 ) * apix / 10000	#convert to micrometers, so x and y axis are in the same units
		#r = Region( img['nx']/2 - options.tilesize, 0, 2*options.tilesize, img['ny'])
	
		
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
		"""
	
	print "\n\n(e2spt_ctf)(calcglobaldefocus) allglobaldefocuses",allglobaldefocuses
	
	
	globaldefocuseslist = []
	defocusstd=0.0
	defocusmean=0.0
	
	if allglobaldefocuses:
		for indx in allglobaldefocuses:
			globaldefocuseslist.append(allglobaldefocuses[indx])
		
		defocusstd = numpy.std(globaldefocuseslist)
		defocusmean = numpy.mean(globaldefocuseslist)
		
		for indx in allglobaldefocuses:
			d = allglobaldefocuses[indx]
			if math.fabs(d-defocusmean) > 2*defocusstd:
				imgfile = imagefilenames[indx][0]
				img = EMData(imgfile,0)
				angle = angles[indx]
				
				print "\n(e2spt_ctf)(calcglobaldefocus) RETRYING global defocus for image %s, indx %d, because its defocus was %f, while the mean defocus was %f with an std of %f" %( imgfile, indx, d, defocusmean, defocusstd )
				
				defocusmin = defocusmean - 2*defocusstd
				defocusmax = defocusmean + 2*defocusstd
				
				imgfile = imagefilenames[indx][0]
				
				fftg, ntiles = tilerfft( options, angle, imgfile, 0, 1, 0, endstart )
				#fftg.write_image(options.path +'/fftg_new' + str(imgindx).zfill(2) +'_'+ str(angle)+ '.hdf',0)
				newglobaldefocus = fitdefocus( fftg, angle, apix, options, 1, 0, defocusmin, defocusmax, defocusstep, 0, str(indx))
				
				allglobaldefocuses[indx]=newglobaldefocus
				
		
	if allglobaldefocuses:	
		angleslist = []
		lines=[]
		globaldefocuseslist = []
		angledefocuspairs={}
	
		for indx in angles:
			if indx in allglobaldefocuses:
				angle = angles[indx]
				defocus = allglobaldefocuses[indx]
		
				globaldefocuseslist.append( defocus )
				angleslist.append( angle )
				
				angledefocuspairs.update({indx:[angle,defocus]})
				
			else:
				print "\n(e2spt_ctf)(calcglobaldefocus) WARNING: missing blobal defocus for image at angle %f" %( angle )
				globaldefocuseslist.append( -1.0 )
				angleslist.append( angle )
			
		lines=[]
		linesplot=[]
		for i in angledefocuspairs:
			angle=angledefocuspairs[i][0]
			defocus=angledefocuspairs[i][1]
			line = str(i)+'\t'+str(defocus)+'\n'
			lines.append(line)
			lineplot = str(angle)+'\t'+str(defocus)+'\n'
			linesplot.append(lineplot)
		
		g=open( options.path + '/defocuses_global.txt','w' )
		g.writelines(lines)
		g.close()
		
		defocusmean = numpy.mean(globaldefocuseslist)
		defocusstd = numpy.std(globaldefocuseslist)
		
		h = open( options.path + '/eucentricity_variation_global.txt', 'w' )
		h.writelines([ 'spread='+str( max(globaldefocuseslist) - min(globaldefocuseslist) ) +'\n', 'mean='+ str(defocusmean) +'\n', 'std='+str(defocusstd)+'\n' ] )
		h.close()
		
		f = open( options.path + '/tiltangle_vs_defocus _global.txt', 'w' )
		f.writelines(linesplot)
		f.close()
		
		ylabel = 'Global defocus (micrometers)'
		xlabel = 'Tilt angle (degrees)'
		plotname = options.path + '/tiltangle_vs_defocus_global.png'
		title = 'Tiltangle vs global defocus'
		
		generalplotter( options, angleslist, globaldefocuseslist, xlabel, ylabel, plotname, title )
		
		return allglobaldefocuses, defocusmean, defocusstd
	else:
		print "\n(e2spt_ctf)(calcglobaldefocus) WARNING! All global defocuses estimations failed!"
	
	
def fitdefocusgradient( options, apix, imagefilenames, angles, globaldefocuses, defocusmean, defocusstd ):
		
	finaldefoci={}
	#print "\n(e2spt_ctf)(fitdefocusgradient) imagefilenames", imagefilenames
	
	#angles.sort()	#use to get image index

	angerrors = {}
	#angerrorlines = []
	imgnum = 0
	
	maxangle = max( [ math.fabs( a ) for a in angles] )
	
	skip = []
	if options.skipimgstrips:
		skip = [ int(i) for i in options.skipimgstrips.split(',') ]
	
	#anglestoexclude=[]
	#print "(e2spt_ctf)(sptctffit) len imagefilenames", len (imagefilenames)
	#sys.exit()
	
	nx=EMData(imagefilenames[0][0],0,True)['nx']
	globalmiddle = ( nx/2.0 ) * apix / 10000	#convert to micrometers, so x and y axis are in the same units
	
	'''
	calculate the maximum allowable depth of focus or variation in defocus across the specimen to achieve 2/3 Nyquist resolution
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
	allowabledzmicrons = allowabledz/10000
	allowabledzpixels = allowabledz/apix

	print "\n(e2spt_ctf)(fitdefocusgradient) the theoretical allowable depth of focus or defocus variation limit in angstroms is %.2f to reach 2/3 nyquist resolution %.2f" %( allowabledz, twothirdsnyquist )

	if options.depthfocus:
		allowabledz = options.depthfocus
		print "however, the allowable depth of focus in micrometers has been set via --depthfocus and is", allowabledz
		#print "which in angstroms is", allowabledz*10000
		allowabledzmicrons = allowabledz/10000
		allowabledzpixels = allowabledz/apix
	
	'''
	iterate through images to fit the defocus gradient of each
	'''
	
	
	defocusstep = 0.1
	
	stripfitstep = options.tilesize
	
	for imgindx in imagefilenames:
	
		globaldefocus = globaldefocuses[imgindx]
		defocusmin = globaldefocus-2.0
		defocusmax = globaldefocus+2.0
		
		angle = imagefilenames[imgindx][1]

		print "\n(e2spt_ctf)(fitdefocusgradient) autofitting defocus gradient for imgindx %d, image %s, angle %f, apix %f, progress = %d/%d" %( imgindx, imagefilenames[imgindx][0],angle,apix,imgindx,len(angles))
		
		imgfile = imagefilenames[imgindx][0]
		imghdr = EMData( imgfile, 0, True )
		
		nx = int( imghdr['nx'] )
		ny = int( imghdr['ny'] )
		
		xs = []
		imgdefocuses = []
		imagemiddle = imghdr['nx']/2
		m=0
		b=0
		defocuscalc = 0	
		nxMicrometers = nx*apix/10000.00
		imgdefocuses = []
		faileddefs = []
		failedmids = []
		#do not strip-fit if an image is bad
		
		if imgindx not in skip:
			'''#
			#Position of the tilt axis in micrometers
			#'''
			pxta = ( nx/2.0 ) * apix/10000
		
			'''#
			#Find the distance dx100dz away from the tilt axis, across which the vertical distance 
			dzx changes by 100 nm, i.e., 0.1 micrometers, given the specified ice thickness
			#'''
			
			stripwidth = nx 						#at zero degrees dzx100 is the micrograph size, since there should be no defocus variation due to tilt
				
			if math.fabs( angle ) > 0.5:
				
				stripwidthold = math.fabs( 0.1 / numpy.sin( math.radians( angle ) ) * 10000/apix )
				print "\nOOOOOOOO for angle %f old strip width in pixels with no thickness was %f" %(angle, stripwidthold)
				
				if options.thickness:
					depthdefocus = math.fabs( options.thickness / numpy.cos( math.radians( angle ) ) )	#icethickness, in pixels, will contribute to defocus variation with tilt, in addition to the tilting itself
				
					stripwidthold += depthdefocus				
				
				stripwidthold*=2 #this factor of two accounts for the fact that the old method only calculated things for half of the micrograph (from the tilt axis)

				
				print "old depthfocus was", math.fabs( options.thickness / numpy.cos( math.radians( angle ) ) )
				print "using this thickness", options.thickness
				print "old strip width in pixels was", stripwidthold
				
				
				#cos/sin is cot, but the function does not exist in numpy
				stripwidth = allowabledzpixels * (numpy.cos( math.radians (angle) ) / numpy.sin( math.radians (angle) ) )  #this is in pixels
				
				#print "however, new strip width is", stripwidth
				depthfocusnew = -1* options.thickness / numpy.sin( math.radians( angle ))
				
				stripwidth += depthfocusnew
				stripwidth = int( math.fabs( stripwidth ))
				
				print "however, new strip width accounting for depth of focus, it is",stripwidth
				
				if stripwidth > nx:
					stripwidth = nx
		
			'''#
			#submicrograph width (the region at "equal/constant defocus" to tile) will be twice dx100dz  
			#'''
			#micrographwidth = int( dx100dz * 2 )
			micrographwidth = int( stripwidth )
		
		
		
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
		
			#print "nmicrosint is %d since (nmicrosint+2)*micrographwidth is %d while nx is %d" %( nmicrosint, (nmicrosint+2)*micrographwidth, nx )
		
			excedent = nx - nmicrosint * micrographwidth
		
			#print "excendent is", excedent
		
			#excedent = nmicros - nmicrosint
		
			#nmicrosroundedup = math.ceil( nmicros )
		
			micrographstarts = []
		
			aux = 0
			if excedent:
				if excedent >= options.tilesize*2 and not options.excludeedges:
					aux = 2
		
		
			adjuststart = 0								#Flag parameter used (below) for highly tilted images in which the region of "constant defocus" is smaller than the tilesize
			if micrographwidth < options.tilesize:
				stripfitstep = micrographwidth
				micrographwidth = options.tilesize
				adjuststart = 1
			
			print "--stripfitstep is!!!!!!!!!!!!!!!!!!!", stripfitstep
		
			include = 1
			if not adjuststart:
				for iz in range( nmicrosint + aux ):		#You'll usually have two extra micrographs to include the excedent at either side of the central region where you can fit an exact odd multiple of submicrograph widths
					#print "iz is", iz
				
					if aux:
				
						if iz == 0:							#first submicrograph starting point
							start = 0
						elif iz == nmicrosint + 1:				#last submicrograph starting point; range goes from 0 to microsint + 2 at most, without including the upper bound
							start = nx - micrographwidth
						else:
							start = (iz-1)*micrographwidth + excedent/2
						
					else:
						start = (iz)*micrographwidth + excedent/2		#Don't bother with edge micrographs if aux is 0
					
					#print "withOUT adjuststart, start to append is", start
				
				
					#if include:
				
					micrographstarts.append( int(start) )
				
			elif adjuststart:
				nmicrosint = int( math.floor( ( nx - options.tilesize ) / stripfitstep ) )
			
				excedentnew = nx - options.tilesize * nmicrosint
			
				init = 0
				if options.excludeedges:
					init = 1
			
				for ii in xrange(init,nmicrosint):
	
					start = ii*stripfitstep
				
					micrographstarts.append( int(start) )
			
				if excedentnew and not options.excludeedges: 	#plus add the final one if there's excedentnew
					start = nx - options.tilesize
					micrographstarts.append( start )
				
			micrographstarts.sort()
			#print "\nfor img %d micrographstarts are" % (imgindx) 
			print micrographstarts	
		
			#micromids = [h+micrographwidth/2 for h in micrographstarts]
			micromids = []
		
			for m in range( len (micrographstarts)):
				'''#
				#The defocus is fitted per submicrograph (regions of pseudo-constant defocus); therefore, things are reset here
				#'''
				
				endstartm = micrographstarts[m] + micrographwidth - options.tilesize + 1
				
				fftc, ntiles = tilerfft( options, angle, imgfile, m, len(micrographstarts), micrographstarts[m], endstartm )
				
				micrographmiddle = float( list(micrographstarts)[m] + micrographwidth/2.0 )
				micrographmiddlemicrometers =  ( list(micrographstarts)[m] + micrographwidth/2.0 ) * apix / 10000.00 #xaxis in micrometers too
				
				if options.verbose > 8:
					print "\n(e2tomo_ctf)(fitgradient) micrograph start is %f, micrograph width is %f, therefore micrograph middle is %f" %(micrographstarts[m],micrographwidth,micrographmiddle)
							
				if options.predictdefocus:
					dx = imagemiddle-micrographmiddle
					shift = dx * math.tan( math.radians(angle) ) * apix /10000
				
					#print "\nSHHHHHIFT is %f for angle %f at middle %f with imagemiddle %f and dx %f" %(shift,angle,micrographmiddle,imagemiddle,dx)
						
					predicteddefocus = globaldefocus+shift
					defocusmin = predicteddefocus-1.0
					defocusmax = predicteddefocus+1.0
						
				stripdefocus = None
				if fftc:
					
					if options.thickness and options.coords:
						icethicknessm = options.thickness * imghdr['apix_x']/10000
					
						defocuswiggley = math.fabs( icethicknessm / math.cos( math.radians(angle) ))
						if defocuserror:
							defocuswiggley += defocuserror
						
						defocuswigglex = 2*math.fabs( -1*(micrographmiddlemicrometers - imghdr['nx']/2.0)*math.sin( math.radians(angle) ) * imghdr['apix_x']/10000)
						
						defocuswiggle = defocuswigglex + defocuswiggley
						
						print "(e2spt_ctf)(fitdefocusgradient) icethicknessm is", icethicknessm
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
						print "WARNING! bad strip for image at index %d; defocus for strip at x %d is None, and the number of good tiles was %d" %( imgindx, micrographstarts[m], ntiles )

					
			
			
				#print "micrographmiddlein micrometers is", micrographmiddle
			
				if stripdefocus:
					#imgdefocuses.append( stripdefocus*10000/apix )		#defocus in pixels
					imgdefocuses.append( stripdefocus )					#defocus in micrometers
				
					micromids.append( micrographmiddlemicrometers )	
					print "\nappended (good) micrographmiddle is", micrographmiddlemicrometers	
				else:
					print "\nappending to failed results, for which the number of good tiles ntiles was", ntiles
					faileddefs.append( (defocusmin+defocusmax)/2 )
					failedmids.append( micrographmiddlemicrometers )
			
			#xs = numpy.array( [i*options.stripfitstep + options.tilesize/2.0 for i in range(len(imgdefocuses))] )
		
			print "micromids are", micromids
			#xs =numpy.array( micromids )
			xs=micromids
			#imgdefocuses = numpy.array( imgdefocuses )
		
			#print 'xs are', xs, type(xs)
		
			#print "\bPLOTTING\n\n"
		
			
		else:
			print "img being skipped for strip-based fitting"
			xs.append(imagemiddle)
			imgdefocuses.append(globaldefocus)
			#anglestoexclude.append(angle)
		
		xs = numpy.array( xs )
		imgdefocuses = numpy.array( imgdefocuses )
		
		
		if xs.any() and imgdefocuses.any():
		
			for y in range( len(imgdefocuses)):
				imgdefocuses[y] = imgdefocuses[y] *-1
			

			slope, b = numpy.polyfit(xs, imgdefocuses, 1)
			defocuscalc = slope * nxMicrometers/2.0 + b
			
			#if m < 0:
			#	if angle > 0:
			#		m *= -1
			#if m > 0:
			#	if angle < 0:
			#		m *= -1
			
			
			anglecalc = math.degrees(numpy.arctan( slope ))
			
			if len(xs) < 2 and len (imgdefocuses) < 2:
				anglecalc = angle
				print "defocus calc is derived from middle strip only! it would have been", defocuscalc
				defocuscalc = imgdefocuses[0]
				print "but is", defocuscalc
				
				
			#if angle not in anglestoexclude:
			angerror = math.fabs( angle - anglecalc )
			
			print "angle, anglecalc, and angerror are", angle, anglecalc, angerror
			angerrors.update( { angle:angerror } )
			
			if angerror > 15.0:
				defocuscalc = globaldefocus
				#middef = imgdefocuses[0]
				#if len(imgdefocuses) > 2:
				#	middef = imgdefocuses[len(imgdefocuses)/2]
				#	defocuscalc = (globaldefocus+middef)/2.0
					
				print "\nWARNING: ERRRRRRRRRRROR; angerror %f > 15.0; using globaldefocus %f, defocuscalc %f" %(angerror, globaldefocus, defocuscalc)
										
			angleindx = angles.keys()[angles.values().index(angle)]

			sptctfplotter( options, nxMicrometers, xs, imgdefocuses, maxangle, angle, angleindx, len(angles), imgindx, slope, b, globaldefocus, globalmiddle, faileddefs, failedmids )		

			#params = {'ampcont':options.ampcont,'apix':apix,'bfactor':options.bfactor,'cs':options.cs,'defocus':math.fabs(defocuscalc),'voltage':options.voltage}	
			#ctf = EMAN2Ctf()
			#ctf.from_dict(params)
			#ctfs.update( {angle:ctf} )
			
			#######ctf.from_dict({'defocus':params['defocus'],'bfactor':params['bfactor'],'ampcont':params['ampcont'],'apix':params['apix'],'voltage':params['voltage'],'cs':params['cs']})	
		
			finaldefoci.update( {imgindx:math.fabs(defocuscalc)} )
			
			
		else:
			print "\nWarning: All defocuses failed for this submicrograph. Nothing to plot."
		
		#if xs.any() and imgdefocuses.any():
			#pass nx in micrometers
		
		
				
		#imgnum += 1
		
	
	#angles.sort()
	
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

	lines=[]
	linesplot=[]
	finaldefocilist=[]
	angleslist=[]
	for i in finaldefoci:
		angle=angles[i]
		angleslist.append(angle)
		
		defocus=finaldefoci[i]
		finaldefocilist.append(defocus)
		
		line = str(i)+'\t'+str(defocus)+'\n'
		lines.append(line)
		lineplot = str(angle)+'\t'+str(defocus)+'\n'
		linesplot.append(lineplot)
	
	g=open( options.path + '/defocuses_fit.txt','w' )
	g.writelines(lines)
	g.close()
	
	defocusmean = numpy.mean(finaldefocilist)
	defocusstd = numpy.std(finaldefocilist)
	
	h = open( options.path + '/eucentricity_variation_fit.txt', 'w' )
	h.writelines([ 'spread='+str( max(finaldefocilist) - min(finaldefocilist) ) +'\n', 'mean='+ str(defocusmean) +'\n', 'std='+str(defocusstd)+'\n' ] )
	h.close()
	
	f = open( options.path + '/tiltangle_vs_defocus_fit.txt', 'w' )
	f.writelines(linesplot)
	f.close()
	
	ylabel = 'Fit defocus (micrometers)'
	xlabel = 'Tilt angle (degrees)'
	plotname = options.path + '/tiltangle_vs_defocus_fit.png'
	title = 'Tiltangle vs fit defocus'
	
	generalplotter( options, angleslist, finaldefocilist, xlabel, ylabel, plotname, title )
	
	return finaldefoci
	

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
	
	if failedxs and failedys:		
		for y in range( len(failedys)):
			failedys[y] = failedys[y] *-1

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
	miny = min(ydata)
	maxy = max(ydata)
	
	if gdefocus:
		if float(gdefocus) > float( maxy ):
			maxy = gdefocus
		if float(gdefocus) < float( miny ):
			miny = gdefocus
		
	if failedys:
		minyf = min( failedys )
		if float(minyf) < float(miny):
			miny = minyf
		maxyf = max( failedys )
		if float(maxyf) > float( maxy ):
			maxy = maxyf
	
	miny -= 1
	maxy += 1
	pylab.ylim([ miny, maxy ])

	#print "max x is", max(xaxis)

	pylab.xlim([ 0, nx + 0.5 ])

	ax.set_xlabel('X coordinate (micrometers)', fontsize=18, fontweight='bold')
	ax.set_ylabel('Defocus (micrometers)', fontsize=18, fontweight='bold')
	
	title ="Tilt image " + str( angleindx ).zfill( len(str( nangles ))) + ", angle=" + str(angle)
	pylab.title( title, fontweight='bold', fontsize=18 )
	
	#plt.plot(xdata, ydata, '.', markersize=10)
	plt.scatter(xdata, ydata, marker='.', s=200,alpha=0.9, color='k')
	
	if failedxs and failedys:
		
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
		
			
	print "\nmax y is", maxy
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
	
	print 'yrange', ylim1,ylim2
	
	print "\nmax x is", max(xaxis)
	
	pylab.xlim([xlim1, xlim2])
	
	print 'xrange', xlim1,xlim2
	
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
	return


def clip2d( img, nx, ny=0, imgycc=None, imgxcc=None ):
	
	sizex = nx
	sizey = nx
	if ny:
		sizey = ny
		
	imgxc = img['nx']/2
	if imgxcc != None:
		imgxc=imgxcc

	imgyc = img['ny']/2
	if imgycc != None:
		imgyc=imgycc

	Rimg =  Region( (2*imgxc - sizex)/2, (2*imgyc - sizey)/2, sizex , sizey )
	img.clip_inplace( Rimg )
	
	return img
	

"""
'''
================================
'''

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

def correctsubtilt( options, subtilts, angles, ctfs, apix, nangles, nimgs, framexsize, maxz, minz ):
	
	if options.reconstructor and options.reconstructor != 'None' and options.reconstructor != 'none': 
		options.reconstructor=parsemodopt(options.reconstructor)

	ii = 0
	globalAvgDefErrors=[]
	
	zfillfactor = len( str( len(subtilts) ) )
	
	#indxstoexclude=[int(i) for i in options.exclude.split(',')]
	
	#ss=0
	for sts in subtilts:
		imghdr = EMData( sts, 0, True )
		
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
			print '''WARNING: The number of angles %d does not coincide with number of images %d. 
			However, the actual angles being used should be obtained directly from the 2d image's header''' % ( nangles,  n )
			#sys.exit(1)
		
		flippedsts = options.path + '/' + os.path.basename( sts ).replace('.hdf','_phflip.hdf')
		#if options.outputstem:
		#	flippedsts = options.path + '/' + options.outputstem + 'ptcl' + str(ii).zfill( zfillfactor ) + '_phflip.hdf'
		
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
		
		#if options.save3d:
		#sts3d = options.path + '/' + os.path.basename( sts ).replace('.hdf','_PHFLIP3D.hdf')
		stack3d = options.path + '/stack_phflip3d.hdf'
		#if options.outputstem:
		#	stack3d = options.path + '/' + options.outputstem +'_phflip3d.hdf'
		
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
"""

if '__main__' == __name__:
	main()
