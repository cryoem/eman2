#!/usr/bin/env python
import os
import EMAN2
from EMAN2 import EMData
#from EMAN2 import EMUtil, EMArgumentParser, EMANVERSION
from sp_applications import header, project3d, mref_ali2d
from sp_utilities import get_im, write_header, model_circle, read_text_row, set_params2D, write_text_row, set_params_proj, compose_transform3, model_blank
from sp_statistics import ccc
from sp_fundamentals import rops_table, fft, rot_shift2D
from sp_projection import prep_vol, prgl
from math import sqrt, degrees, atan2
from sp_filter import filt_table
import sp_global_def
from datetime import datetime
from sp_logger import Logger, BaseLogger_Files, BaseLogger_Print
import shutil
import numpy as np
from sp_alignment import Numrinit, ringwe, search_range, align2d, align2d_scf, align2d_direct3
from sp_global_def import Util, sxprint
import glob
from sp_pixel_error import angle_diff
import argparse

# Set default values (global variables written in ALL CAPS)
TIMESTAMP_LENGTH = 23  # chars

USAGE = """
PURPOSE:
Compare projections of a 3D model to a stack of 2D projections.

To use angles stored in the image header:
%s <input_imgs> <input_volume>
General, required command-line parameters:
  1. Input image stack
  2. Input volume to be projected
Outputs:
  docangles.txt : Projection angles applied to input model
  comp-proj-reproj.hdf : Stack of reprojections (numbered 0,2,4...) and averages (numbered 1,3,5...)
  
General options:
%s <input_imgs> <input_volume> <output_directory> --classangles <angles_file> --select <img_selection_file> --prjmethod <interpolation_method> --verbose --display
Parameters:
  --prjmethod : Interpolation method : trilinear (default), gridding, nn (nearest neighbor)
  --verbose : Increase verbosity
  --display : Automatically open montage in e2display
   
To use angles from a VIPER or RVIPER run:
%s <input_imgs> <input_volume> <output_directory> --mode viper --classangles <angles_file> --classselect <img_selection_file>
Parameters:
  --classangles : Projection angles (For RVIPER, this file will have a name like main003/run002/rotated_reduced_params.txt)
  --classselect : Image selection file (For RVIPER, not all images will necessarily be assigned alignment angles, in which case the number of angles in the doc file above will differ from the number of images in the input stack. This file will have a name like main003/index_keep_image.txt)
 
To perform an internal round of projection matching against the input model:
%s <input_imgs> <input_volume> <output_directory> --mode projmatch --delta <angular_increment> --matchshift <shift_range> --matchrad <outer_radius> --matchstep <ring_step> --symmetry <optional_symmetry>
Parameters:
  --delta : Angular-sampling increment
  --matchshift : Maximum shift to allow during translation alignment (default 0)
  --matchrad : Outer alignment radius (default Automatic)
  --matchstep : Alignment radius step size (default 1)
  --symmetry : To limit angular projections (default c1)

To use the average projection angles from MERIDIEN refinement:
%s <input_imgs> <input_volume> <output_directory> --mode meridien --partangles <refinement_params> --partselect <substack_select> --refineshift <shift_range> --outliers <max_angle_diff> --refinerad <outer_radius> --refinestep <ring_step> --align <2d_alignment_method>
Parameters:
  --partangles : Input refinement parameters
  --partselect : Input substack selection file if particles removed before refinement (e.g., Substack/isac_substack_particle_id_list.txt)
  --outliers : Particles differing from average Euler angle by more than this threshold (in degrees) will be excluded from average calculation (default keep all)
  --refinerad : Outer alignment radius (default Automatic)
  --refineshift : Maximum shift to allow during translation alignment (default 0)
  --refinestep : Alignment radius step size (default 1)
  --align : Alignment method: apsh (default) or scf

Modified 2019-03-23

""" % ((__file__,)*5)
	
def parse_command_line():
	"""
	Parse the command line.  Adapted from sxmask.py

	Arguments:
	None:

	Returns:
	Parsed arguments object
	"""

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument(
		'classavgs', 
		type=str, 
		help='Input class averages')
	
	parser.add_argument(
		'vol3d', 
		type=str, 
		help='Input 3D reconstruction')
	
	parser.add_argument(
		'outdir', 
		type=str, 
		default='.', 
		help='Output directory (default .)')
	
	parser.add_argument(
		'--mode', 
		type=str, 
		default='viper', 
		help='Options: viper (default), projmatch, meridien')
	
	parser.add_argument(
		'--prjmethod', 
		type=str, 
		default='trilinear', 
		help='Projection method: trilinear (default), gridding, nn (nearest neighbor)')
	
	parser.add_argument(
		'--display', 
		action="store_true", 
		help='Automatically open montage in e2display')
	
	parser.add_argument(
		'--verbose', "-v", 
		action="store_true", 
		help='Increase verbosity')

	group_viper = parser.add_argument_group(
		title='VIPER options (--mode=viper)',
		description='Projection angles for input images will be extracted from initial-reconstruction results.')
	
	group_viper.add_argument(
		'--classangles', 
		type=str, 
		help='Mode viper, angles file, which will be imported into the header of the class-average stack')
	
	group_viper.add_argument(
		'--classselect', 
		type=str, 
		help='Selection file for included classes.')
	
	group_projmatch = parser.add_argument_group(
		title='Projection-matching options (--mode=projmatch)',
		description='Input image will be aligned to re-projections of the input volume.')
	
	group_projmatch.add_argument(
		'--delta', 
		type=float, 
		help='Angular search increment for reference projections')
	
	group_projmatch.add_argument(
		'--symmetry', 
		type=str, 
		default='c1', 
		help='Symmetry (default c1)')
	
	group_projmatch.add_argument(
		'--matchshift', 
		type=int, 
		default=2, 
		help='Maximum translational shift during multireference alignment (pixels, default 2)')
	
	group_projmatch.add_argument(
		'--matchrad', 
		type=int, 
		default=-1, 
		help='Outer radius for alignment (pixels, default automatic)')
	
	group_projmatch.add_argument(
		'--matchstep', 
		type=int, 
		default=1, 
		help='Radius increment during alignment (pixels, default 1)')
	
	group_meridien = parser.add_argument_group(
		title='MERIDIEN options (--mode=meridien)',
		description='Projection angles from existing refinement will be averages. Outliers can optionally be omitted.')
	
	group_meridien.add_argument(
		'--partangles', 
		type=str, 
		help='Mode meridien, alignment parameters for experimental images')
	
	group_meridien.add_argument(
		'--partselect', 
		type=str, 
		help='Selection file for included particles.')
	
	group_meridien.add_argument(
		'--classdocs', 
		type=str, 
		help='Mode meridien, file pattern for class-to-particle lists')
	
	group_meridien.add_argument(
		'--outliers', 
		type=float, 
		help='Mode meridien, maximum angular difference from initial estimate')
	
	group_meridien.add_argument(
		'--refinerad', 
		type=int, 
		default=-1, 
		help='Outer radius for alignment (pixels, default automatic)')
	
	group_meridien.add_argument(
		'--refineshift', 
		type=int, 
		default=2, 
		help='Maximum translational shift during multireference alignment (pixels, default 2)')
	
	group_meridien.add_argument(
		'--refinestep', 
		type=int, 
		default=1, 
		help='Radius increment during alignment (pixels, default 1)')
	
	group_meridien.add_argument(
		'--align', 
		type=str, 
		help='Mode meridien, 2D alignment method: apsh (default) or scf')
	
	return parser.parse_args()

def main_proj_compare(classavgstack, reconfile, outdir, options, mode='viper', prjmethod='trilinear', 
			 classangles=None, partangles=None, selectdoc=None, 
			 verbose=False, displayYN=False):
	"""
	Main function overseeing various projection-comparison modes.
	
	Arguments:
		classavgstack : Input image stack
		reconfile : Map of which to generate projections (an optionally perform alignment)
		outdir : Output directory
		mode : Mode, viper (pre-existing angles for each input image), projmatch (angles from internal projection-matching)
		verbose : (boolean) Whether to write additional information to screen
		options : (list) Command-line options, run 'sxproj_compare.py -h' for an exhaustive list
		classangles : Angles and shifts for each input class average
		partangles : Angles and shifts for each particle (mode meridien)
		selectdoc : Selection file for included images
		prjmethod : Interpolation method to use
		displayYN : (boolean) Whether to automatically open montage
	"""
	
	# Expand path for outputs
	refprojstack = os.path.join(outdir, 'refproj.hdf')
	refanglesdoc = os.path.join(outdir, 'refangles.txt')
	outaligndoc = os.path.join(outdir, 'docalign2d.txt')

	# If not an input, will create an output, in modes projmatch
	if classangles == None:
		classangles = os.path.join(outdir, 'docangles.txt')
		
		# You need either input angles (mode viper) or to calculate them on the fly (mode projmatch)
		if mode=='viper':
			sp_global_def.ERROR("\nERROR!! Input alignment parameters not specified.", __file__, 1)
			sxprint('Type %s --help to see available options\n' % os.path.basename(__file__))
			exit()
	
	# Check if inputs exist
	check(classavgstack, verbose=verbose)
	check(reconfile, verbose=verbose)
	if verbose: sxprint('')
	
	# Check that dimensions of images and volume agree (maybe rescale volume)
	voldim = EMAN2.EMData(reconfile).get_xsize()
	imgdim = EMAN2.EMData(classavgstack,0).get_xsize()
	if voldim != imgdim:
			sp_global_def.ERROR("\nERROR!! Dimension of input volume doesn't match that of image stack: %s vs. %s" % 
					(voldim, imgdim), __file__, 1)
			
			scale = float(imgdim)/voldim  # only approximate, since full-sized particle radius is arbitrary
			msg  = 'The command to resize the volume will be of the form:\n'
			msg += 'e2proc3d.py %s resized_vol.hdf --scale=%1.5f --clip=%s,%s,%s\n' % (reconfile, scale, imgdim, imgdim, imgdim)
			msg += 'Check the file in the ISAC directory named "README_shrink_ratio.txt" for confirmation.\n'
			sxprint(msg)
			exit()
	
	#  Here if you want to be fancy, there should be an option to chose the projection method,
	#  the mechanism can be copied from sxproject3d.py  PAP
	if prjmethod=='trilinear':
		method_num = 1
	elif prjmethod=='gridding':
		method_num = -1
	elif prjmethod=='nn':
		method_num = 0
	else:
		sp_global_def.ERROR("\nERROR!! Valid projection methods are: trilinear (default), gridding, and nn (nearest neighbor).", __file__, 1)
		sxprint('Usage:\n%s' % USAGE)
		exit()
	
	# Set output directory and log file name
	log, verbose = prepare_outdir_log(outdir, verbose)

	# In case class averages include discarded images, apply selection file
	if mode == 'viper':
		if selectdoc:
			goodavgs, extension = os.path.splitext(os.path.basename(classavgstack))
			newclasses = os.path.join(outdir, goodavgs + "_kept" + extension)
			
			# e2proc2d appends to existing files, so rename existing output
			if os.path.exists(newclasses):
				renamefile = newclasses + '.bak'
				print_log_msg("Selected-classes stack %s exists, renaming to %s" % (newclasses, renamefile), log, verbose)
				print_log_msg("mv %s %s\n" % (newclasses, renamefile), log, verbose)
				os.rename(newclasses, renamefile)
			
			print_log_msg('Creating subset of %s to %s based on selection list %s' % (classavgstack, newclasses, selectdoc), log, verbose)
			cmd = "e2proc2d.py %s %s --list=%s" % (classavgstack, newclasses, selectdoc)
			print_log_msg(cmd, log, verbose)
			os.system(cmd)
			sxprint('')
			
			# Update class-averages
			classavgstack = newclasses
	
	# align de novo to reference map
	if mode=='projmatch':
		# Generate reference projections
		print_log_msg('Projecting %s to output %s using an increment of %s degrees using %s symmetry' % (reconfile, refprojstack, options.delta, options.symmetry), log, verbose)
		cmd = 'sxproject3d.py %s %s --delta=%s --method=S --phiEqpsi=Minus --symmetry=%s' % (reconfile, refprojstack, options.delta, options.symmetry)
		if options.prjmethod == 'trilinear': cmd += ' --trilinear'
		cmd += '\n'
		print_log_msg(cmd, log, verbose)
		project3d(reconfile, refprojstack, delta=options.delta, symmetry=options.symmetry)
		
		# Export projection angles
		print_log_msg("Exporting projection angles from %s to %s" % (refprojstack, refanglesdoc), log, verbose)
		cmd = "sp_header.py %s --params=xform.projection --import=%s\n" % (refprojstack, refanglesdoc)
		print_log_msg(cmd, log, verbose)
		header(refprojstack, 'xform.projection', fexport=refanglesdoc)
		
		# Perform multi-reference alignment
		if options.align=='ali2d':
			projdir = os.path.join(outdir, 'Projdir')  # used if input angles no provided
			if os.path.isdir(projdir):
					print_log_msg('Removing pre-existing directory %s' % projdir, log, verbose)
					print_log_msg('rm -r %s\n' % projdir, log, verbose)
					shutil.rmtree(projdir)  # os.rmdir only removes empty directories
			
			# Zero out alignment parameters in header
			print_log_msg('Zeroing out alignment parameters in header of %s' % classavgstack, log, verbose)
			cmd = 'sxheader.py %s --params xform.align2d --zero\n' % classavgstack
			print_log_msg(cmd, log, verbose)
			header(classavgstack, 'xform.align2d', zero=True)
			
			# Perform multi-reference alignment
			msg = 'Aligning images in %s to projections %s with a radius of %s and a maximum allowed shift of %s' % (classavgstack, refprojstack, options.matchrad, options.matchshift)
			print_log_msg(msg, log, verbose)
			cmd = 'sxmref_ali2d.py %s %s %s --ou=%s --xr=%s --yr=%s\n' % (classavgstack, refprojstack, projdir, options.matchrad, options.matchshift, options.matchshift)
			print_log_msg(cmd, log, verbose)
			mref_ali2d(classavgstack, refprojstack, projdir, ou=options.matchrad, xrng=options.matchshift, yrng=options.matchshift)
			
			# Export alignment parameters
			print_log_msg('Exporting angles from %s into %s' % (classavgstack, classangles), log, verbose)
			cmd = "sp_header.py %s --params=xform.align2d --export=%s\n" % (classavgstack, classangles)
			print_log_msg(cmd, log, verbose)
			header(classavgstack, 'xform.align2d', fexport=classangles)	
		
		# By default, use AP SH
		else:
			apsh(refprojstack, classavgstack, outangles=classangles, refanglesdoc=refanglesdoc, outaligndoc=outaligndoc, 
				outerradius=options.matchrad, maxshift=options.matchshift, ringstep=options.matchstep, log=log, verbose=verbose)
		
		# Diagnostic
		alignlist = read_text_row(classangles)  # contain 2D alignment parameters
		nimg1   = EMAN2.EMUtil.get_image_count(classavgstack)
		assert len(alignlist) == nimg1, "MRK_DEBUG"
		
	# Get alignment parameters from MERIDIEN 
	if mode=='meridien':
		continueTF = True  # Will proceed unless some information is missing
		
		if not partangles:
			sp_global_def.ERROR("\nERROR!! Input alignment parameters not provided.", __file__, 1)
			continueTF = False
	
		if not continueTF:
			sxprint('Type %s --help to see available options\n' % os.path.basename(__file__))
			exit()
		
		if not options.classdocs or options.outliers:
			classdir = os.path.join(outdir, 'Byclass')
			if not os.path.isdir(classdir): os.makedirs(classdir)
			
			if options.outliers: 
				goodclassparttemplate = os.path.join(classdir, 'goodpartsclass{0:03d}.txt')
			else:
				goodclassparttemplate = None
		
			if not options.classdocs:
				classmap = os.path.join(classdir, 'classmap.txt')
				classdoc = os.path.join(classdir, 'docclass{0:03d}.txt')
				options.classdocs = os.path.join(classdir, 'docclass*.txt')
				
				# Separate particles by class
				vomq(classavgstack, classmap, classdoc, log=log, verbose=verbose)
			
		mode_meridien(reconfile, classavgstack, options.classdocs, partangles, selectdoc, options.refineshift, options.refinerad, 
				classangles, outaligndoc, interpolation_method=method_num, outliers=options.outliers,
				goodclassparttemplate=goodclassparttemplate, alignopt=options.align, ringstep=options.refinestep, 
				log=log, verbose=verbose)
		
	# Import Euler angles
	print_log_msg("Importing parameter information into %s from %s" % (classavgstack, classangles), log, verbose)
	cmd = "sp_header.py %s --params=xform.projection --import=%s\n" % (classavgstack, classangles)
	print_log_msg(cmd, log, verbose)
	header(classavgstack, 'xform.projection', fimport=classangles)
	
	# Make comparison stack between class averages (images 0,2,4,...) and re-projections (images 1,3,5,...)
	compstack = compare_projs(reconfile, classavgstack, classangles, outdir, interpolation_method=method_num, log=log, verbose=verbose)

	# Optionally pop up e2display
	if displayYN:
		sxprint('Opening montage')
		cmd = "e2display.py %s\n" % compstack
		sxprint(cmd)
		os.system(cmd)
	
	sxprint("Done!")
	
def check(file, verbose=True):
	"""
	Checks whether file exists.
	
	Arguments:
		file : File to look for
		verbose : (boolean) Whether to write to screen
	"""
	
	if os.path.exists(file):
		if verbose: sxprint("Found %s" % file)
	else:
		sxprint("ERROR!! %s doesn't exist!\n" % file)
		exit()

def prepare_outdir_log(outdir='.', verbose=False, is_main=True):
	"""
	Prepares output directory and sets up log file.
	
	Arguments:
		outdir : Output directory
		verbose : (boolean) Whether to write to screen
		is_main : (boolean) If using multiple cores, some tasks only need to be performed once, not by all cores
	Returns:
		log : instance of Logger class
		verbose : (boolean) New version of Logger can write to screen simultaneously
	"""
	
	# Create directory if it doesn't exist
	if is_main:
		if os.path.isdir(outdir):
			sxprint("Writing to output directory: %s" % outdir)
		else:
			sxprint("Created output directory: %s" % outdir)
			os.makedirs(outdir)  # os.mkdir() can only operate one directory deep
		sp_global_def.write_command(outdir)

	logname = "log_" + datetime.now().strftime("%Y%m%d_%H%M%S") +  ".txt"
	logname = os.path.join(outdir, logname)
	
	#if global_def.LOGFILE: 
		#global_def.LOGFILE = logname
		#print('LOGFILE', global_def.LOGFILE)
		#exit()
	
	# May be using old version of logger.py
	try:
		if verbose:
			log = Logger(base_logger=BaseLogger_Files(), base_logger2=BaseLogger_Print(), file_name=logname)
			verbose = False  # logger output will be echoed to screen
		else:
			log = Logger(base_logger=BaseLogger_Files(), file_name=logname)
	except TypeError:
		if is_main: sxprint("WARNING: Using old sp_logger.py library")
		log = Logger(base_logger=BaseLogger_Files())#, file_name=logname)
		logname = 'log.txt'
		
	if is_main: sxprint("Writing log file to %s\n" % logname)
	
	if is_main:
		progbase = os.path.basename(__file__).split('.')[0].upper()
		#print('progbase', progbase)
		#exit()
		length = len(progbase) + 4
		
		log.add("\n" +
				" "*TIMESTAMP_LENGTH + "*"*length + "\n" +
				" "*TIMESTAMP_LENGTH + "* " + progbase + " *\n" +
				" "*TIMESTAMP_LENGTH + "*"*length)
	
	return log, verbose

def print_log_msg(msg, log=None, verbose=False, is_main=True):
	"""
	Prints messages to log file and, optionally, to the screen.
	
	Arguments:
		msg : Message to write
		log : Instance of Logger class
		verbose : (boolean) Whether to write to screen
		is_main : (boolean) If using MPI, some tasks only need to be performed once, not by all cores
	"""
	
	if is_main:
		if verbose: sxprint(msg)
		if log: log.add(msg)

def apsh(refimgs, imgs2align, outangles=None, refanglesdoc=None, outaligndoc=None, outerradius=-1, 
				maxshift=0, ringstep=1, mode="F", log=None, verbose=False):
	"""
	Generates polar representations of a series of images to be used as alignment references.
	
	Arguments:
		refimgs : Input reference image stack (filename or EMData object)
		imgs2align : Image stack to be aligned (filename or EMData object)
		outangles : Output Euler angles doc file
		refanglesdoc : Input Euler angles for reference projections
		outaligndoc : Output 2D alignment doc file
		outerradius : Outer alignment radius
		maxshift : Maximum shift allowed
		ringstep : Alignment radius step size
		mode : Mode, full circle ("F") vs. half circle ("H")
		log : Logger object
		verbose : (boolean) Whether to write additional information to screen
	"""
	
	# Generate polar representation(s) of reference(s)
	alignrings, polarrefs = mref2polar(refimgs, outerradius=outerradius, ringstep=ringstep, 
			log=log, verbose=verbose)
			 
	# Read image stack (as a filename or already an EMDataobject)
	if isinstance(imgs2align, str): 
		imagestacklist = EMData.read_images(imgs2align)
	else:
		imagestacklist = [imgs2align]
	
	# Get number of images
	numimg = len(imagestacklist)
	
	# Get image dimensions (assuming square, and that images and references have the same dimension)
	idim = imagestacklist[0]['nx']

	# Calculate image center
	halfdim = idim/2 + 1
	
	# Set constants
	currshift = 0
	scale = 1
	
	# Initialize output angles
	outangleslist = []
	outalignlist = []
	
	if outerradius <= 0:
		outerradius = halfdim - 3
		
	# Set search range
	txrng = tyrng = search_range(idim, outerradius, currshift, maxshift)
	
	print_log_msg('Running multireference alignment allowing a maximum shift of %s\n' % maxshift, log, verbose)
	
	# Loop through images
	for imgindex in range(numimg):
		currimg = imagestacklist[imgindex]
		
		# Perform multi-reference alignment (adapted from alignment.mref_ali2d)
		best2dparamslist = [angt, sxst, syst, mirrorfloat, bestreffloat, peakt] = Util.multiref_polar_ali_2d(
			currimg, polarrefs, txrng, tyrng, ringstep, mode, alignrings, halfdim, halfdim)
		bestref = int(bestreffloat)
		mirrorflag = int(mirrorfloat)
		
		# Store parameters
		params2dlist = [angt, sxst, syst, mirrorflag, scale]
		outalignlist.append(params2dlist)
	
		if refanglesdoc:
			refangleslist = read_text_row(refanglesdoc)
			besteulers = refangleslist[bestref]
		else:
			besteulers = [0]*5
		
		# Check for mirroring
		if mirrorflag == 1:
			tempeulers = list(compose_transform3(besteulers[0],besteulers[1],besteulers[2], besteulers[3],besteulers[4],0, 1,
						0,180,0, 0,0,0, 1))
			combinedparams = list(compose_transform3(tempeulers[0],tempeulers[1],tempeulers[2], tempeulers[3],tempeulers[4],0, 1, 
						0,0,-angt, 0,0,0, 1))
		else:
			combinedparams = list(compose_transform3(besteulers[0],besteulers[1],besteulers[2], besteulers[3],besteulers[4],0, 1, 
						0,0,-angt, 0,0,0, 1))
		# compose_transform3: returns phi,theta,psi, tx,ty,tz, scale
		
		outangleslist.append(combinedparams)
		
		# Set transformations as image attribute
		set_params2D(currimg, params2dlist, xform="xform.align2d")  # sometimes I get a vector error with sxheader
		set_params_proj(currimg, besteulers, xform="xform.projection")  # use shifts
		
	if outangles or outaligndoc:
		msg = ''
		if outangles : 
			write_text_row(outangleslist, outangles)
			msg += 'Wrote alignment angles to %s\n' % outangles
			print_log_msg(msg, log, verbose)
		if outaligndoc : 
			write_text_row(outalignlist, outaligndoc)
			msg += 'Wrote 2D alignment parameters to %s\n' % outaligndoc
			print_log_msg(msg, log, verbose)
		
	return outalignlist
	
def mref2polar(refimgs, firstring=1, outerradius=-1, ringstep=1, mode="F", normbysquare=0, log=None, verbose=False):
	"""
	Generates polar representations of a series of images to be used as alignment references.
	
	Arguments:
		refimgs : Input reference image stack (filename or EMData object)
		firstring : Inner alignment radius
		outerradius : Outer alignment radius
		ringstep : Alignment radius step size
		mode : Mode, full circle ("F") vs. half circle ("H)
		normbysquare : If other than 0, normalization by setting the norm to 1
		log : Logger object
		verbose : (boolean) Whether to write additional information to screen
	Returns:
		alignringlist : List of alignment-ring data
		polarreflist : List of polar representation of refernences
	"""
	
	# Read reference stack
	if isinstance(refimgs, str):
		referencelist = EMData.read_images(refimgs)
	else:
		referencelist = [refimgs]  # For single image
	
	numrefs = len(referencelist)
	polarreflist = []
	
	# Get image dimensions (assuming square, and that images and references have the same dimension)
	#print('referencelist', type(refimgs), type(referencelist), type(referencelist[0]))
	#exit()
	idim = referencelist[0]['nx']

	# Calculate image center
	halfdim = idim/2 + 1
	
	if outerradius <= 0:
		outerradius = halfdim - 3
		#print('outerradius1', outerradius)
		
	# Prepare alignment rings
	alignringlist = Numrinit(firstring, outerradius, ringstep, mode)
	
	# Calculate ring weights
	ringweightlist = ringwe(alignringlist, mode)
	
	print_log_msg('Converting %s references to polar coordinates from radius %s to %s with step %s and mode "%s"' % 
			   (numrefs, firstring, outerradius, ringstep, mode), log, verbose)
	
	# Loop through reference images (adapted from sxisac2)
	for refindex in range(numrefs):
		# Convert to polar
		cimage = Util.Polar2Dm(referencelist[refindex] , halfdim, halfdim, alignringlist, mode)
		
		# Fourier transform of rings
		Util.Frngs(cimage, alignringlist)

		# Apply weights to rings
		Util.Applyws(cimage, alignringlist, ringweightlist)

		# Normalize
		Util.Normalize_ring(cimage, alignringlist, normbysquare)  
		#normbysquare: if other than 0, normalizes by setting the norm to 1
		
		# Copy to reference stack
		polarreflist.append(cimage.copy())
		
	return alignringlist, polarreflist

def vomq(classavgstack, classmap, classdoc, log=None, verbose=False):
	"""
	Separate particles according to class assignment.
	
	Arguments:
		classavgstack : Input image stack
		classmap : Output class-to-particle lookup table. Each (long) line contains particles assigned to a class, one file for all classes
		classdoc : Output lists of particles assigned to a class, one file per class
		mode : Mode, viper (pre-existing angles for each input image), projmatch (angles from internal projection-matching)
		log : instance of Logger class
		verbose : (boolean) Whether to write additional information to screen
	"""
	
	# Generate class-to-particle lookup table
	print_log_msg("Exporting members of stack %s to class map %s" % (classavgstack, classmap), log, verbose)
	cmd = "sp_header.py %s --params=members --export=%s" % (classavgstack, classmap) 
	print_log_msg(cmd, log, verbose)
	header(classavgstack, 'members', fexport=classmap)
	
	counter = 0
	
	# Loop through classes
	with open(classmap) as r:
		for idx, line in enumerate(r.readlines()):
			with open(classdoc.format(idx), 'w') as w:
				w.write('\n'.join(line[1:-3].split(', ')))
			counter += 1

	print_log_msg("Wrote %s class selection files\n" % counter, log, verbose)

def mode_meridien(reconfile, classavgstack, classdocs, partangles, selectdoc, maxshift, outerrad, 
					outanglesdoc, outaligndoc, interpolation_method=1, outliers=None, 
					goodclassparttemplate=None, alignopt='apsh', ringstep=1, 
					log=None, verbose=False):
	
	# Resample reference
	recondata = EMAN2.EMData(reconfile)
	idim = recondata['nx']
	reconprep = prep_vol(recondata, npad=2, interpolation_method=interpolation_method)
	
	# Initialize output angles
	outangleslist = []
	outalignlist = []

	# Read class lists
	classdoclist = glob.glob(classdocs)
	partangleslist = read_text_row(partangles)
	
	# Loop through class lists
	for classdoc in classdoclist:  # [classdoclist[32]]:  # 
		# Strip out three-digit filenumber
		classexample = os.path.splitext(classdoc)
		classnum = int(classexample[0][-3:])
		
		# Initial average
		[avg_phi_init, avg_theta_init] = average_angles(partangleslist, classdoc, selectdoc=selectdoc)
		
		# Look for outliers
		if outliers:
			[avg_phi_final, avg_theta_final] = average_angles(partangleslist, classdoc, selectdoc=selectdoc,
					init_angles=[avg_phi_init, avg_theta_init], threshold=outliers, 
					goodpartdoc=goodclassparttemplate.format(classnum), log=log, verbose=verbose)
		else:
			[avg_phi_final, avg_theta_final] = [avg_phi_init, avg_theta_init]
		
		# Compute re-projection
		refprjreal = prgl(reconprep, [avg_phi_final,avg_theta_final,0,0,0], interpolation_method=1, return_real=True)
		
		# Align to class average
		classavg = get_im(classavgstack, classnum)
		
		# Alignment using self-correlation function
		if alignopt == 'scf':
			ang_align2d, sxs, sys, mirrorflag, peak = align2d_scf(classavg, refprjreal, maxshift, maxshift, ou=outerrad)
		
		# Weird results
		elif alignopt == 'align2d':
			# Set search range
			currshift = 0
			txrng = tyrng = search_range(idim, outerrad, currshift, maxshift)
			
			# Perform alignment
			ang_align2d, sxs, sys, mirrorflag, peak = align2d(classavg, refprjreal, txrng, tyrng, last_ring=outerrad)  
		
		# Direct3 (angles seemed to be quantized)
		elif alignopt == 'direct3':
			[[ang_align2d, sxs, sys, mirrorflag, peak]] = align2d_direct3([classavg], refprjreal, maxshift, maxshift, ou=outerrad)
		
		# APSH-like alignment (default)
		else:
			[[ang_align2d, sxs, sys, mirrorflag, scale]] = apsh(refprjreal, classavg, 
				outerradius=outerrad, maxshift=maxshift, ringstep=ringstep)
			
		outalignlist.append([ang_align2d, sxs, sys, mirrorflag, 1])
		msg = "Particle list %s: ang_align2d=%s sx=%s sy=%s mirror=%s\n" % (classdoc, ang_align2d, sxs, sys, mirrorflag)
		print_log_msg(msg, log, verbose)
		
		# Check for mirroring
		if mirrorflag == 1:
			tempeulers = list(compose_transform3(avg_phi_final,avg_theta_final,0, 0,0,0, 1,
						0,180,0, 0,0, 0,1))
			combinedparams = list(compose_transform3(tempeulers[0],tempeulers[1],tempeulers[2], tempeulers[3],tempeulers[4],0, 1, 
						0,0,-ang_align2d, 0,0,0, 1))
		else:
			combinedparams = list(compose_transform3(avg_phi_final,avg_theta_final,0, 0,0,0, 1,
						0,0,-ang_align2d, 0,0, 0,1))
		# compose_transform3: returns phi,theta,psi, tx,ty,tz, scale
		
		outangleslist.append(combinedparams)
	# End class-loop
	
	write_text_row(outangleslist, outanglesdoc)
	write_text_row(outalignlist, outaligndoc)
	print_log_msg('Wrote alignment parameters to %s and %s\n' % (outanglesdoc, outaligndoc), log, verbose)
	
	del recondata  # Clean up
		
	
def average_angles(alignlist, partdoc, selectdoc=None, init_angles=None, threshold=None, goodpartdoc=None, 
				   log=None, verbose=False):
	"""
	Computes a vector average of a set of particles' Euler angles phi and theta.
	
	Arguments:
		alignlist : Alignment parameter doc file, i.e., from MERIDIEN refinment
		partdoc : List of particle indices whose angles should be averaged
		selectdoc : Input substack selection file if particles removed before refinement (e.g., Substack/isac_substack_particle_id_list.txt)
		init_angles : List (2 elements) with initial phi and theta angles, for excluding outliers
		threshold : Angular threshold (degrees) beyond which particles exceeding this angular difference from init_angles will be excluded
		goodpartdoc : Output list of retained particles if a threshold was specified
		log : Logger object
		verbose : (boolean) Whether to write additional information to screen
	Returns:
		list of 2 elements:
			avg_phi
			avg_theta
	"""
	
	# Read alignment parameters
	if isinstance(alignlist, str): alignlist = read_text_row(outaligndoc)
	# (If loading the same parameters repeatedly, better to read the file once externally and pass only the list.)
	
	# Read class list
	partlist = read_text_row(partdoc)
			
	if selectdoc: 
		selectlist = read_text_row(selectdoc)
	else:
		selectlist = None
		
	sum_phi = np.array([0.0,0.0])
	sum_theta = np.array([0.0,0.0])
	totparts = 0
	num_outliers = 0
	goodpartlist = []
	goodpartcounter = 0
	
	# Loop through particles
	for totpartnum in partlist:
		if selectlist:
			goodpartnum = selectlist.index(totpartnum)
		else:
			goodpartnum = totpartnum[0]

		try:
			phi_deg = alignlist[goodpartnum][0]
			theta_deg = alignlist[goodpartnum][1]
			phi_rad = np.deg2rad(phi_deg)
			theta_rad =  np.deg2rad(theta_deg)
		except IndexError:
			msg = "\nERROR!! %s tries to access particle #%s" % (partdoc, goodpartnum)
			numalignparts = len(alignlist)
			msg += "\nAlignment doc file has only %s entries" % (numalignparts)
			msg += "\nMaybe try substack selection file with flag '--select <substack_select>'?"
			sp_global_def.ERROR(msg, __file__, 1)
			exit()
		
		if init_angles:
			angdif = angle_diff(init_angles, [phi_deg, theta_deg])
			if angdif > 180: angdif = 360.0 - angdif
		
		totparts += 1
		
		# Exclude particles exceeding optional threshold
		if threshold==None or angdif < threshold:
			sum_phi += (np.cos(phi_rad), np.sin(phi_rad))
			sum_theta += (np.cos(theta_rad), np.sin(theta_rad))
			goodpartlist.append(goodpartnum)
			goodpartcounter += 1
		else:
			num_outliers += 1
		
	# Compute final average
	avg_phi = degrees(atan2(sum_phi[1],sum_phi[0]))
	avg_theta = degrees(atan2(sum_theta[1],sum_theta[0]))
	
	# Clean up, might reuse
	del alignlist
	del partlist
	del selectlist
	
	msg = "Particle list %s: average angles (%s, %s)" % (partdoc, avg_phi, avg_theta)
	print_log_msg(msg, log, verbose)

	if threshold:
		msg = "Found %s out of %s outliers exceeding an angle difference of %s degrees from initial estimate" % (num_outliers, totparts, threshold)
		print_log_msg(msg, log, verbose)
		
		if goodpartdoc:
			if goodpartcounter > 0:
				write_text_row(goodpartlist, goodpartdoc)
				msg = "Wrote %s particles to %s" % (goodpartcounter, goodpartdoc)
				print_log_msg(msg, log, verbose)
			else:
				msg = "WARNING!! Kept 0 particles from class %s" % partdoc
				print_log_msg(msg, log, verbose)
				[avg_phi, avg_theta] = init_angles
				
	return [avg_phi, avg_theta]
	
def compare_projs(reconfile, classavgstack, inputanglesdoc, outdir, interpolation_method=1, log=None, verbose=False):
	"""
	Make comparison stack between class averages (even-numbered (starts from 0)) and re-projections (odd-numbered).
	
	Arguments:
		reconfile : Input volume from which to generate re-projections
		classavgstack ; Input image stack
		inputanglesdoc : Input Euler angles doc
		outdir ; Output directory
		interpolation_method : Interpolation method: nearest neighbor (nn, 0), trilinear (1, default), gridding (-1)
		log : Logger object
		verbose : (boolean) Whether to write additional information to screen
	Returns:
		compstack : Stack of comparisons between input image stack (even-numbered (starts from 0)) and input volume (odd-numbered)
	"""
	
	recondata = EMAN2.EMData(reconfile)
	nx = recondata.get_xsize()

	# Resample reference
	reconprep = prep_vol(recondata, npad=2, interpolation_method=interpolation_method)

	ccclist = []
	
	#  Here you need actual radius to compute proper ccc's, but if you do, you have to deal with translations, PAP
	mask = model_circle(nx//2-2,nx,nx)
	mask.write_image(os.path.join(outdir, 'maskalign.hdf'))
	compstack = os.path.join(outdir, 'comp-proj-reproj.hdf')
	
	# Number of images may have changed
	nimg1   = EMAN2.EMUtil.get_image_count(classavgstack)
	angleslist = read_text_row(inputanglesdoc)
	
	for imgnum in range(nimg1):
		# Get class average
		classimg = get_im(classavgstack, imgnum)
		
		# Compute re-projection
		prjimg = prgl(reconprep, angleslist[imgnum], interpolation_method=1, return_real=False)
		
		# Calculate 1D power spectra
		rops_dst = rops_table(classimg*mask)  
		rops_src = rops_table(prjimg)
		
		#  Set power spectrum of reprojection to the data.
		#  Since data has an envelope, it would make more sense to set data to reconstruction,
		#  but to do it one would have to know the actual resolution of the data. 
		#  you can check sxprocess.py --adjpw to see how this is done properly  PAP
		table = [0.0]*len(rops_dst)  # initialize table
		for j in range( len(rops_dst) ):
			table[j] = sqrt( rops_dst[j]/rops_src[j] )
		prjimg = fft(filt_table(prjimg, table))  # match FFT amplitudes of re-projection and class average

		cccoeff = ccc(prjimg, classimg, mask)
		#print imgnum, cccoeff
		classimg.set_attr_dict({'cross-corr':cccoeff})
		prjimg.set_attr_dict({'cross-corr':cccoeff})
		
		montagestack = []
		montagestack.append(prjimg)
		montagestack.append(classimg)
		comparison_pair = montage2(montagestack, ncol=2, marginwidth=1)
		comparison_pair.write_image(compstack,imgnum)
		
		ccclist.append(cccoeff)
	del angleslist
	meanccc = sum(ccclist)/nimg1
	print_log_msg("Average CCC is %s\n" % meanccc, log, verbose)
	
	nimg2 = EMAN2.EMUtil.get_image_count(compstack)
	
	for imgnum in range(nimg2):  # xrange will be deprecated in Python3
		prjimg = get_im(compstack,imgnum)
		meanccc1 = prjimg.get_attr_default('mean-cross-corr', -1.0)
		prjimg.set_attr_dict({'mean-cross-corr':meanccc})
		write_header(compstack,prjimg,imgnum)
	
	return compstack
	
def montage2(inputstack, ncol, marginwidth=0, bkgd=0, outfile=None):
	"""
	Generates montage of images into one image.
	Adapted from sxmontage.py
	
	Arguments:
		inputstack : Stack of input images to merge into montage
		ncol : Number of images per row
		marginwidth : Margin width, pixels
		bkgd : Background value of montage
		outfile : Optional output file with montage output
	Returns:
		montage : EMData object of image montage
	"""
	
	if isinstance(inputstack, str): inputstack = EMData.read_images(inputstack)
	
	# Get single-image dimensions
	nx = inputstack[0].get_xsize()
	ny = inputstack[0].get_ysize()
	
	# Get number of images and calculate montage dimensions
	numimgs = len(inputstack)
	numrows = (numimgs-1)/ncol + 1
	
	# Create blank image
	montage_xdim = (nx + marginwidth)*ncol
	montage_ydim = (ny + marginwidth)*numrows
	montage = model_blank(montage_xdim, montage_ydim, 1, bkgd)
	
	# Loop through images
	for imgnum in range(numimgs):
		# Horizontal grid position is image# modulo NCOL
		colnum = imgnum % ncol
		
		# Montage is numbered from the top down
		rownum = numrows - 1 - imgnum/ncol
		
		xoffset = colnum*(nx+marginwidth)
		yoffset = rownum*(ny+marginwidth)
		insert_image(inputstack[imgnum], montage, xoffset, yoffset)
	
	if outfile: montage.write_image(outfile)
	
	return montage
	
def insert_image(smallimg, largeimage, xoffset, yoffset):
	"""
	Inserts small image into large image.
	Adapted from sxmontage.py
	
	Arguments:
		smallimg : Small image to insert into large image
		largeimage : Large image (OVERWRITTEN!) into which small image will be inserted
		xoffset : Top-left x-coordinate of large image where small image will be inserted
		yoffset : Top-left y-coordinate of large image where small image will be inserted
	"""
	
	# Get small-image dimensions
	nx = smallimg.get_xsize()
	ny = smallimg.get_ysize()
	
	for xcoord in range(nx):
		for ycoord in range(ny):
			getpixel = smallimg.get_value_at(xcoord, ycoord)
			largeimage.set_value_at(xoffset+xcoord, yoffset+ycoord, getpixel)

if __name__ == "__main__":
	sp_global_def.print_timestamp( "Start" )
	options = parse_command_line()
	
	##print args, options  # (Everything is in options.)
	#print options  
	#print('LOGFILE',global_def.LOGFILE)
	#exit()
	
	# If output directory not specified, write to same directory as class averages
	if not options.outdir:
		outdir = os.path.dirname(os.path.realpath(options.classavgs))
	else:
		outdir = options.outdir

	if options.mode == 'viper':
		selectdoc = options.classselect
	elif options.mode == 'projmatch':
		selectdoc = None
	elif options.mode == 'meridien':
		selectdoc = options.partselect
	else:
		sp_global_def.ERROR("\nERROR!! Valid mode not specified. Valid modes are: viper, projmatch, and meridien.", __file__, 1)
		sxprint('Type %s --help to see available options\n' % os.path.basename(__file__))
		exit()

	main_proj_compare(options.classavgs, options.vol3d, outdir, options, mode=options.mode, prjmethod=options.prjmethod, 
			  classangles=options.classangles, partangles=options.partangles, selectdoc=selectdoc, 
			  verbose=options.verbose, displayYN=options.display)
	sp_global_def.print_timestamp( "Finish" )
