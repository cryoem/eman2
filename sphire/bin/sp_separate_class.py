#!/usr/bin/env python
from __future__ import print_function
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#
#

import os
from EMAN2 import *  # EMUtil, EMArgumentParser, EMANVERSION
from sp_applications import header
from datetime import datetime
from sp_logger import Logger, BaseLogger_Files, BaseLogger_Print
from sp_utilities import get_params2D, read_text_row
from sp_global_def import ERROR
from sp_fundamentals import rot_shift2D, resample
from sp_filter import filt_gaussl

# Set default values (global variables written in ALL CAPS)
TIMESTAMP_LENGTH = 23
DOCFILEDIR = "Docs"
CLASSMAPFILE = "classmap.txt"
BIGSTACKCOPY = 'stack_all'
CLASSDOCPREFIX = 'docclass'
STACKFILEDIR = "Stacks"
CLASSSTACKPREFIX = 'stkclass_'
ALIGNSTACKPREFIX = 'stkalign_'
FILTSTACKPREFIX = 'stkflt_'
COMBINEDPARAMS = 'params_combined.txt'
TESTPARAMSOUT = 'params_testout.txt'

USAGE = """ 
PURPOSE:
Separate particles according to class assignment.

General usage:
  %s <input_avgs> <input_images> <output_directory> 
Required command-line parameters:
  1. Input class averages
  2. Input image stack
  3. Output directory
Outputs:
  %s/%s : Class-to-particle lookup table, one file for all classes 
  %s/%s???.txt : List of particles for each class, one file per class
  EMAN2DB/%s???.bdb : Virtual stacks of particles for each class
  %s/%s???.mrcs : (Optional) stacks of aligned particles for each class
  %s/%s???.mrcs : (Optional) stacks of filtered/shrunken particles for each class

To apply a low-pass (tangent) filter:
  %s <input_avgs> <input_images> <output_directory> --filt <filter_radius> --apix <pixel_size> --verbose
Parameters:
  --filtrad : Low-pass filter radius, Angstroms (or if, apix not provided, pixels^-1)
  --apix : Pixel size of input images (NOT class averages), Angstroms (optional, if not provided, filter radius assumed to be pixels^-1)
  --verbose : Increase verbosity

To downsample the output images:
  %s <input_avgs> <input_images> <output_directory> --filt <filter_radius> --shrink <shrink_factor> --apix <pixel_size> --verbose
Parameter:
  --shrink : Downsampling factor (e.g., 4 -> 1/4 original size)

To apply alignments from ISAC to output image stacks:
  %s <input_avgs> <input_images> <output_directory> --align_isac_dir <ISAC_directory> --verbose
Parameters:
  --align_isac_dir : If applying alignments, directory for ISAC output

To apply ISAC alignments, filter, and shrink:
  %s <input_avgs> <input_images> <output_directory> --align_isac_dir <ISAC_directory> --filt <filter_radius> --shrink <shrink_factor> --verbose

Modified 2019-07-18

""" % (__file__, 
	   DOCFILEDIR, CLASSMAPFILE, DOCFILEDIR, CLASSDOCPREFIX, CLASSSTACKPREFIX, STACKFILEDIR, ALIGNSTACKPREFIX, STACKFILEDIR, FILTSTACKPREFIX, 
	   __file__, __file__, __file__, __file__, )

def separate_class(classavgstack, instack, options, outdir='.', verbose=False):
	"""
	Main function overseeing various projection-comparison modes.
	
	Arguments:
		classavgstack : Input image stack
		classmap : Class-to-particle lookup table. Each (long) line contains particles assigned to a class, one file for all classes
		options : (list) Command-line options, run 'sxproj_compare.py -h' for an exhaustive list
		outdir : Output directory
		verbose : (boolean) Whether to write additional information to screen
	"""
	
	# Set output directory and log file name
	prepare_outdir(outdir, verbose)
	log, verbose = prepare_log(outdir, verbose)
	
	# Expand paths for outputs
	prepare_outdir(os.path.join(outdir, DOCFILEDIR) )
	classmap = os.path.join(outdir, DOCFILEDIR, CLASSMAPFILE)
	classdoc = os.path.join(outdir, DOCFILEDIR, CLASSDOCPREFIX + '{0:03d}.txt')
	stackcp = 'bdb:' + os.path.join(outdir, BIGSTACKCOPY)
	outbdb = os.path.join(outdir, CLASSSTACKPREFIX + '{0:03d}')
	if options.align_isac_dir or options.filtrad or options.shrink:
		prepare_outdir(os.path.join(outdir, STACKFILEDIR) )
		outali = os.path.join(outdir, STACKFILEDIR, ALIGNSTACKPREFIX + '{0:03d}')
		outflt = os.path.join(outdir, STACKFILEDIR, FILTSTACKPREFIX + '{0:03d}' + '.mrcs')
	num_classes = EMUtil.get_image_count(classavgstack)
	
	# Generate class-to-particle lookup table and class-selection lists
	vomq(classavgstack, classmap, classdoc, log=log, verbose=verbose)
	
	if options.filtrad: 
		if options.apix: 
			filtrad = options.apix/options.filtrad
		else:
			filtrad = options.filtrad
		print_log_msg("Will low-pass filter to %s px^-1" % options.filtrad, log, verbose)
	
	if options.shrink: print_log_msg("Will downsample stacks by a factor of %s" % options.shrink, log, verbose)
	
	tot_parts = 0
	
	if options.align_isac_dir:
		init_params_file = os.path.join(options.align_isac_dir, "2dalignment", "initial2Dparams.txt")
		all_params_file = os.path.join(options.align_isac_dir, "all_parameters.txt")
		
		if options.debug:
			num_tot_images = EMUtil.get_image_count(instack)
			print('num_tot_images', num_tot_images)
			
			num_init_params = read_text_row(init_params_file)
			print('num_init_params', len(num_init_params) )
			
		# Copy image stack
		print_log_msg("Copying %s to virtual stack %s" % (instack, stackcp), log, verbose )
		cmd = "e2bdb.py %s --makevstack %s" % (instack, stackcp)
		print_log_msg(cmd, log, verbose)
		os.system(cmd)
		
		# Combine alignment parameters
		combined_params_file = os.path.join(outdir, COMBINEDPARAMS)
		combine_isac_params(options.align_isac_dir, init_params_file, all_params_file, combined_params_file, log, verbose)
		
		# Import alignment parameters
		cmd = "sp_header.py %s --params=xform.align2d --import=%s\n" % (stackcp, combined_params_file) 
		print_log_msg(cmd, log, verbose)
		header(stackcp, 'xform.align2d', fimport=combined_params_file)
		
		if options.debug:
			test_params_file = os.path.join(outdir, TESTPARAMSOUT)
			print_log_msg("Writing imported parameters to %s" % test_params_file)
			header(stackcp, 'xform.align2d', fexport=test_params_file)
		
		stack2split = stackcp
	
	# If not aligning images
	else:
		stack2split = instack
		
	print_log_msg("Writing class stacks", log, verbose)
	if options.align_isac_dir:
		print_log_msg('Writing aligned images', log, verbose)
	
	# Loop through classes
	for class_num in xrange(num_classes):
		
		# Optional alignment
		if options.align_isac_dir:
			
			# Write class stack
			cmd = "e2bdb.py %s --makevstack bdb:%s --list %s" % (stack2split, outbdb.format(class_num), classdoc.format(class_num))
			print_log_msg(cmd, log, verbose)
			os.system(cmd)
			num_class_imgs = EMUtil.get_image_count('bdb:' + outbdb.format(class_num) )
			tot_parts += num_class_imgs

			# Read class stack
			class_stack = EMData.read_images('bdb:' + outbdb.format(class_num) )
			box_size = class_stack[0].get_attr('nx')
			if options.shrink:
				sub_rate = float(1)/options.shrink
				box_size = int(float(box_size)/options.shrink + 0.5)
				if options.debug and class_num == 0:
					print('sub_rate', sub_rate, type(sub_rate), box_size, options.shrink)
			
			# Set filenames
			local_mrcs_path = outali.format(class_num) + ".mrcs"
			local_mrcs = EMData(box_size, box_size, num_class_imgs)
			
			# Loop through images
			for img_num in range(len(class_stack) ):
				img_orig = class_stack[img_num]
				alpha, sx, sy, mirror, scale = get_params2D(img_orig)
				img_ali = rot_shift2D(img_orig, alpha, sx, sy, mirror, scale, "quadratic")
				if options.debug and class_num == 0 and img_num == 0:
					print('img_orig.get_attr_dict0', img_orig.get_attr_dict() )
					img_ali_dict = img_ali.get_attr_dict()
					print('img_ali.get_attr_dict1', img_ali.get_attr_dict() )

				if options.filtrad:
					img_ali = filt_gaussl(img_ali, filtrad)
				
				if options.shrink:
					img_ali = resample(img_ali, sub_rate)
				
				local_mrcs.insert_clip(img_ali, (0, 0, img_num) )
				
			local_mrcs.write_image(local_mrcs_path)
			
		# No aligned images
		else:
			
			# Write class stack
			cmd = "e2bdb.py %s --makevstack bdb:%s --list %s" % (stack2split, outbdb.format(class_num), classdoc.format(class_num))
			print_log_msg(cmd, log, verbose)
			
			os.system(cmd)
			num_class_imgs = EMUtil.get_image_count('bdb:' + outbdb.format(class_num) )
			tot_parts += num_class_imgs

			# Optional filtered stack
			if options.filtrad or options.shrink:
				
				cmd = "e2proc2d.py bdb:%s %s --inplace" % (outbdb.format(class_num), outflt.format(class_num))
				# --inplace overwrites existing images
				
				if options.filtrad:
					cmd = cmd + " --process=filter.lowpass.gauss:cutoff_freq=%s" % filtrad
				
				if options.shrink:
					cmd = cmd + " --meanshrink=%s" % options.shrink
					
				print_log_msg(cmd, log, verbose)
				os.system(cmd)
			
	print_log_msg("\nDone! Separated %s particles from %s classes\n" % (tot_parts, num_classes), log, verbose)  #=True)
	
def prepare_outdir(outdir='.', verbose=False, main=True):
	"""
	Create directory if it doesn't exist.
	
	Arguments:
		outdir : Output directory
		verbose : (boolean) Whether to write additional information to screen
		main : (boolean) Whether main MPI process
	"""
	
	# If using MPI, only need to check once
	if main:
		if os.path.isdir(outdir):
			if verbose: print("Writing to output directory: %s" % outdir)
		else:
			if verbose: print("Created output directory: %s" % outdir)
			os.makedirs(outdir)  # os.mkdir() can only operate one directory deep
			
def prepare_log(outdir='.', verbose=False, main=True):
	"""
	Prepare log file.
	
	Arguments:
		outdir : Output directory
		verbose : (boolean) Whether to write additional information to screen
		main : (boolean) Whether main MPI process
	Returns:
		log : Instance of Logger class
		verbose : (boolean) Updates flag to False if Logger class can mirror output to screen
	"""
	
	TIMESTAMP_LENGTH = 23
	
	logname = "log_" + datetime.now().strftime("%Y%m%d_%H%M%S") +  ".txt"
	logname = os.path.join(outdir, logname)
	
	# May be using old version of logger.py
	try:
		if verbose:
			log = Logger(base_logger=BaseLogger_Files(), base_logger2=BaseLogger_Print(), file_name=logname)
			verbose = False  # logger output will be echoed to screen
		else:
			log = Logger(base_logger=BaseLogger_Files(), file_name=logname)
	except TypeError:
		print("WARNING: Using old sp_logger.py library")
		log = Logger(base_logger=BaseLogger_Files())#, file_name=logname)
		logname = 'log.txt'
		
	print("Writing log file to %s" % logname)
	
	if main:
		progbase = os.path.basename(__file__).split('.')[0].upper()

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
		msg : message to write
		log : instance of Logger class
		verbose : (boolean) whether to write to screen
		is_main : (boolean) if using multiple cores, some tasks only need to be performed once, not by all cores
	"""
	
	if is_main:
		if verbose: print(msg)
		if log: log.add(msg)

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

def combine_isac_params(isac_dir, init_params_file, all_params_file, combined_params_file, log=False, verbose=False):
	"""
	Combines initial and all_params from ISAC.
	
	Arguments:
		isac_dir : ISAC directory
		init_params_file : Prealignment parameters file, ISAC_DIR/2dalignment/initial2Dparams.txt, will be all zeroes if prealignment skipped
		all_params_file : Alignment parameters relative to prealignment, ISAC_DIR/all_parameters.txt
		combined_params_file : Output combined alignment parameters
		log : instance of Logger class
		verbose : (boolean) -- whether to write to screen
	"""
	
	from sp_utilities import combine_params2, read_text_row, write_text_row
	
	# Combine alignment parameters
	init_params_list = read_text_row(init_params_file)
	all_params = read_text_row(all_params_file)
	isac_shrink_path = os.path.join(isac_dir, "README_shrink_ratio.txt")
	isac_shrink_file = open(isac_shrink_path, "r")
	isac_shrink_lines = isac_shrink_file.readlines()
	isac_shrink_ratio = float(isac_shrink_lines[5])
	msg = "Combining alignment parameters from %s and %s, dividing by %s, and writing to %s" % \
		(init_params_file, all_params_file, isac_shrink_ratio, combined_params_file)
	print_log_msg(msg, log, verbose)
	combined_params = []
	
	# Loop through images
	for im in range(len(all_params) ): 
		P = combine_params2(
			init_params_list[im][0], init_params_list[im][1], init_params_list[im][2], init_params_list[im][3],
			all_params[im][0], all_params[im][1]/isac_shrink_ratio, all_params[im][2]/isac_shrink_ratio, all_params[im][3] )
		combined_params.append([P[0], P[1], P[2], P[3], 1.0])
	write_text_row(combined_params, combined_params_file)
	
	return combined_params
	

if __name__ == "__main__":
	# Command-line arguments
	parser = EMArgumentParser(usage=USAGE,version=EMANVERSION)
	parser.add_argument('classavgs', help='Input class averages')
	parser.add_argument('instack', help='Input image stack')
	parser.add_argument('outdir', type=str, help='Output directory')
	parser.add_argument('--align_isac_dir', type=str, default=None, help='ISAC directory, for aligning images')
	parser.add_argument('--filtrad', type=float, help='For optional filtered images, low-pass filter radius (1/px or, if pixel size specified, Angstroms)')
	parser.add_argument('--apix', type=float, default=None, help='Pixel size, Angstroms (might be downsampled by ISAC)')
	parser.add_argument('--shrink', type=int, help='Optional downsampling factor')
	parser.add_argument('--verbose', "-v", action="store_true", help='Increase verbosity')
	parser.add_argument('--debug', action="store_true", help='Debug mode')
	
	(options, args) = parser.parse_args()
	#print('args',args)
	#print('options', options)
	#exit()
	
	separate_class(options.classavgs, options.instack, options, outdir=options.outdir)
