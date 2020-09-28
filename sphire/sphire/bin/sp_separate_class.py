#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div
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

import EMAN2
import EMAN2_cppwrap
import EMAN2_meta
import datetime
import os
from EMAN2 import *  # EMUtil, EMArgumentParser, EMANVERSION
from ..libpy.sp_applications import header, pca, prepare_2d_forPCA
from datetime import datetime
from ..libpy.sp_logger import Logger, BaseLogger_Files, BaseLogger_Print
from ..libpy.sp_utilities import get_params2D, read_text_row, get_im, write_text_row, montage_scale
from ..libpy.sp_global_def import ERROR, write_command
from ..libpy.sp_fundamentals import rot_shift2D, resample
from ..libpy.sp_filter import filt_gaussl
from operator import itemgetter
from ..libpy.sp_statistics import ave_var
from EMAN2db import db_open_dict, db_close_dict
from ..libpy import sp_global_def

# Set default values (global variables written in ALL CAPS)
TIMESTAMP_LENGTH = 23
DOCFILEDIR = "Docs"
CLASSMAPFILE = "classmap.txt"
BIGSTACKCOPY = 'stack_all'
CLASSDOCPREFIX = 'docclass'
STACKFILEDIR = "Stacks"
CLASSORIGPREFIX = 'stkorig_'
CLASSFINALPREFIX = 'stkfilt_'
ALIGNSTACKPREFIX = 'stkalign_'
FILTSTACKPREFIX = 'stkflt_'
COMBINEDPARAMS = 'params_combined.txt'
TESTPARAMSOUT = 'params_testout.txt'

CLASSDOCPREFIX = "docclass"
CLASSSTACKPREFIX = "stkclass_"
FILTSTACKPREFIX = "stkflt_"
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

To compute eigenimages:
  %s <input_avgs> <input_images> <output_directory> --align_isac_dir <ISAC_directory> --nvec <number_of_factors> --verbose
Modified 2019-12-03

""" % (
    __file__,
    CLASSMAPFILE,
    CLASSDOCPREFIX,
    CLASSSTACKPREFIX,
    FILTSTACKPREFIX,
    __file__, __file__, __file__, __file__, __file__,
)


def separate_class(classavgstack, instack, options, outdir=".", verbose=False):
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
	write_command(outdir)
	log, verbose = prepare_log(outdir, verbose)
	prepare_outdir(os.path.join(outdir, DOCFILEDIR))

	# Set filename for copy of input BDB
	bdb_path = os.path.join(outdir, 'EMAN2DB', BIGSTACKCOPY + '.bdb')
	if os.path.exists(bdb_path): os.remove(bdb_path)  # will otherwise merge with pre-existing file
	stackcp = 'bdb:' + outdir + '#' + BIGSTACKCOPY

	# Expand paths for outputs
	classmap = os.path.join(outdir, DOCFILEDIR, CLASSMAPFILE)
	classdoc_template = os.path.join(outdir, DOCFILEDIR, CLASSDOCPREFIX + '{0:03d}.txt')
	class_init_bdb_template = os.path.join(outdir, CLASSORIGPREFIX + '{0:03d}')
	class_filt_bdb_template = os.path.join(outdir, CLASSFINALPREFIX + '{0:03d}')
	if options.align_isac_dir or options.filtrad or options.shrink:
		prepare_outdir(os.path.join(outdir, STACKFILEDIR))
		outali = os.path.join(outdir, STACKFILEDIR, ALIGNSTACKPREFIX + '{0:03d}' + options.format)
		outflt = os.path.join(outdir, STACKFILEDIR, FILTSTACKPREFIX + '{0:03d}' + options.format)
		if options.align_isac_dir:
			chains_params = os.path.join(options.align_isac_dir, 'chains_params.txt')
			classed_imgs = os.path.join(options.align_isac_dir, 'processed_images.txt')

	num_classes = EMUtil.get_image_count(classavgstack)

	# Generate class-to-particle lookup table and class-selection lists
	vomq(classavgstack, classmap, classdoc, log=log, verbose=verbose)

	print_log_msg("Writing each class to a separate stack", log, verbose)
	if options.shrink:
		print_log_msg(
			"Will downsample stacks by a factor of %s" % options.shrink, log, verbose
		)

	if options.filtrad:
		if options.apix:
			filtrad = old_div(options.apix , options.filtrad)
		else:
			filtrad = options.filtrad
		print_log_msg("Will low-pass filter to %s px^-1" % options.filtrad, log, verbose)

	if options.nvec != None and options.align_isac_dir == None:
		sp_global_def.ERROR("\nERROR!! To compute eigenimages, need to specify --align_isac_dir", __file__, 1)
		exit()

	if options.nvec:
		print_log_msg('Writing %s eigenimages per class' % options.nvec, log, verbose)

	tot_parts = 0

	if options.align_isac_dir:
		if options.debug:
			num_tot_images = EMUtil.get_image_count(instack)
			print('num_tot_images', num_tot_images)

		if os.path.basename(classavgstack) == 'ordered_class_averages.hdf' and os.path.exists(chains_params):
			# Make substack with processed images
			print_log_msg(
				"Making substack %s from original stack %s using subset in %s" % (stackcp, instack, classed_imgs), log,
				verbose)
			cmd = "e2bdb.py %s --makevstack %s --list %s" % (instack, stackcp, classed_imgs)
		else:
			# Simply copy image stack
			print_log_msg("Copying %s to virtual stack %s" % (instack, stackcp), log, verbose)
			cmd = "e2bdb.py %s --makevstack %s" % (instack, stackcp)

		print_log_msg(cmd, log, verbose)
		os.system(cmd)

		# Combine alignment parameters
		combined_params_file = os.path.join(outdir, COMBINEDPARAMS)
		combine_isac_params(options.align_isac_dir, classavgstack, chains_params, classed_imgs, classdoc_template,
							combined_params_file, log, verbose)

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

	print_log_msg("Writing %s class stacks" % num_classes, log, verbose)
	if options.align_isac_dir:
		print_log_msg('Writing aligned images', log, verbose)

	# Loop through classes
	for class_num in xrange(num_classes):

		# No aligned images
		if not options.align_isac_dir:
			# Write class stack
			cmd = "e2bdb.py %s --makevstack bdb:%s --list %s" % (
				stack2split, class_init_bdb_template.format(class_num), classdoc_template.format(class_num))
			print_log_msg(cmd, log, verbose)

			os.system(cmd)
			num_class_imgs = EMUtil.get_image_count('bdb:' + class_init_bdb_template.format(class_num))
			tot_parts += num_class_imgs

			# Optional filtered stack
			if options.filtrad or options.shrink:

				cmd = "e2proc2d.py bdb:%s %s --inplace" % (
					class_init_bdb_template.format(class_num), outflt.format(class_num))
				# --inplace overwrites existing images

				if options.filtrad:
					cmd = cmd + " --process=filter.lowpass.gauss:cutoff_freq=%s" % filtrad

				if options.shrink:
					cmd = cmd + " --meanshrink=%s" % options.shrink

				print_log_msg(cmd, log, verbose)
				os.system(cmd)

		# Optionally apply alignment
		else:
			# Write class stack
			class_init_bdb_name = 'bdb:' + class_init_bdb_template.format(class_num)
			cmd = "e2bdb.py %s --makevstack %s --list %s" % (
				stack2split, class_init_bdb_name, classdoc_template.format(class_num))
			print_log_msg(cmd, log, verbose)
			os.system(cmd)

			# Set filenames
			aligned_stack_path = outali.format(class_num)
			class_filt_bdb_name = 'bdb:' + class_filt_bdb_template.format(class_num)

			# PCA works better if application of alignment parameters is done internally.
			# So, alignment will be applied afterward.

			if options.debug and class_num == 0:
				verbosity = True
			else:
				verbosity = False

			if options.nvec != None:
				num_class_imgs = filter_shrink(
					class_init_bdb_name,
					aligned_stack_path,
					class_filt_bdb_name,
					alignYN=False,
					filtrad=filtrad,
					shrink=options.shrink,
					verbose=verbosity
				)

				tmp_classavg, class_stack_list = prepare_2d_forPCA(class_filt_bdb_name, mode='a', CTF=False)
				eig_list = pca(class_stack_list, nvec=options.nvec)
				montage_file = os.path.join(outdir, 'stkeigen.hdf')
				avg_img, var_img = ave_var(class_filt_bdb_name)  # needs to be a BDB, aligned_stack_obj didn't work

			# Not computing eigenimages
			else:
				eig_list = []
				montage_file = os.path.join(outdir, 'stkavgvar.hdf')
				avg_img, var_img = ave_var(class_init_bdb_name)  # needs to be a BDB, aligned_stack_obj didn't work

			montage_list = [avg_img] + [var_img] + eig_list
			montage_row = montage_scale(montage_list, scale=True)
			montage_row.write_image(montage_file, class_num)

			# Apply alignments
			num_class_imgs = filter_shrink(
				class_init_bdb_name,
				aligned_stack_path,
				class_filt_bdb_name,
				alignYN=True,
				filtrad=filtrad,
				shrink=options.shrink,
				verbose=verbosity
			)

			tot_parts += num_class_imgs

	print_log_msg("Done! Separated %s particles from %s classes\n" % (tot_parts, num_classes), log, verbose)  # =True)


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
		progbase = os.path.basename(__file__).split(".")[0].upper()

		length = len(progbase) + 4

		log.add(
			"\n"
			+ " " * TIMESTAMP_LENGTH
			+ "*" * length
			+ "\n"
			+ " " * TIMESTAMP_LENGTH
			+ "* "
			+ progbase
			+ " *\n"
			+ " " * TIMESTAMP_LENGTH
			+ "*" * length
		)

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
        if verbose:
            sp_global_def.sxprint(msg)
        if log:
            log.add(msg)

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

    class_counter = 0

    # Loop through classes
    with open(classmap) as r:
        for class_num, line in enumerate(r.readlines()):
            class_list = line[1:-3].split(', ')
            class_counter += 1
            part_counter += len(class_list)
            write_text_row(class_list, classdoc.format(class_num))
    print_log_msg("Wrote %s class selection files\n" % counter, log, verbose)


def combine_isac_params(isac_dir, classavgstack, chains_params_file, old_combined_parts, classdoc, combined_params_file, log=False, verbose=False):
	"""
	Combines initial and all_params from ISAC.

	Arguments:
		isac_dir : ISAC directory
		classavgstack : Input image stack
		chains_params_file : Input alignment parameters applied to averages in sp_chains
		old_combined_parts
		classdoc
		combined_params_file : Output combined alignment parameters
		log : instance of Logger class
		verbose : (boolean) Whether to write to screen
	"""

	from ..libpy.sp_utilities import combine_params2, read_text_row

	# File-handling
	init_params_file = os.path.join(isac_dir, "2dalignment", "initial2Dparams.txt")
	all_params_file = os.path.join(isac_dir, "all_parameters.txt")
	init_params_list = read_text_row(init_params_file)
	all_params = read_text_row(all_params_file)
	isac_shrink_path = os.path.join(isac_dir, "README_shrink_ratio.txt")
	isac_shrink_file = open(isac_shrink_path, "r")
	isac_shrink_lines = isac_shrink_file.readlines()
	isac_shrink_ratio = float(isac_shrink_lines[5])

	"""
	Three cases:
		1) Using class_averages.hdf
		2) Using ordered_class_averages.hdf, but chains_params.txt doesn't exist
		3) Using ordered_class_averages.hdf and chains_params.txt 
	"""

	msg = "Combining alignment parameters from %s and %s, dividing by %s, and writing to %s" % \
		(init_params_file, all_params_file, isac_shrink_ratio, combined_params_file)

	# Check if using ordered class averages and whether chains_params exists
	if os.path.basename(classavgstack) == 'ordered_class_averages.hdf':
		if not os.path.exists(chains_params_file):
			msg += "WARNING: '%s' does not exist. " % chains_params_file
			msg += "         Using '%s' but alignment parameters correspond to 'class_averages.hdf'.\n" % classavgstack
		else:
			msg = "Combining alignment parameters from %s, %s, and %s, dividing by %s, and writing to %s" % \
				(init_params_file, all_params_file, chains_params_file, isac_shrink_ratio, combined_params_file)

	print_log_msg(msg, log, verbose)

	if os.path.basename(classavgstack) == 'ordered_class_averages.hdf' and os.path.exists(chains_params_file):
		chains_params_list = read_text_row(chains_params_file)
		old_combined_list = read_text_row(old_combined_parts)
		num_classes = EMUtil.get_image_count(classavgstack)
		tmp_combined = []

		# Loop through classes
		for class_num in range(num_classes):
			# Extract members
			image = get_im(classavgstack, class_num)
			members = sorted(image.get_attr("members") )
			old_class_list = read_text_row(classdoc.format(class_num) )
			new_class_list = []

			# Loop through particles
			for idx, im in enumerate(members):
				tmp_par = combine_params2(
					init_params_list[im][0], init_params_list[im][1], init_params_list[im][2], init_params_list[im][3],
					all_params[im][0], old_div( all_params[im][1], isac_shrink_ratio ), old_div(all_params[im][2], isac_shrink_ratio) , all_params[im][3] )

				# Combine with class-average parameters
				P = combine_params2(tmp_par[0], tmp_par[1], tmp_par[2], tmp_par[3],
					chains_params_list[class_num][2], chains_params_list[class_num][3], chains_params_list[class_num][4], chains_params_list[class_num][5])

				tmp_combined.append([ im, P[0], P[1], P[2], P[3] ])

				# Need to update class number in class docs
				old_part_num = old_class_list[idx]

				try:
					new_part_num = old_combined_list.index(old_part_num)
				except ValueError:
					print("Couldn't find particle: class_num %s, old_part_num %s, new_part_num %s" % (class_num, old_part_num[0], new_part_num) )

				new_class_list.append(new_part_num)
			# End particle-loop

			# Overwrite pre-existing class doc
			write_text_row(new_class_list, classdoc.format(class_num) )
		# End class-loop

		# Sort by particle number
		combined_params = sorted( tmp_combined, key=itemgetter(0) )

		# Remove first column
		for row in combined_params: del row[0]

	# Not applying alignments of ordered_class_averages
	else:
		combined_params = []

		# Loop through images
		for im in range(len(all_params) ):
			P = combine_params2(
				init_params_list[im][0], init_params_list[im][1], init_params_list[im][2], init_params_list[im][3],
				all_params[im][0], old_div(all_params[im][1], isac_shrink_ratio), old_div(all_params[im][2],isac_shrink_ratio), all_params[im][3] )
			combined_params.append([P[0], P[1], P[2], P[3], 1.0])

	write_text_row(combined_params, combined_params_file)
	print_log_msg('Wrote %s entries to %s\n' % (len(combined_params), combined_params_file), log, verbose )

	return combined_params

def filter_shrink(input_bdb_name, output_stack_path, output_bdb_name, alignYN=False, filtrad=None, shrink=None, verbose=False):
	"""
	Filters and shrinks image stack.

	Arguments:
		input_bdb_name : Input BDB stack  (in the form bdb:DIRECTORY#STACK)
		output_stack_path : Name for output image stack (MRCS/HDF/etc.)
		output_bdb_name : Name for output BDB stack (in the form bdb:DIRECTORY#STACK)
		alignYN: (boolean) Whether to apply alignments
		filtrad : Filter radius, reciprocal pixels
		shrink : Downsampling factor
		combined_params_file : Output combined alignment parameters
		verbose : (boolean) Whether to write to screen
	"""

	input_stack = EMData.read_images(input_bdb_name)
	num_imgs = EMUtil.get_image_count(input_bdb_name)

	box_size = input_stack[0].get_attr('nx')

	if options.shrink:
		sub_rate = old_div(float(1), options.shrink)
		box_size = int(old_div(float(box_size), options.shrink) + 0.5)
		if verbose : print('sub_rate', sub_rate, type(sub_rate), box_size, options.shrink)

	# Initialize stack & BDB
	aligned_stack_obj = EMData(box_size, box_size, num_imgs)
	new_bdb_dict = db_open_dict(output_bdb_name)

	# Loop through images
	for img_num in range( len(input_stack) ):
		img_orig = input_stack[img_num]

		try:
			alpha, sx, sy, mirror, scale = get_params2D(img_orig)
		except RuntimeError:
			print('\nERROR! Exiting with RuntimeError')
			img_prev = input_stack[img_num-1]
			print('\nPrevious particle: %s %s' % (img_num-1, img_prev.get_attr_dict() ) )
			print('\nCurrent particle: %s %s' % (img_num, img_orig.get_attr_dict() ) )
			exit()

		# Optionally apply alignment parameters
		if alignYN:
			img_ali = rot_shift2D(img_orig, alpha, sx, sy, mirror, scale, "quadratic")
		else:
			img_ali = img_orig

		if verbose and img_num == 0:
			print('\nimg_orig.get_attr_dict0', img_orig.get_attr_dict() )
			img_ali_dict = img_ali.get_attr_dict()
			print('\nimg_ali.get_attr_dict1\n', img_ali.get_attr_dict() )

		if filtrad:
			img_ali = filt_gaussl(img_ali, filtrad)

		if shrink:
			img_ali = resample(img_ali, sub_rate)
			#### (Maybe update resample_ratio)

		img_ali_dict = img_ali.get_attr_dict()
		img_ali_dict["data_path"] = os.path.join('..', STACKFILEDIR, os.path.basename(output_stack_path) )
		img_ali_dict["ptcl_source_coord_id"] = img_num
		new_bdb_dict[img_num] = img_ali_dict

		aligned_stack_obj.insert_clip(img_ali, (0, 0, img_num) )
	# End image-loop

	aligned_stack_obj.write_image(output_stack_path)
	db_close_dict(output_bdb_name)

	return num_imgs

def main():
	# Command-line arguments
	parser = EMArgumentParser(usage=USAGE, version=EMANVERSION)
	parser.add_argument('classavgs', help='Input class averages')
	parser.add_argument('instack', help='Input image stack')
	parser.add_argument('outdir', type=str, help='Output directory')
	parser.add_argument('--align_isac_dir', type=str, default=None,
						help='ISAC directory, for aligning images')
	parser.add_argument('--filtrad', type=float,
						help='For optional filtered images, low-pass filter radius (1/px or, if pixel size specified, Angstroms)')
	parser.add_argument('--apix', type=float, default=None,
						help='Pixel size, Angstroms (might be downsampled by ISAC)')
	parser.add_argument('--shrink', type=int, help='Optional downsampling factor')
	parser.add_argument('--nvec', type=int, help='Number of eigenimages to compute')
	parser.add_argument('--format', type=str, default='.mrcs',
						help='Format of optional output aligned-imaged stacks')
	parser.add_argument('--verbose', "-v", action="store_true", help='Increase verbosity')
	parser.add_argument('--debug', action="store_true", help='Debug mode')

	(options, args) = parser.parse_args()
	# print('args',args)
	# print('options', options)
	# exit()

	separate_class(options.classavgs, options.instack, options, outdir=options.outdir)

if __name__ == "__main__":
	main()
