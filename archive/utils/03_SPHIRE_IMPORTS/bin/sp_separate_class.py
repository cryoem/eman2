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

import EMAN2
import EMAN2_cppwrap
import EMAN2_meta
import datetime
import os
import sp_applications
import sp_global_def
import sp_logger
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT from EMAN2 import EMUtil, EMArgumentParser, EMANVERSION
pass#IMPORTIMPORTIMPORT from sp_applications import header
pass#IMPORTIMPORTIMPORT from datetime import datetime
pass#IMPORTIMPORTIMPORT from sp_logger import Logger, BaseLogger_Files, BaseLogger_Print
pass#IMPORTIMPORTIMPORT import sp_global_def
pass#IMPORTIMPORTIMPORT from sp_global_def import *
#from global_def import ERROR

# Set default values (global variables written in ALL CAPS)
TIMESTAMP_LENGTH = 23
CLASSMAPFILE = "classmap.txt"
CLASSDOCPREFIX = 'docclass'
CLASSSTACKPREFIX = 'stkclass_'
FILTSTACKPREFIX = 'stkflt_'
USAGE = """ 
PURPOSE:
Separate particles according to class assignment:

%s <input_avgs> <input_images> <output_directory> 
General, required command-line parameters:
  1. Input class averages
  2. Input image stack
  3. Output directory
Outputs:
  %s : Class-to-particle lookup table, one file for all classes 
  %s???.txt : List of particles for each class, one file per class
  EMAN2DB/%s???.bdb : Virtual stacks of particles for each class
  EMAN2DB/%s???.bdb : (Optional) virtual stacks of filtered particles for each class

General options:
%s <input_avgs> <input_images> <output_directory> --filt <filter_radius> --shrink <shrink_factor> --pxsz <pixel_size> --verbose
Parameters:
  --filtrad : low-pass filter radius, Angstroms (or if, pxsz not provided, pixels^-1)
  --pxsz : Pixel size, Angstroms
  --shrink : Downsampling factor (e.g., 6 -> 1/6 original size)
  --verbose : Increase verbosity

Modified 2019-03-19

""" % (__file__, CLASSMAPFILE, CLASSDOCPREFIX, CLASSSTACKPREFIX, FILTSTACKPREFIX, __file__)

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
	classmap = os.path.join(outdir, CLASSMAPFILE)
	classdoc = os.path.join(outdir, CLASSDOCPREFIX + '{0:03d}.txt')
	outbdb = os.path.join(outdir, CLASSSTACKPREFIX + '{0:03d}')
	outflt = os.path.join(outdir, FILTSTACKPREFIX + '{0:03d}')
	num_classes = EMAN2_cppwrap.EMUtil.get_image_count(classavgstack)
	
	# Generate class-to-particle lookup table and class-selection lists
	vomq(classavgstack, classmap, classdoc, log=log, verbose=verbose)
	
	print_log_msg("Writing each class to a separate stack", log, verbose)
	if options.shrink: print_log_msg("Will downsample stacks by a factor of %s" % options.shrink, log, verbose)
	
	if options.filtrad: 
		if options.pxsz: 
			filtrad = options.pxsz/options.filtrad
		else:
			filtrad = options.filtrad
		print_log_msg("Will low-pass filter to %s px^-1" % options.shrink, log, verbose)
	
	tot_parts = 0
	
	# Loop through classes
	for class_num in xrange(num_classes):
		# Write BDB stacks
		cmd = "e2bdb.py %s --makevstack bdb:%s --list %s" % (instack, outbdb.format(class_num), classdoc.format(class_num))
		print_log_msg(cmd, log, verbose)
		os.system(cmd)
		num_class_imgs = EMAN2_cppwrap.EMUtil.get_image_count('bdb:' + outbdb.format(class_num))
		tot_parts += num_class_imgs
		
		# Optional filtered stack
		if options.filtrad or options.shrink:
			cmd = "e2proc2d.py bdb:%s bdb:%s --inplace" % (outbdb.format(class_num), outflt.format(class_num))
			# --inplace overwrites existing images
			
			if options.filtrad:
				cmd = cmd + " --process=filter.lowpass.gauss:cutoff_freq=%s" % filtrad
			
			if options.shrink:
				cmd = cmd + " --meanshrink=%s" % options.shrink
				
			print_log_msg(cmd, log, verbose)
			os.system(cmd)
	
	print_log_msg("\nDone! Separated %s particles from %s classes\n" % (tot_parts, num_classes), log, verbose=True)
	
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
			if verbose: sp_global_def.sxprint("Writing to output directory: %s" % outdir)
		else:
			if verbose: sp_global_def.sxprint("Created output directory: %s" % outdir)
			os.makedirs(outdir)  # os.mkdir() can only operate one directory deep
		sp_global_def.write_command(outdir)
			
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
	
	logname = "log_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") +  ".txt"
	logname = os.path.join(outdir, logname)
	
	# May be using old version of logger.py
	try:
		if verbose:
			log = sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), base_logger2=sp_logger.BaseLogger_Print(), file_name=logname)
			verbose = False  # logger output will be echoed to screen
		else:
			log = sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), file_name=logname)
	except TypeError:
		sp_global_def.sxprint("WARNING: Using old sp_logger.py library")
		log = sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files())#, file_name=logname)
		logname = 'log.txt'
		
	sp_global_def.sxprint("Writing log file to %s" % logname)
	
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
		if verbose: sp_global_def.sxprint(msg)
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
	sp_applications.header(classavgstack, 'members', fexport=classmap)
	
	counter = 0
	
	# Loop through classes
	with open(classmap) as r:
		for idx, line in enumerate(r.readlines()):
			with open(classdoc.format(idx), 'w') as w:
				w.write('\n'.join(line[1:-3].split(', ')))
			counter += 1

	print_log_msg("Wrote %s class selection files\n" % counter, log, verbose)


if __name__ == "__main__":
	# Command-line arguments
	parser = EMAN2.EMArgumentParser(usage=USAGE,version=EMAN2_meta.EMANVERSION)
	parser.add_argument('classavgs', help='Input class averages')
	parser.add_argument('instack', help='Input image stack')
	parser.add_argument('outdir', type=str, help='Output directory')
	parser.add_argument('--filtrad', type=float, help='For optional filtered images, low-pass filter radius (1/px or, if pixel size specified, Angstroms)')
	parser.add_argument('--pxsz', type=float, default=None, help='Pixel size, Angstroms (might be downsampled by ISAC)')
	parser.add_argument('--shrink', type=int, help='Optional downsampling factor')
	parser.add_argument('--verbose', "-v", action="store_true", help='Increase verbosity')
	
	(options, args) = parser.parse_args()
	#print('args',args)
	#print('options', options)
	#exit()
	
	separate_class(options.classavgs, options.instack, options, outdir=options.outdir)
