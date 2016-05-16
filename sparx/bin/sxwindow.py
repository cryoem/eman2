#!/usr/bin/env python
#
# Author: T. Durmaz 08/29/2014 (tunay.durmaz@uth.tmc.edu)
# Copyright (c) 2014 The University of Texas - Houston Medical School
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

from __future__ import print_function
import os, sys
import json
import glob

from optparse import OptionParser, SUPPRESS_HELP
from EMAN2 import *
from EMAN2db import *
from EMAN2jsondb import *
from sparx import *
from applications import MPI_start_end

def check_options(options, progname):

	error_status = None
	if options.coordinates_format.lower() not in ["sparx", "eman1", "eman2", "spider"]:
		error_status = ("Invalid option value: --coordinates_format=%s. Please run %s -h for help." % (options.coordinates_format, progname), getframeinfo(currentframe()))
	if_error_then_all_processes_exit_program(error_status)

	error_status = None
	if options.import_ctf:
		if os.path.exists(options.import_ctf) == False:
			error_status = ("Specified CTER CTF File is not found. Please check --import_ctf option. Run %s -h for help." % (progname), getframeinfo(currentframe()))
	else:
		# assert (not options.import_ctf)
		if options.limit_ctf:
			error_status = ("--limit_ctf option requires valid CTER CTF File (--import_ctf). Please run %s -h for help." % (progname), getframeinfo(currentframe()))
	if_error_then_all_processes_exit_program(error_status)

	error_status = None
	if (options.resample_ratio <= 0.0 or options.resample_ratio > 1.0):
		error_status = ("Invalid option value: --resample_ratio=%s. Please run %s -h for help." % (options.resample_ratio, progname), getframeinfo(currentframe()))
	if_error_then_all_processes_exit_program(error_status)

	
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """  input_micrograph_pattern  input_coordinates_pattern  output_directory  --coordinates_format  --box_size=box_size  --invert  --import_ctf=ctf_file  --limit_ctf  --resample_ratio=resample_ratio  --defocus_error=defocus_error  --astigmatism_error=astigmatism_error
	
Window particles from a micrograph. The coordinates of the particles should be given as input.
Specify name pattern of input micrographs and coordinates files with wild card (*) enclosed by single quotes (') or double quotes (") 
(Note: sxgui.py automatically adds single quotes (')). Use the wild card (*) to specify the place of micrograph ID (e.g. serial number, time stamp, and etc). 
BDB files can not be selected as input micrographs.
	
	sxwindow.py  ./mic*.hdf  info/mic*_info.json  particles  --coordinates_format=eman2  --box_size=64  --invert  --import_ctf=outdir_cter/partres/partres.txt
	
"""
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--coordinates_format",  type="string",        default="eman1",   help="format of input coordinates files: 'sparx', 'eman1', 'eman2', or 'spider'. the coordinates of sparx, eman2, and spider format is particle center. the coordinates of eman1 format is particle box conner associated with the original box size. (default eman1)")
	parser.add_option("--box_size",            type="int",           default=256,       help="x and y dimension of square area to be windowed (in pixels): pixel size after resampling is assumed when resample_ratio < 1.0 (default 256)")
	parser.add_option("--invert",              action="store_true",  default=False,     help="invert image contrast: recommended for cryo data (default False)")
	parser.add_option("--import_ctf",          type="string",        default="",        help="file name of sxcter output: normally partres.txt (default none)") 
	parser.add_option("--limit_ctf",           action="store_true",  default=False,     help="filter micrographs based on the CTF limit: this option requires --import_ctf. (default False)")	
	parser.add_option("--resample_ratio",      type="float",         default=1.0,       help="ratio of new to old image size (or old to new pixel size) for resampling: Valid range is 0.0 < resample_ratio <= 1.0. (default 1.0)")
	parser.add_option("--defocus_error",       type="float",         default=1000000.0, help="defocus errror limit: exclude micrographs whose relative defocus error as estimated by sxcter is larger than defocus_error percent. the error is computed as (std dev defocus)/defocus*100%. (default 1000000.0)" )
	parser.add_option("--astigmatism_error",   type="float",         default=360.0,     help="astigmatism error limit: Set to zero astigmatism for micrographs whose astigmatism angular error as estimated by sxcter is larger than astigmatism_error degrees. (default 360.0)")

	### detect if program is running under MPI
	RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
	
	main_node = 0
	
	if RUNNING_UNDER_MPI:
		from mpi import mpi_init
		from mpi import MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
		
		
		mpi_init(0, [])
		myid = mpi_comm_rank(MPI_COMM_WORLD)
		number_of_processes = mpi_comm_size(MPI_COMM_WORLD)
	else:
		number_of_processes = 1
		myid = 0
	
	(options, args) = parser.parse_args(sys.argv[1:])
	
	error_status = None
	while True:
		if len(args) != 3:
			error_status = ("Please check usage for number of arguments.\n Usage: " + usage + "\n" + "Please run %s -h for help." % (progname), getframeinfo(currentframe()))
			break
		
		mic_pattern = args[0]
		if mic_pattern[:len("bdb:")].lower() == "bdb":
			error_status = ("BDB file can not be selected as input micrographs. Please convert the format, and restart the program. Run %s -h for help." % (progname), getframeinfo(currentframe()))
			break

		if mic_pattern.find("*") == -1:
			error_status = ("Input micrograph file name pattern must contain wild card (*). Please check input_micrograph_pattern argument. Run %s -h for help." % (progname), getframeinfo(currentframe()))
			break

		coords_pattern = args[1]
		if coords_pattern.find("*") == -1:
			error_status = ("Input coordinates file name pattern must contain wild card (*). Please check input_coordinates_pattern argument. Run %s -h for help." % (progname), getframeinfo(currentframe()))
			break

		out_dir = args[2]
		if myid == main_node:
			if os.path.exists(out_dir):
				error_status = ("Output directory exists. Please change the name and restart the program.", getframeinfo(currentframe()))
				break

		break
	if_error_then_all_processes_exit_program(error_status)
	
	# Check invalid conditions of options
	check_options(options, progname)
	
	mic_name_list = None
	error_status = None
	if myid == main_node:
		mic_name_list = glob.glob(mic_pattern)
		if len(mic_name_list) == 0:
			error_status = ("No micrograph file is found. Please check input_micrograph_pattern argument. Run %s -h for help." % (progname), getframeinfo(currentframe()))
	if_error_then_all_processes_exit_program(error_status)
	if RUNNING_UNDER_MPI:
		mic_name_list = wrap_mpi_bcast(mic_name_list, main_node)
	
	coords_name_list = None
	error_status = None
	if myid == main_node:
		coords_name_list = glob.glob(coords_pattern)
		if len(coords_name_list) == 0:
			error_status = ("No coordinates file is found. Please check input_coordinates_pattern argument. Run %s -h for help." % (progname), getframeinfo(currentframe()))
	if_error_then_all_processes_exit_program(error_status)
	if RUNNING_UNDER_MPI:
		coords_name_list = wrap_mpi_bcast(coords_name_list, main_node)
	
##################################################################################################################################################################################################################	
##################################################################################################################################################################################################################	
##################################################################################################################################################################################################################	

	# all processes must have access to indices
	if options.import_ctf:
		i_enum = -1
		i_enum += 1; idx_cter_def          = i_enum # defocus [um]; index must be same as ctf object format
		i_enum += 1; idx_cter_cs           = i_enum # Cs [mm]; index must be same as ctf object format
		i_enum += 1; idx_cter_vol          = i_enum # voltage[kV]; index must be same as ctf object format
		i_enum += 1; idx_cter_apix         = i_enum # pixel size [A]; index must be same as ctf object format
		i_enum += 1; idx_cter_bfactor      = i_enum # B-factor [A^2]; index must be same as ctf object format
		i_enum += 1; idx_cter_ac           = i_enum # amplitude contrast [%]; index must be same as ctf object format
		i_enum += 1; idx_cter_astig_amp    = i_enum # astigmatism amplitude [um]; index must be same as ctf object format
		i_enum += 1; idx_cter_astig_ang    = i_enum # astigmatism angle [degree]; index must be same as ctf object format
		i_enum += 1; idx_cter_sd_def       = i_enum # std dev of defocus [um]
		i_enum += 1; idx_cter_sd_astig_amp = i_enum # std dev of ast amp [A]
		i_enum += 1; idx_cter_sd_astig_ang = i_enum # std dev of ast angle [degree]
		i_enum += 1; idx_cter_cv_def       = i_enum # coefficient of variation of defocus [%]
		i_enum += 1; idx_cter_cv_astig_amp = i_enum # coefficient of variation of ast amp [%]
		i_enum += 1; idx_cter_spectra_diff = i_enum # average of differences between with- and without-astig. experimental 1D spectra at extrema
		i_enum += 1; idx_cter_error_def    = i_enum # frequency at which signal drops by 50% due to estimated error of defocus alone [1/A]
		i_enum += 1; idx_cter_error_astig  = i_enum # frequency at which signal drops by 50% due to estimated error of defocus and astigmatism [1/A]
		i_enum += 1; idx_cter_error_ctf    = i_enum # limit frequency by CTF error [1/A]
		i_enum += 1; idx_cter_mic_name     = i_enum # micrograph name
		i_enum += 1; n_idx_cter            = i_enum
	
	
	# Prepare loop variables
	mic_basename_pattern = os.path.splitext(os.path.basename(mic_pattern))[0]
	coords_format = options.coordinates_format.lower()
	box_size = options.box_size
	box_half = box_size // 2
	mask2d = model_circle(box_size//2, box_size, box_size) # Create circular 2D mask to Util.infomask of particle images
	resample_ratio = options.resample_ratio
	
	n_mic_process = 0
	n_mic_reject_no_coords = 0
	n_mic_reject_no_cter_entry = 0
	n_global_coords_detect = 0
	n_global_coords_process = 0
	n_global_coords_reject_out_of_boundary = 0
	
	serial_id_list = []
	error_status = None
	## not a real while, an if with the opportunity to use break when errors need to be reported
	while myid == main_node:

		# Create list of micrograph serial ID
		# Break micrograph name pattern into prefix and suffix to find the head index of the micrograph serial id
		mic_tokens = mic_pattern.split('*')
		# assert (len(mic_tokens) == 2)
		serial_id_head_index = len(mic_tokens[0])
		# Loop through micrograph names
		for mic_name in mic_name_list:
			# Find the tail index of the serial id and extract serial id from the micrograph name
			serial_id_tail_index = mic_name.index(mic_tokens[1])
			serial_id = mic_name[serial_id_head_index:serial_id_tail_index]
			serial_id_list.append(serial_id)
		# assert (len(serial_id_list) == len(mic_name))
		del mic_name_list # Do not need this anymore
		
		# Load CTFs if necessary
		if options.import_ctf:
			
			ctf_list = read_text_row(options.import_ctf)
			# print("Detected CTF entries : %6d ..." % (len(ctf_list)))
			
			if len(ctf_list) == 0:
				error_status = ("No CTF entry is found in %s. Please check --import_ctf option. Run %s -h for help." % (options.import_ctf, progname), getframeinfo(currentframe()))
				break
			
			if (len(ctf_list[0]) != n_idx_cter):
				error_status = ("Number of columns (%d) must be %d in %s. The format might be old. Please run sxcter.py again." % (len(ctf_list[0]), n_idx_cter, options.import_ctf), getframeinfo(currentframe()))
				break
			
			ctf_dict={}
			n_reject_defocus_error = 0
			ctf_error_limit = [options.defocus_error/100.0, options.astigmatism_error]
			for ctf_params in ctf_list:
				assert(len(ctf_params) == n_idx_cter)
				# mic_basename is name of micrograph minus the path and extension
				mic_basename = os.path.splitext(os.path.basename(ctf_params[idx_cter_mic_name]))[0]
				if(ctf_params[idx_cter_sd_def] / ctf_params[idx_cter_def] > ctf_error_limit[0]):
					print("Defocus error %f exceeds the threshold. Micrograph %s is rejected." % (ctf_params[idx_cter_sd_def] / ctf_params[idx_cter_def], mic_basename))
					n_reject_defocus_error += 1
				else:
					if(ctf_params[idx_cter_sd_astig_ang] > ctf_error_limit[1]):
						ctf_params[idx_cter_astig_amp] = 0.0
						ctf_params[idx_cter_astig_ang] = 0.0
					ctf_dict[mic_basename] = ctf_params
			del ctf_list # Do not need this anymore
		
		break
		
	if_error_then_all_processes_exit_program(error_status)

	if options.import_ctf:
		if options.limit_ctf:
			cutoff_histogram = []  #@ming compute the histogram for micrographs cut of by ctf_params limit.
	
##################################################################################################################################################################################################################	
##################################################################################################################################################################################################################	
##################################################################################################################################################################################################################	
	
	restricted_serial_id_list = []
	if myid == main_node:
		# Loop over serial IDs of micrographs
		for serial_id in serial_id_list:
			# mic_basename is name of micrograph minus the path and extension
			# Here, assuming micrograph and coordinates have the same file basename
			mic_basename = mic_basename_pattern.replace("*", serial_id)
			mic_name = mic_pattern.replace("*", serial_id)
			coords_name = coords_pattern.replace("*", serial_id)
			
			########### # CHECKS: BEGIN
			if coords_name not in coords_name_list:
				print("    Cannot read %s. Skipping %s ..." % (coords_name, mic_basename))
				n_mic_reject_no_coords += 1
				continue
			
			# IF mic is in CTER results
			if options.import_ctf:
				if mic_basename not in ctf_dict:
					print("    Is not listed in CTER results. Skipping %s ..." % (mic_basename))
					n_mic_reject_no_cter_entry += 1
					continue
				else:
					ctf_params = ctf_dict[mic_basename]
			# CHECKS: END
			
			n_mic_process += 1
			
			restricted_serial_id_list.append(serial_id)
		# restricted_serial_id_list = restricted_serial_id_list[:128]  ## for testing against the nonMPI version

	
	if myid != main_node:
		if options.import_ctf:
			ctf_dict = None

	error_status = None
	if len(restricted_serial_id_list) < number_of_processes:
		error_status = ('Number of processes (%d) supplied by --np in mpirun cannot be greater than %d (number of micrographs that satisfy all criteria to be processed) ' % (number_of_processes, len(restricted_serial_id_list)), getframeinfo(currentframe()))
	if_error_then_all_processes_exit_program(error_status)

	## keep a copy of the original output directory where the final bdb will be created
	original_out_dir = out_dir
	if RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD)
		restricted_serial_id_list = wrap_mpi_bcast(restricted_serial_id_list, main_node)
		mic_start, mic_end = MPI_start_end(len(restricted_serial_id_list), number_of_processes, myid)
		restricted_serial_id_list_not_sliced = restricted_serial_id_list
		restricted_serial_id_list = restricted_serial_id_list[mic_start:mic_end]
	
		if options.import_ctf:
			ctf_dict = wrap_mpi_bcast(ctf_dict, main_node)

		# generate subdirectories of out_dir, one for each process
		out_dir = os.path.join(out_dir,"%03d"%myid)
	
	if myid == main_node:
		print("Micrographs processed by main process (including percent complete):")

	len_processed_by_main_node_divided_by_100 = len(restricted_serial_id_list)/100.0

##################################################################################################################################################################################################################	
##################################################################################################################################################################################################################	
##################################################################################################################################################################################################################	
#####  Starting main parallel execution

	for my_idx, serial_id in enumerate(restricted_serial_id_list):
		mic_basename = mic_basename_pattern.replace("*", serial_id)
		mic_name = mic_pattern.replace("*", serial_id)
		coords_name = coords_pattern.replace("*", serial_id)

		if myid == main_node:
			print(mic_name, " ---> % 2.2f%%"%(my_idx/len_processed_by_main_node_divided_by_100))
		mic_img = get_im(mic_name)

		ctf_params = ctf_dict[mic_basename]

		# Read coordinates according to the specified format and 
		# make the coordinates the center of particle image 
		if coords_format == "sparx":
			coords_list = read_text_row(coords_name)
		elif coords_format == "eman1":
			coords_list = read_text_row(coords_name)
			for i in xrange(len(coords_list)):
				coords_list[i] = [(coords_list[i][0] + coords_list[i][2] // 2), (coords_list[i][1] + coords_list[i][3] // 2)]
		elif coords_format == "eman2":
			coords_list = js_open_dict(coords_name)["boxes"]
			for i in xrange(len(coords_list)):
				coords_list[i] = [coords_list[i][0], coords_list[i][1]]
		elif coords_format == "spider":
			coords_list = read_text_row(coords_name)
			for i in xrange(len(coords_list)):
				coords_list[i] = [coords_list[i][2], coords_list[i][3]]
			# else: assert (False) # Unreachable code
		
		# Calculate the new pixel size
		if options.import_ctf:
			pixel_size_origin = ctf_params[idx_cter_apix]
			
			if resample_ratio < 1.0:
				# assert (resample_ratio > 0.0)
				new_pixel_size = pixel_size_origin / resample_ratio
				print("Resample micrograph to pixel size %6.4f and window segments from resampled micrograph." % new_pixel_size)
			else:
				# assert (resample_ratio == 1.0)
				new_pixel_size = pixel_size_origin
		
			# Set ctf along with new pixel size in resampled micrograph
			ctf_params[idx_cter_apix] = new_pixel_size
		else:
			# assert (not options.import_ctf)
			if resample_ratio < 1.0:
				# assert (resample_ratio > 0.0)
				print("Resample micrograph with ratio %6.4f and window segments from resampled micrograph." % resample_ratio)
			# else:
			#	assert (resample_ratio == 1.0)
		
		# Apply filters to micrograph
		fftip(mic_img)
		if options.limit_ctf:
			# assert (options.import_ctf)
			# Cut off frequency components higher than CTF limit 
			q1, q2 = ctflimit(box_size, ctf_params[idx_cter_def], ctf_params[idx_cter_cs], ctf_params[idx_cter_vol], new_pixel_size)
			
			# This is absolute frequency of CTF limit in scale of original micrograph
			if resample_ratio < 1.0:
				# assert (resample_ratio > 0.0)
				q1 = resample_ratio * q1 / float(box_size) # q1 = (pixel_size_origin / new_pixel_size) * q1/float(box_size)
			else:
				# assert (resample_ratio == 1.0) -> pixel_size_origin == new_pixel_size -> pixel_size_origin / new_pixel_size == 1.0
				q1 = q1 / float(box_size)
			
			if q1 < 0.5:
				mic_img = filt_tanl(mic_img, q1, 0.01)
				cutoff_histogram.append(q1)
		
		# Cut off frequency components lower than the box size can express 
		mic_img = fft(filt_gaussh(mic_img, resample_ratio / box_size))
		
		
		# Resample micrograph, map coordinates, and window segments from resampled micrograph using new coordinates
		# after resampling by resample_ratio, new pixel size will be pixel_size/resample_ratio = new_pixel_size
		# NOTE: 2015/04/13 Toshio Moriya
		# resample() efficiently takes care of the case resample_ratio = 1.0 but
		# it does not set apix_*. Even though it sets apix_* when resample_ratio < 1.0 ...
		mic_img = resample(mic_img, resample_ratio)
		
		if options.invert:
			mic_stats = Util.infomask(mic_img, None, True) # mic_stat[0:mean, 1:SD, 2:min, 3:max]
			Util.mul_scalar(mic_img, -1.0)
			mic_img += 2 * mic_stats[0]
		
		if options.import_ctf:
			from utilities import generate_ctf
			ctf_obj = generate_ctf(ctf_params) # indexes 0 to 7 (idx_cter_def to idx_cter_astig_ang) must be same in cter format & ctf object format.
		
		# Prepare loop variables
		nx = mic_img.get_xsize() 
		ny = mic_img.get_ysize()
		x0 = nx//2
		y0 = ny//2

		n_coords_reject_out_of_boundary = 0
		local_stack_name  = "bdb:%s#" % out_dir + mic_basename + '_ptcls'
		local_particle_id = 0 # can be different from coordinates_id
		# Loop over coordinates
		for coords_id in xrange(len(coords_list)):
			
			x = int(coords_list[coords_id][0])
			y = int(coords_list[coords_id][1])
			
			if resample_ratio < 1.0:
				# assert (resample_ratio > 0.0)
				x = int(x * resample_ratio)	
				y = int(y * resample_ratio)
			# else:
			# 	assert(resample_ratio == 1.0)
				
			if( (0 <= x - box_half) and ( x + box_half <= nx ) and (0 <= y - box_half) and ( y + box_half <= ny ) ):
				particle_img = Util.window(mic_img, box_size, box_size, 1, x-x0, y-y0)
			else:
				print("Coordinates ID = %04d (x = %4d, y = %4d, box_size = %4d) is out of micrograph bound, skipping ..." % (coords_id, x, y, box_size))
				n_coords_reject_out_of_boundary += 1
				continue
			
			particle_img = ramp(particle_img)
			particle_stats = Util.infomask(particle_img, mask2d, False) # particle_stats[0:mean, 1:SD, 2:min, 3:max]
			particle_img -= particle_stats[0]
			particle_img /= particle_stats[1]
			
			# NOTE: 2015/04/09 Toshio Moriya
			# ptcl_source_image might be redundant information ...
			# Consider re-organizing header entries...
			particle_img.set_attr("ptcl_source_image", mic_name)
			particle_img.set_attr("ptcl_source_coord_id", coords_id)
			particle_img.set_attr("ptcl_source_coord", [int(coords_list[coords_id][0]), int(coords_list[coords_id][1])])
			particle_img.set_attr("resample_ratio", resample_ratio)
			
			# NOTE: 2015/04/13 Toshio Moriya
			# apix_* attributes are updated by resample() only when resample_ratio != 1.0
			# Let's make sure header info is consistent by setting apix_* = 1.0 
			# regardless of options, so it is not passed down the processing line
			particle_img.set_attr("apix_x", 1.0)
			particle_img.set_attr("apix_y", 1.0)
			particle_img.set_attr("apix_z", 1.0)
			if options.import_ctf:
				particle_img.set_attr("ctf",ctf_obj)
				particle_img.set_attr("ctf_applied", 0)
				particle_img.set_attr("pixel_size_origin", pixel_size_origin)
				# particle_img.set_attr("apix_x", new_pixel_size)
				# particle_img.set_attr("apix_y", new_pixel_size)
				# particle_img.set_attr("apix_z", new_pixel_size)
			# NOTE: 2015/04/13 Toshio Moriya 
			# Pawel Comment: Micrograph is not supposed to have CTF header info.
			# So, let's assume it does not exist & ignore its presence.
			# Note that resample() "correctly" updates pixel size of CTF header info if it exists
			# elif (particle_img.has_ctff()):
			# 	assert(not options.import_ctf)
			# 	ctf_origin = particle_img.get_attr("ctf_obj")
			# 	pixel_size_origin = round(ctf_origin.apix, 5) # Because SXCTER ouputs up to 5 digits 
			# 	particle_img.set_attr("apix_x",pixel_size_origin)
			# 	particle_img.set_attr("apix_y",pixel_size_origin)
			# 	particle_img.set_attr("apix_z",pixel_size_origin)	
			
			# print("local_stack_name, local_particle_id", local_stack_name, local_particle_id)
			particle_img.write_image(local_stack_name, local_particle_id)
			local_particle_id += 1
		
		n_global_coords_detect += len(coords_list)
		n_global_coords_process += local_particle_id
		n_global_coords_reject_out_of_boundary += n_coords_reject_out_of_boundary
		
#		# MRK_DEBUG: Toshio Moriya 2016/05/03
#		# Following codes are for debugging bdb. Delete in future
#		result = db_check_dict(local_stack_name)
#		print('# MRK_DEBUG: result = db_check_dict(local_stack_name): %s' % (result))
#		result = db_list_dicts('bdb:%s' % out_dir)
#		print('# MRK_DEBUG: result = db_list_dicts(out_dir): %s' % (result))
#		result = db_get_image_info(local_stack_name)
#		print('# MRK_DEBUG: result = db_get_image_info(local_stack_name)', result)
		
		# Release the data base of local stack from this process
		# so that the subprocess can access to the data base
		db_close_dict(local_stack_name)
		
#		# MRK_DEBUG: Toshio Moriya 2016/05/03
#		# Following codes are for debugging bdb. Delete in future
#		cmd_line = "e2iminfo.py %s" % (local_stack_name)
#		print('# MRK_DEBUG: Executing the command: %s' % (cmd_line))
#		cmdexecute(cmd_line)
		
#		# MRK_DEBUG: Toshio Moriya 2016/05/03
#		# Following codes are for debugging bdb. Delete in future
#		cmd_line = "e2iminfo.py bdb:%s#data" % (out_dir)
#		print('# MRK_DEBUG: Executing the command: %s' % (cmd_line))
#		cmdexecute(cmd_line)
		
	if RUNNING_UNDER_MPI:
		if options.import_ctf:
			if options.limit_ctf:
				cutoff_histogram = wrap_mpi_gatherv(cutoff_histogram, main_node)

	if myid == main_node:
		if options.limit_ctf:
			# Print out the summary of CTF-limit filtering
			print(" ")
			print("Global summary of CTF-limit filtering (--limit_ctf) ...")
			print("Percentage of filtered micrographs: %8.2f\n" % (len(cutoff_histogram) * 100.0 / len(restricted_serial_id_list_not_sliced)))

			n_bins = 10
			if len(cutoff_histogram) >= n_bins:
				from statistics import hist_list
				cutoff_region, cutoff_counts = hist_list(cutoff_histogram, n_bins)
				print("      Histogram of cut-off frequency")
				print("      cut-off       counts")
				for bin_id in xrange(n_bins):
					print(" %14.7f     %7d" % (cutoff_region[bin_id], cutoff_counts[bin_id]))
			else:
				print("The number of filtered micrographs (%d) is less than the number of bins (%d). No histogram is produced." % (len(cutoff_histogram), n_bins))
	
	n_mic_process = mpi_reduce(n_mic_process, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	n_mic_reject_no_coords = mpi_reduce(n_mic_reject_no_coords, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	n_mic_reject_no_cter_entry = mpi_reduce(n_mic_reject_no_cter_entry, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	n_global_coords_detect = mpi_reduce(n_global_coords_detect, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	n_global_coords_process = mpi_reduce(n_global_coords_process, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	n_global_coords_reject_out_of_boundary = mpi_reduce(n_global_coords_reject_out_of_boundary, 1, MPI_INT, MPI_SUM, main_node, MPI_COMM_WORLD)
	
	# Print out the summary of all micrographs
	if main_node == myid:
		print(" ")
		print("Global summary of micrographs ...")
		print("Detected                        : %6d" % (len(restricted_serial_id_list_not_sliced)))
		print("Processed                       : %6d" % (n_mic_process))
		print("Rejected by no coordinates file : %6d" % (n_mic_reject_no_coords))
		print("Rejected by no CTER entry       : %6d" % (n_mic_reject_no_cter_entry))
		print(" ")
		print("Global summary of coordinates ...")
		print("Detected                        : %6d" % (n_global_coords_detect))
		print("Processed                       : %6d" % (n_global_coords_process))
		print("Rejected by out of boundary     : %6d" % (n_global_coords_reject_out_of_boundary))
		# print(" ")
		# print("DONE!!!")
	
	mpi_barrier(MPI_COMM_WORLD)
	
	if main_node == myid:
	
		import time
		time.sleep(1)
		print("\n Creating bdb:%s/data\n"%original_out_dir)
		for proc_i in range(number_of_processes):
			mic_start, mic_end = MPI_start_end(len(restricted_serial_id_list_not_sliced), number_of_processes, proc_i)
			for serial_id in restricted_serial_id_list_not_sliced[mic_start:mic_end]:
				e2bdb_command = "e2bdb.py "
				mic_basename = mic_basename_pattern.replace("*", serial_id)
				if number_of_processes > 1:
					e2bdb_command += "bdb:" + os.path.join(original_out_dir,"%03d/"%proc_i) + mic_basename + "_ptcls "
				else:
					e2bdb_command += "bdb:" + os.path.join(original_out_dir, mic_basename + "_ptcls ") 
				
				e2bdb_command += " --appendvstack=bdb:%s/data  1>/dev/null"%original_out_dir
				cmdexecute(e2bdb_command, printing_on_success = False)
				
		print("Done!\n")
				
	if RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD)
		from mpi import mpi_finalize
		mpi_finalize()

	sys.stdout.flush()
	sys.exit(0)

if __name__=="__main__":
	main()
