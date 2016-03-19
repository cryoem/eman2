#!/usr/bin/env python

#
# Author: Horatiu Voicu, 2016-02-23--16-04-47-612 (horatiu.voicu@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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

import os
import sys


import global_def
from   global_def import *
from optparse import OptionParser, SUPPRESS_HELP

from EMAN2 import *
from sparx import *

import datetime

from mpi   import  *
from math  import  *

from utilities import send_string_to_all

def shrink_or_enlarge_or_stay_the_same(images, shrink_ratio, target_nx, target_radius, newx, nima):
	from fundamentals import resample
	from utilities import pad
	
	if newx >= target_nx  :
		msk = model_circle(target_radius, target_nx, target_nx)
	else:
		msk = model_circle(newx//2-2, newx,  newx)
		
	do_not_stay_the_same_size = shrink_ratio != 1.0 
	cut_to_window = newx > target_nx 
	need_padding = newx < target_nx 

	for im in xrange(nima):
		if do_not_stay_the_same_size : images[im]  = resample(images[im], shrink_ratio)
		if cut_to_window : images[im] = Util.window(images[im], target_nx, target_nx, 1)
		
		p = Util.infomask(images[im], msk, False)
		images[im] -= p[0]
		p = Util.infomask(images[im], msk, True)
		images[im] /= p[1]
		
		if need_padding : images[im] = pad(images[im], target_nx, target_nx, 1, 0.0)

def main(args):
	
	from alignment import align2d

	progname = os.path.basename(sys.argv[0])
	usage = ( progname + " stack_file isac_directory --radius=particle_radius")
	
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--radius",                type="int",           help="particle radius: there is no default, a sensible number has to be provided, units - pixels (default required int)")
	parser.add_option("--CTF",                   action="store_true",  default=False,      help="apply phase-flip for CTF correction: if set the data will be phase-flipped using CTF information included in image headers (default False)")

	##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	parser.add_option("--return_options", action="store_true", dest="return_options", default=False, help = SUPPRESS_HELP)

	required_option_list = ['radius']
	(options, args) = parser.parse_args(args)

	if options.return_options:
		return parser
	
	if len(args) > 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		sys.exit()
	
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	global_def.BATCH = True

	main_node = 0
	mpi_init(0, [])
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	comm = MPI_COMM_WORLD

	CTF = options.CTF

	error_status = None
	# Making sure all required options appeared.
	for required_option in required_option_list:
		if not options.__dict__[required_option]:
			error_status = ("\n ==%s== mandatory option is missing.\n"%required_option + "Please run '" + progname + " -h' for detailed options", getframeinfo(currentframe()))

	if_error_then_all_processes_exit_program(error_status)
	
	radi  = options.radius
	if(radi < 1):  ERROR("Particle radius geater than zero has to be provided!","sxisac",1,myid)

	command_line_provided_stack_filename = args[0]

	# in this directory 'class_averages.hdf' must be found!
	masterdir = args[1]
	if masterdir == "" or masterdir == None: ERROR("Isac directory has to be provided!","sxisac",1,myid)
	# masterdir = send_string_to_all(masterdir)

	class_averages_file_name = os.path.join(masterdir, "class_averages.hdf")
	if not os.path.exists(class_averages_file_name): ERROR("Isac directory does not contain class_averages.hdf!","sxisac_post_processing",1,myid)

	if (myid == main_node):
		isac_post_processing_output_directory = os.path.join(masterdir,"post_processing_{}".format(datetime.datetime.now().strftime('%Y-%m-%d--%H-%M-%S')))
		cmd = "{} {}".format("mkdir", isac_post_processing_output_directory)
		cmdexecute(cmd)
	else: isac_post_processing_output_directory = ""
	
	isac_post_processing_output_directory = send_string_to_all(isac_post_processing_output_directory)

	size_adjusted_class_averages_file_name = os.path.join(isac_post_processing_output_directory, "radius_adjusted_class_averages.hdf")

# pseudo-code:
########################################################################################################################
### start: 1. take class averages, distribute them among processors and scale them back to original size 

	if myid == main_node:	print "stack:::::::", class_averages_file_name ; total_class_averages_nima = EMUtil.get_image_count(class_averages_file_name)
	else:					total_class_averages_nima = 0
	total_class_averages_nima = bcast_number_to_all(total_class_averages_nima, source_node = main_node)
	list_of_particles = range(total_class_averages_nima)

	image_start, image_end = MPI_start_end(total_class_averages_nima, nproc, myid)
	list_of_particles = list_of_particles[image_start: image_end]
	
	class_averages_nima = len(list_of_particles)

	error_status = None
	if myid == main_node:
		try:
			shrink_ratio_file = os.path.join(masterdir, "README_shrink_ratio.txt")
			fp = open(shrink_ratio_file, "r")
			for line in fp:
				if "--------" in line:
					line = fp.next()
					shrink_ratio = 1/float(line)
					# line = fp.readline()
					# target_radius = int(line)
					break
			else:
				error_status = ("Could not find shrink ratio in '%s', exiting."%shrink_ratio_file, getframeinfo(currentframe()))
		except:
			error_status = ("Could not obtain shrink ratio from isac output directory!", getframeinfo(currentframe()))
	else:
		shrink_ratio = 0

	if_error_then_all_processes_exit_program(error_status)
	
	shrink_ratio = float(mpi_bcast(shrink_ratio, 1, MPI_FLOAT, main_node, MPI_COMM_WORLD)[0])
	
	class_average_images = EMData.read_images(class_averages_file_name, list_of_particles)
	nx = class_average_images[0].get_xsize()
	ref_image = EMData.read_images(command_line_provided_stack_filename, [0])[0]
	target_nx = ref_image.get_xsize()
	newx = int(nx*shrink_ratio + 0.5)

	shrink_or_enlarge_or_stay_the_same(class_average_images, shrink_ratio, target_nx, radi, newx, class_averages_nima)  ## radi is target_radius
	
	mpi_barrier(MPI_COMM_WORLD)
	
### end: 1. take class averages, scale them back to original size
########################################################################################################################

########################################################################################################################
### start: 2. align original images to class averages from isac using exhaustive search (class averages from isac have original image numbers)

	original_images_grouped_by_class_averages = [None]*class_averages_nima
	for class_avg_img_iter in xrange(class_averages_nima):
		original_images_id_list_associated_with_this_class = class_average_images[class_avg_img_iter].get_attr("members")
		original_images_grouped_by_class_averages[class_avg_img_iter] = EMData.read_images(command_line_provided_stack_filename, original_images_id_list_associated_with_this_class)

	if myid == main_node:
		ima = EMData()
		ima.read_image(command_line_provided_stack_filename, list_of_particles[0], True)
		nx = ima.get_xsize()
		del ima
	else:
		nx = 0
	nx = bcast_number_to_all(nx, source_node = main_node)

	#  must know the total number of images on each processor for gathering at root the 2D parameters
	total_number_of_original_images_on_this_processor = 0
	for class_avg_img_iter in xrange(class_averages_nima):
		total_number_of_original_images_on_this_processor += len(original_images_grouped_by_class_averages[class_avg_img_iter])

	list_used_for_sending_parameters = [[] for i in xrange(total_number_of_original_images_on_this_processor)]
	list_count = 0
 
	## ensure that exhaustive search is performed
	range_lim = target_nx//2 + 1 - radi

	make_a_list_to_send_to_master__original_images_id_list_associated_with_this_class = []
	for class_avg_img_iter in xrange(class_averages_nima):
		this_class_average_params = align2d_direct3(original_images_grouped_by_class_averages[class_avg_img_iter], class_average_images[class_avg_img_iter], xrng=range_lim, yrng=range_lim, psimax=180, psistep=1, ou = radi, CTF = CTF)
		original_images_id_list_associated_with_this_class = class_average_images[class_avg_img_iter].get_attr("members")
		make_a_list_to_send_to_master__original_images_id_list_associated_with_this_class.append(original_images_id_list_associated_with_this_class)
		for img_iter in xrange(len(original_images_grouped_by_class_averages[class_avg_img_iter])):
			original_images_grouped_by_class_averages[class_avg_img_iter][img_iter] = rot_shift2D(original_images_grouped_by_class_averages[class_avg_img_iter][img_iter], *this_class_average_params[img_iter][:4])
			list_used_for_sending_parameters[list_count].extend(this_class_average_params[img_iter][:4] + [original_images_id_list_associated_with_this_class[img_iter]])    
			list_count += 1

	mpi_barrier(MPI_COMM_WORLD)	

### end: 2. align original images to class averages from isac using well chosen range parameters (class averages from isac have original image numbers)
########################################################################################################################

########################################################################################################################
### start: 3. Output alignment parameters and members for each class average

	twoD_params_info_list = wrap_mpi_gatherv(list_used_for_sending_parameters, main_node, comm)
	make_a_list_to_send_to_master__original_images_id_list_associated_with_this_class = wrap_mpi_gatherv(make_a_list_to_send_to_master__original_images_id_list_associated_with_this_class, main_node, comm)

	if myid == main_node:
		write_text_row(twoD_params_info_list, os.path.join(isac_post_processing_output_directory,"twoD_params_info_list.txt"))
		write_text_row(make_a_list_to_send_to_master__original_images_id_list_associated_with_this_class, os.path.join(isac_post_processing_output_directory,"original_image_number_in_each_class_average.txt"))

### end: 3. Output alignment parameters
########################################################################################################################

########################################################################################################################
### start: 4. Output the new averages using the newly aligned images

	mask = model_circle(radi, nx, nx)

	new_class_average_images = [EMData() for i in range(class_averages_nima) ]
	for class_avg_img_iter in xrange(class_averages_nima):
		new_class_average_images[class_avg_img_iter], _ = avgvar(original_images_grouped_by_class_averages[class_avg_img_iter], mode='b', interp='linear')

		p = Util.infomask(new_class_average_images[class_avg_img_iter], mask, False)
		new_class_average_images[class_avg_img_iter] -= p[0]
		p = Util.infomask(new_class_average_images[class_avg_img_iter], mask, True)
		new_class_average_images[class_avg_img_iter] /= p[1]

		new_class_average_images[class_avg_img_iter].set_attr_dict(class_average_images[class_avg_img_iter].get_attr_dict())
		new_class_average_images[class_avg_img_iter].write_image(size_adjusted_class_averages_file_name[:-4] + "_%04d"%(image_start + class_avg_img_iter)   +".hdf",0)

### end: 4. Output the new averages using the newly aligned images
########################################################################################################################

	mpi_finalize()

if __name__=="__main__":
	main(sys.argv[1:])


