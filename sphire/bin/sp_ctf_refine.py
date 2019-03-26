#! /usr/bin/env python
"""
CTF Refinement with error assessment for SPHIRE

#
# Author: Thorsten Wagner 02/18/2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Copyright (C) 2019 Max planck institute for molecular physiology, Dortmund
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
"""
# pylint: disable=C0330
import argparse
import sp_statistics
import time
import multiprocessing
import os
import itertools
import numpy as np
import sys
import copy

from tqdm import tqdm
from scipy import ndimage

import sp_ctf_refine_plotting
import sp_ctf_refine_io
import sp_filter as fltr
import EMAN2
import sp_global_def
from sp_projection import prgl


def create_ctf_list(current_ctf, def_search_range, def_step_size):
	"""
	Creates a set of CTFs with different defoci based on a defocus search range arround the
	defocus of the current ctf estimate.

	:param current_ctf: Current ctf
	:param def_search_range: Defocus search range
	:param def_step_size: Defocus step size
	:return: Set of CTFs
	"""
	current_defocus = current_ctf.defocus
	ctfs = []
	lower_defocus = current_defocus - def_search_range
	if lower_defocus < 0:
		lower_defocus = 0

	upper_defocus = current_defocus + def_search_range + def_step_size

	for defocus in np.arange(lower_defocus, upper_defocus, def_step_size):
		ctf_parameter = EMAN2.EMAN2Ctf()
		ctf_parameter.from_dict(current_ctf.to_dict())
		ctf_parameter.defocus = defocus
		ctfs.append(ctf_parameter)

	return ctfs


def calc_2d_projection(
		particle_volume, projection_angle_and_shift, interpolation=1, return_real=True
):
	"""
	:param particle_volume: input volume
	:param projection_angle_and_shift: input parameters given as a list [phi, theta, psi, s2x, s2y],
	projection in calculated using the three Eulerian angles and then shifted by sx,sy
	:param interpolation: interpolation_method,1 = triliniear
	:param return_real: True - return real; False - return FT of a projection.
	:return: 2D Projection
	"""
	projection = prgl(
		particle_volume, projection_angle_and_shift, interpolation, return_real
	)
	return projection


def calc_similarity_real(particle_image, projection_ctf_applied_2d, mask=None):
	"""
	Calculates the similarity based on projections
	:param particle_image: Experimental 2D particle image
	:param projection_ctf_applied_2d:
	:param mask: Mask for ccc
	:return: cross correlation coefficient
	"""
	ccc = sp_statistics.ccc(particle_image, projection_ctf_applied_2d, mask)
	return ccc


def calc_similarity_complex(particle_image, projection_ctf_applied_2d):
	"""
	Calculates the similarty based of fourier shell correlation
	:param particle_image: Experimental 2D particle image
	:param projection_ctf_applied_2d: CTF filtered 2D Projection
	:return: fourier shell correlation
	"""
	fsc = sp_statistics.fsc(particle_image, projection_ctf_applied_2d)
	data = np.array(fsc[1])
	return np.average(data[2:50])


def apply_ctf(projection_2d, ctf):
	"""
	Applies a CTF to a projection
	:param projection_2d: 2D Projection
	:param ctf: CTF
	:return: CTF filtered projection
	"""
	projection_filtered = fltr.filt_ctf(projection_2d, ctf)
	return projection_filtered


def create_mask(size, index_particle_pixels, fraction_size=0.25):
	mask = np.zeros(size)

	rand_pixel_selection = np.random.choice(
		range(len(index_particle_pixels[0])),
		size=int(len(index_particle_pixels[0]) * fraction_size),
		replace=False,
	)
	rand_pixel_selection = (
		index_particle_pixels[0][rand_pixel_selection],
		index_particle_pixels[1][rand_pixel_selection],
	)
	mask[rand_pixel_selection] = 1.0
	mask = EMAN2.EMNumPy.numpy2em(mask)
	return mask


def find_index_best_projection(particle_image, reprojections, mask):
	"""
	Compares the particle images with many reprojections to find the one with highest similarity.
	:param particle_image: Particle image
	:param reprojections: Projections with ctfs applied
	:param ctf_set: List of CTFSs which were used to create the ctfs corrupted CTFs.
	:param mask: Mask used during calculating of the CCC
	:return: Index in projections_ctf of best projection
	"""
	similarity_values = [
		calc_similarity_real(particle_image, projection_ctf_applied_2d=pro, mask=mask)
		for pro in reprojections
	]

	index_max = np.argmax(similarity_values)
	return index_max


def optimize(particle_image, projection_2d, ctf_list, masks):
	# #########
	# Calculate CTF corrupted projections
	# #########
	projections_ctf_applied_fft = [apply_ctf(projection_2d, ctf) for ctf in ctf_list]
	projections_ctf = [proj.do_ift() for proj in projections_ctf_applied_fft]

	# #########
	# Calculate the best ctf for each mask
	# #########
	index_best_projections = [
		find_index_best_projection(particle_image, projections_ctf, mask)
		for mask in masks
	]
	best_ctf_no_boostrap = ctf_list[index_best_projections[-1]]

	# #########
	# Estimate error and difference/error ratio
	# #########
	best_ctfs_bootstrap = [ctf_list[i] for i in index_best_projections[:-1]]

	return best_ctfs_bootstrap, best_ctf_no_boostrap


def refine_defocus_with_error_est(
		particle_volume,
		particle_image,
		projection_angle_and_shift,
		current_ctf,
		def_search_range,
		def_step_size,
):
	"""
	Runs defocus refinemenet for a single particle. It will create a reprojection of the particle
	given the known projection angle and shift parameters. It will use this reprojection to
	create multiple version of the projection each filtered with a different ctf. The CTFs used have
	different defoci in range given by the current defocus +- defocus search range.

	:param particle_volume: Particle volume
	:param particle_image: 2D particle image
	:param projection_angle_and_shift: Projection parameters
	:param current_ctf: Estimated CTF for the 2D particle image (micrograph)
	:param def_search_range: Defocus search range in micrometer
	:param def_step_size: Defocus step size in micrometer
	:return: CTF with highest correlation, the estimated
	error (standard deviation) and significance value
	"""

	# #########
	# Calculate projection
	# #########
	projection_2d = calc_2d_projection(
		particle_volume, projection_angle_and_shift, return_real=False
	)

	# #########
	# Calculate particle mask
	# #########
	mask2d = np.copy(projection_2d.do_ift().get_2dview())

	mask2d = ndimage.gaussian_filter(mask2d, 3)
	object_pixel = np.where(mask2d > 0.05)

	# #########
	# Calculate all masks for bootstrapping
	# #########
	num_bootstrap_runs = 15

	# Create maks for bootstrapping
	masks = [
		create_mask(particle_image.get_2dview().shape, object_pixel, fraction_size=0.25)
		for _ in range(num_bootstrap_runs)
	]
	# Append one mask that use all of the object_pixels to calculate the final defocus
	masks.append(
		create_mask(particle_image.get_2dview().shape, object_pixel, fraction_size=1.0)
	)

	# ########
	# Try to find the CTF region of interest, by unique ctf estimates of the bootstrapping
	# ########
	ctf_list_corse = create_ctf_list(current_ctf, def_search_range, def_step_size * 10)
	best_ctfs_bootstrap, best_ctf_no_boostrap = optimize(particle_image, projection_2d,
														 ctf_list_corse, masks)
	defocus_bootstrap = [ctf.defocus for ctf in best_ctfs_bootstrap]
	unique_defocuses = set(defocus_bootstrap)
	ctf_list_fine = []
	current_ctf_cpy = copy.copy(current_ctf)
	for unique_def in unique_defocuses:
		current_ctf_cpy.defocus = unique_def
		ctf_list_fine.extend(create_ctf_list(current_ctf_cpy, def_step_size * 10, def_step_size))

	# ########
	# Final optimization
	# ########
	best_ctfs_bootstrap, best_ctf_no_boostrap = optimize(particle_image, projection_2d,
														 ctf_list_fine, masks)

	# #########
	# Estimate error and difference/error ratiopyt
	# #########

	defocus_bootstrap = [ctf.defocus for ctf in best_ctfs_bootstrap]
	error_sd = np.std(defocus_bootstrap)
	drratio = np.abs(best_ctf_no_boostrap.defocus - current_ctf.defocus) / (
			error_sd + 0.0005
	)

	# #########
	# Cleaning!
	# #########
	del projection_2d
	del ctf_list_corse

	return best_ctf_no_boostrap, error_sd, drratio


# Python2 multiprocessing trick. Should replaced with change to python 3.6.
STACK_FILE_PATH = None
PROJECTION_PARAMETERS = None
VOLUME1 = None
VOLUME2 = None
DEFOCUS_SEARCH_RANGE = None
DEFOCUS_STEP_SIZE = None
MASK_VOLUME = None
RESOLUTION = None
PIXEL_SIZE = None
HALF_MAP_ASSIGNEMENTS = None


def add_proj_error(particle_projection_params, index, error):
	"""
	Just for testing purposes. Introduces errors to meridien parameters
	:param particle_projection_params:
	:param index:
	:param error:
	:return:
	"""

	new = particle_projection_params[index] + np.random.random() * error
	if index in (0, 2):
		if new < 0:
			new = 0
		elif new > 360:
			new = 360
	elif index == 1:
		if new < 0:
			new = 0
		elif new > 180:
			new = 180

	return new


def refine_set(particle_indices):
	"""
	Runs the defocus refinement for a set of particles.
	:param particle_indices: Indicies of the particles inside the stack
	:return: List of refined CTFs.
	"""
	try:
		refined_ctfs = []
		particle_micrographs = {}
		for particle_index in particle_indices:
			particle = sp_ctf_refine_io.read_particle(STACK_FILE_PATH, particle_index)
			# particle.write_image("/home/twagner/temp/refine/"+str(particle_index)+"particle.hdf")
			if RESOLUTION is not None and PIXEL_SIZE is not None:
				particle = particle.process(
					"filter.lowpass.gauss",
					{"cutoff_freq": 1.0 / RESOLUTION, "apix": PIXEL_SIZE},
				)

				# In case of phase plate a high pass might be necessary
				# if PARTICLE_DIAMETER is not None
				# particle = particle.process("filter.highpass.gauss",{"cutoff_freq":1.0/160})

			particle_projection_params = PROJECTION_PARAMETERS[particle_index].tolist()
			# As i shift the volume and not the particle, multiply by -1
			particle_projection_params[3] = -1 * particle_projection_params[3]
			particle_projection_params[4] = -1 * particle_projection_params[4]

			# particle_projection_params[3] = add_proj_error(particle_projection_params, 3, 2)
			# particle_projection_params[4] = add_proj_error(particle_projection_params, 4, 2)
			# particle_projection_params[2] = add_proj_error(particle_projection_params, 2, 2)

			particle_ctf = particle.get_attr("ctf")
			volume = VOLUME1
			if HALF_MAP_ASSIGNEMENTS is not None:
				if HALF_MAP_ASSIGNEMENTS[particle_index] != 0:
					volume = VOLUME2

			best_ctf, error, drratio = refine_defocus_with_error_est(
				particle_volume=volume,
				particle_image=particle,
				projection_angle_and_shift=particle_projection_params,
				current_ctf=particle_ctf,
				def_search_range=DEFOCUS_SEARCH_RANGE,
				def_step_size=DEFOCUS_STEP_SIZE,
			)

			particle_micrograph_filename = particle.get_attr("ptcl_source_image")
			if particle_micrograph_filename in particle_micrographs:
				particle_micrographs[particle_micrograph_filename]["indices"].append(
					particle_index
				)
				particle_micrographs[particle_micrograph_filename]["defocus"].append(
					best_ctf.defocus
				)
				particle_micrographs[particle_micrograph_filename]["diff"].append(
					particle_ctf.defocus - best_ctf.defocus
				)
				particle_micrographs[particle_micrograph_filename]["error"].append(
					error
				)
				particle_micrographs[particle_micrograph_filename]["drratio"].append(
					drratio
				)
			else:
				particle_micrographs[particle_micrograph_filename] = {
					"indices": [particle_index],
					"diff": [particle_ctf.defocus - best_ctf.defocus],
					"error": [error],
					"defocus": [best_ctf.defocus],
					"drratio": [drratio],
				}

			refined_ctfs.append(best_ctf)

			del particle_projection_params
			del volume
			del particle_ctf
			del particle
	except Exception as err:

		import traceback

		sp_global_def.sxprint(
			"Exception happend for particles in range ",
			np.min(particle_indices),
			"-",
			np.max(particle_indices),
		)

		traceback.print_exc()

		raise err
	return refined_ctfs, particle_micrographs


def indices_to_chunks(stack_path, chunk_size, number_of_particles_to_read=None):
	"""
	Divides the particles into chunks
	:param stack_path: stack file path
	:param chunk_size: Size of each chunk
	:param number_of_particles_to_read: Number of particles that sould be read from the stack
	:return: Chunks of indices
	"""
	number_of_particles = number_of_particles_to_read
	if number_of_particles_to_read is None:
		number_of_particles = EMAN2.EMUtil.get_image_count(stack_path)
	particle_inidices = range(number_of_particles)
	particle_chunks = [
		particle_inidices[x: x + chunk_size]
		for x in range(0, number_of_particles, chunk_size)
	]
	return particle_chunks, number_of_particles


def get_half_map_assigments_per_particle(
		stack_path, chunk_path, number_of_particles_to_read=None
):
	"""
	Given the chunk file, it will create an array where the index corresponds to the particle index
	and the value to half map (0 or 1)
	:param stack_path: particle bdb stack path
	:param chunk_path: path to chunk file
	:param number_of_particles_to_read: Number of particles to read
	:return: half map assignments. Index corresponds to the particle index
	and the value to half map
	"""

	number_of_particles = number_of_particles_to_read
	if number_of_particles_to_read is None:
		number_of_particles = EMAN2.EMUtil.get_image_count(stack_path)

	chunk = np.genfromtxt(chunk_path).astype(np.int64)

	if chunk[0] == 1:
		half_map_assignments = np.zeros(number_of_particles).astype(np.int64)
		half_map_assignments[chunk[:number_of_particles]] = 1
	else:
		half_map_assignments = np.ones(number_of_particles).astype(np.int64)
		half_map_assignments[:number_of_particles] = 0

	return half_map_assignments


def print_progress(refinement_result):
	"""
	Prints the progress of the refinement
	:param refinement_result: Object of the asynchronos pool.map
	:return: None
	"""
	start_num_chunks = refinement_result._number_left
	chunks_in_progress = start_num_chunks

	with tqdm(total=start_num_chunks, file=sys.stdout) as pbar:
		while True:
			num_unprocessed_chunks = refinement_result._number_left
			chunks_done = chunks_in_progress - num_unprocessed_chunks
			if chunks_done != 0:
				chunks_in_progress = num_unprocessed_chunks
				pbar.update(chunks_done)
			time.sleep(0.5)

			if num_unprocessed_chunks == 0:
				break


def calc_statistics(micrograph_indices):
	"""
	:param micrograph_indices: Dictonary where the keys are the micrographname and the values
	are list of particle indices.
	:return: Structured array with micrograph statistics.
	"""
	mic_stats = np.empty(
		len(micrograph_indices),
		dtype=[
			("Micrograph", "S35"),
			("MEAN", "<f8"),
			("MEAN_ABS", "<f8"),
			("STD", "<f8"),
			("STD_ABS", "<f8"),
			("p25", "<f8"),
			("p50", "<f8"),
			("p75", "<f8"),
			("MIN", "<f8"),
			("MAX", "<f8"),
			("MEAN_ERROR", "<f8"),
		],
	)

	for k, mic_name in enumerate(micrograph_indices.keys()):
		diff_abs_mean = np.mean(np.abs(micrograph_indices[mic_name]["diff"]))
		diff_abs_std = np.std(np.abs(micrograph_indices[mic_name]["diff"]))
		diff_mean = np.mean(micrograph_indices[mic_name]["diff"])
		diff_std = np.std(micrograph_indices[mic_name]["diff"])
		diff_p25 = np.percentile(micrograph_indices[mic_name]["diff"], q=25)
		diff_p50 = np.percentile(micrograph_indices[mic_name]["diff"], q=50)
		diff_p75 = np.percentile(micrograph_indices[mic_name]["diff"], q=75)
		diff_max = np.max(micrograph_indices[mic_name]["diff"])
		diff_min = np.min(micrograph_indices[mic_name]["diff"])
		mean_error = np.mean(micrograph_indices[mic_name]["error"])

		mic_stats[k] = (
			os.path.basename(mic_name),
			diff_mean,
			diff_abs_mean,
			diff_std,
			diff_abs_std,
			diff_p25,
			diff_p50,
			diff_p75,
			diff_min,
			diff_max,
			mean_error,
		)
	return mic_stats


def calculate_result_ranges(refinement_results_per_micrograph):
	"""
	Calculates the ranges of the error and
	:param refinement_results_per_micrograph:
	:return:
	"""
	all_errors = [
		refinement_results_per_micrograph[mic_name]["error"]
		for mic_name in refinement_results_per_micrograph
	]
	all_errors = list(itertools.chain(*all_errors))
	max_error = np.percentile(all_errors, 95)
	min_error = np.percentile(all_errors, 5)

	all_ratios = [
		refinement_results_per_micrograph[mic_name]["drratio"]
		for mic_name in refinement_results_per_micrograph
	]
	all_ratios = list(itertools.chain(*all_ratios))
	max_ratio = np.percentile(all_ratios, 95)
	min_ratio = np.percentile(all_ratios, 5)

	return (min_error, max_error), (min_ratio, max_ratio)


def get_refinement_results_matrix(refinement_results_per_micrograph):
	"""
	Builds a refinement results matrix. Column 1: Indices, columns 2: errors, column 3: new defocus,
	column 4: significance
	:param refinement_results_per_micrograph: Results per micrograph list.
	:return: Refinement results matrix
	"""

	all_errors = [
		refinement_results_per_micrograph[mic_name]["error"]
		for mic_name in refinement_results_per_micrograph
	]
	all_errors = list(itertools.chain(*all_errors))

	all_indices = [
		refinement_results_per_micrograph[mic_name]["indices"]
		for mic_name in refinement_results_per_micrograph
	]
	all_indices = list(itertools.chain(*all_indices))

	all_defocus = [
		refinement_results_per_micrograph[mic_name]["defocus"]
		for mic_name in refinement_results_per_micrograph
	]
	all_defocus = list(itertools.chain(*all_defocus))

	all_drratio = [
		refinement_results_per_micrograph[mic_name]["drratio"]
		for mic_name in refinement_results_per_micrograph
	]
	all_drratio = list(itertools.chain(*all_drratio))

	index_error_matrix = np.column_stack(
		(all_indices, all_errors, all_defocus, all_drratio)
	)
	return index_error_matrix


def merge_ctf_refinement_results(refinement_results):
	"""
	Merges the refinement result chunks and sorts them according to the micrograph
	:param refinement_results: Resulting object of multiprocessing map_async
	:return: Refined ctfs as list and all results per micrograph
	"""
	refined_ctfs_as_list = []
	refinement_results_per_micrograph = {}
	for chunk_refinement_result in refinement_results:
		chunk_ctfs = chunk_refinement_result[0]
		refined_ctfs_as_list.extend(chunk_ctfs)
		particle_indicies_per_mic_dict = chunk_refinement_result[1]

		for mic_name in particle_indicies_per_mic_dict.keys():
			if mic_name not in refinement_results_per_micrograph.keys():
				refinement_results_per_micrograph[mic_name] = {
					"indices": [],
					"diff": [],
					"error": [],
					"defocus": [],
					"drratio": [],
				}

			refinement_results_per_micrograph[mic_name]["indices"].extend(
				particle_indicies_per_mic_dict[mic_name]["indices"]
			)

			refinement_results_per_micrograph[mic_name]["diff"].extend(
				particle_indicies_per_mic_dict[mic_name]["diff"]
			)

			refinement_results_per_micrograph[mic_name]["error"].extend(
				particle_indicies_per_mic_dict[mic_name]["error"]
			)

			refinement_results_per_micrograph[mic_name]["defocus"].extend(
				particle_indicies_per_mic_dict[mic_name]["defocus"]
			)

			refinement_results_per_micrograph[mic_name]["drratio"].extend(
				particle_indicies_per_mic_dict[mic_name]["drratio"]
			)
	return refined_ctfs_as_list, refinement_results_per_micrograph


def setup_argparser():
	example_txt_meridien = """Example:

    sxrefinectf.py bdb:/path/to/stack/input/STACKNAME bdb_file_name /path/to/output/folder/ /path/to/meridien/refinement/folder/
                  -m /pth/adp_mask/mask.hdf -r 0.2 -d 0.01 -res 2.2 -apix 1.07

    """

	example_txt_manual = """Example:

    sxrefinectf.py -v1 /path/to/volume/vol_0_unfil_XXX.hdf -v2 /path/to/volume/vol_0_unfil_XXX.hdf 
                -c /path/to/meridien/mainXXX/chunk_0_XXX.txt -s bdb:/path/to/stack/STACKNAME 
                -p /path/to/merdien/refinement/final_params_XXX.txt -o bdb:/your/output/path/bdb_file_name 
                -os /statistics/out/path/ -m /pth/adp_mask/mask.hdf -r 0.2 -d 0.01 -res 2.2 -apix 1.07

    """

	argparser = argparse.ArgumentParser(
		description="Refine your CTF",
		formatter_class=argparse.RawDescriptionHelpFormatter,
		add_help=False,
	)

	argparser.add_argument("inputstack", help="Path to your particle stack")

	argparser.add_argument("outputdir", help="Path to the output directory")

	argparser.add_argument(
		"-r", "--range", default=0.15, type=float, help="Defocus search range (in microns)"
	)

	argparser.add_argument(
		"-d", "--delta", default=0.0025, type=float, help="Defocus step size (in microns)"
	)

	argparser.add_argument(
		"-res", "--resolution", type=float, help="Nominal resolution (in angstrom)"
	)

	argparser.add_argument(
		"-apix", "--pixelsize", type=float, help="Pixel size (in angstrom)"
	)

	argparser.add_argument("-m", "--mask", help="Path to adaptive mask for the volume")

	argparser.add_argument(
		"-num",
		"--number_particles",
		type=int,
		help="Number of particles to process. Option is mainly used " "for debugging.",
	)

	child_argparser = argparse.ArgumentParser(
		description="Refine your CTF",
		formatter_class=argparse.RawDescriptionHelpFormatter,
		add_help=True,
	)

	subparsers = child_argparser.add_subparsers(
		help="You can either run it from meridien or specify everything manually"
	)
	parser_meridien = subparsers.add_parser(
		"meridien",
		help="Run it by specifing meridien folder",
		parents=[argparser],
		epilog=example_txt_meridien,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	parser_manual = subparsers.add_parser(
		"manual",
		help="Run it by specifing parameters manually",
		parents=[argparser],
		epilog=example_txt_manual,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)

	parser_manual.add_argument("volume", help="Path to your first half map")

	parser_manual.add_argument("params_path", help="Path to your params file")

	parser_manual.add_argument(
		"-v2",
		"--volume2",
		help="Path to your second half map. Only necessary if chunk file is available",
	)

	# TODO: Use both chunks might be necessay in the future.
	parser_manual.add_argument("-c", "--chunk", help="Path to one of the chunk files")

	parser_meridien.add_argument("meridien_path", help="Path to meridien refinement folder")

	return argparser, child_argparser


def _main_():
	_, child_argparser = setup_argparser()
	# These global variables are kind of ugly, but are necessary in python2. Will be removed in
	# python 3.
	global STACK_FILE_PATH
	global PROJECTION_PARAMETERS
	global VOLUME1
	global VOLUME2
	global DEFOCUS_SEARCH_RANGE
	global DEFOCUS_STEP_SIZE
	global MASK_VOLUME
	global RESOLUTION
	global PIXEL_SIZE
	global HALF_MAP_ASSIGNEMENTS

	args = child_argparser.parse_args()

	if "meridien" in sys.argv[1]:
		meridien_path = args.meridien_path
		files = sp_ctf_refine_io.read_meridien_data(meridien_path)
		volume1_file_path = files["first_halfmap"]
		volume2_file_path = files["second_halfmap"]
		chunk_file_path = files["chunk1"]
		params_file_path = files["final_params"]
	else:
		volume1_file_path = args.volume
		volume2_file_path = args.volume2
		params_file_path = args.params_path
		chunk_file_path = args.chunk
		if volume2_file_path is None and chunk_file_path is not None:
			sp_global_def.ERROR(
				"If chunk file is specified, you need to specify a second volume (-v2)"
			)

	STACK_FILE_PATH = args.inputstack  # "bdb:name"

	mask_file_path = args.mask

	DEFOCUS_SEARCH_RANGE = args.range

	DEFOCUS_STEP_SIZE = args.delta

	output_folder = args.outputdir

	if os.path.exists(output_folder):
		sp_global_def.ERROR("Output folder already exists. Stop execution.")
	else:
		os.makedirs(output_folder)
		sp_global_def.write_command(output_folder)

	output_virtual_stack_path = "bdb:" + os.path.join(output_folder, "ctf_refined")

	output_stats_path = os.path.join(output_folder, "statistics")

	number_of_particles_to_read = args.number_particles

	volume_nominal_resolution = args.resolution
	RESOLUTION = volume_nominal_resolution
	PIXEL_SIZE = args.pixelsize

	PROJECTION_PARAMETERS = sp_ctf_refine_io.read_meridien_params(params_file_path)

	num_cpu = multiprocessing.cpu_count() - 1

	volume1, volume2, MASK_VOLUME = sp_ctf_refine_io.read_volume(
		path_vol_1=volume1_file_path,
		path_vol_2=volume2_file_path,
		path_mask=mask_file_path,
		resolution=volume_nominal_resolution,
		pixel_size=PIXEL_SIZE,
	)
	VOLUME1 = volume1
	VOLUME2 = volume2

	particle_chunks, number_of_particles = indices_to_chunks(
		STACK_FILE_PATH,
		chunk_size=100,
		number_of_particles_to_read=number_of_particles_to_read,
	)

	if chunk_file_path:
		HALF_MAP_ASSIGNEMENTS = get_half_map_assigments_per_particle(
			stack_path=STACK_FILE_PATH,
			chunk_path=chunk_file_path,
			number_of_particles_to_read=number_of_particles_to_read,
		)

	sp_global_def.sxprint("####Start refinement####")
	start = time.time()

	# for chunk in particle_chunks:
	#   refine_set(chunk)

	#####################################################################
	# Python 2 workaround to get rid of a memory leak.
	# Pool eats up memory and only give it back after it returns.
	# So let it return more often by deviding chunks into chunks :-)
	#####################################################################
	num_chunks = len(particle_chunks)
	refinement_results = []

	with tqdm(total=num_chunks, file=sys.stdout) as pbar:
		for i in range(0, num_chunks, num_cpu):
			subset_chunk = particle_chunks[i: (i + num_cpu)]
			pool = multiprocessing.Pool(num_cpu)
			refinement_result = pool.map_async(refine_set, subset_chunk, chunksize=1)
			pool.close()
			# print_progress(refinement_result)
			pool.join()
			for res in refinement_result.get():
				refinement_results.append(res)
			pbar.update(len(subset_chunk))
	#####################################################################

	end = time.time()
	sp_global_def.sxprint("Time for ", number_of_particles, " Particles:", end - start)

	# Ouput results
	refined_ctfs_as_list, refinement_results_per_micrograph = merge_ctf_refinement_results(
		refinement_results=refinement_results
	)

	sp_ctf_refine_io.write_virtual_bdb_stack(
		output_stack_path=output_virtual_stack_path,
		origin_stack_path=STACK_FILE_PATH,
		refined_ctfs=refined_ctfs_as_list,
		number_of_particles=number_of_particles,
	)

	sp_global_def.sxprint("Write statistics...")
	# WRITE STATISTICS
	if not os.path.exists(output_stats_path):
		os.makedirs(output_stats_path)

	refinement_stats_per_micrograh = calc_statistics(refinement_results_per_micrograph)
	sp_ctf_refine_io.write_statistics(output_stats_path, refinement_stats_per_micrograh)

	# Estimate error Range
	min_max_error, min_max_ratio = calculate_result_ranges(
		refinement_results_per_micrograph
	)
	# Save particle plots
	sp_global_def.sxprint("Write images...")
	path_output_img = os.path.join(output_stats_path, "img/")
	if not os.path.exists(path_output_img):
		os.makedirs(path_output_img)
	sp_ctf_refine_plotting.create_and_save_particle_plots(
		path_output_img=path_output_img,
		stack_file_path=STACK_FILE_PATH,
		refinement_results_per_micrograph=refinement_results_per_micrograph,
		min_max_error=min_max_error,
		min_max_ratio=min_max_ratio,
	)

	sp_global_def.sxprint("Write other...")
	refinement_results_matrix = get_refinement_results_matrix(
		refinement_results_per_micrograph
	)

	particle_error_path = os.path.join(output_stats_path, "particle_results.txt")
	refinement_results_matrix = refinement_results_matrix[
		refinement_results_matrix[:, 0].argsort()
	]
	np.savetxt(
		particle_error_path,
		refinement_results_matrix,
		delimiter=",",
		fmt=["%d", "%f", "%f", "%f"],
	)

	sp_global_def.sxprint("Done")


if __name__ == "__main__":
	sp_global_def.BATCH = True
	sp_global_def.print_timestamp("Start")
	_main_()
	sp_global_def.print_timestamp("Finish")
