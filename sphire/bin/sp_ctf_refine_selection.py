#! /usr/bin/env python
"""
Create substacks based on errors estimated with CTF Refinement in SPHIRE

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
import argparse
import numpy as np
import EMAN2db
import EMAN2
import sp_ctf_refine_io
import sp_global_def


def setup_argparser():
	argparser = argparse.ArgumentParser(
		description="Error assessment for CTF refinement",
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)

	argparser.add_argument("-r", "--resultsfile", help="Path to your results file")
	argparser.add_argument(
		"-m", "--mode", default="UPDATE", choices=["EXTRACT", "UPDATE"], help="Mode"
	)
	argparser.add_argument(
		"-f", "--field", default="ERROR", choices=["ERROR", "DR_RATIO"], help="Mode"
	)
	argparser.add_argument(
		"-lt",
		"--lower_then",
		type=float,
		default=0.1,
		help="Select particles with field values lower than this threshold",
	)
	argparser.add_argument(
		"-gt",
		"--greater_then",
		type=float,
		default=0,
		help="Select particles with field values greater than this threshold",
	)
	argparser.add_argument("-o", "--output", required=True, help="Path output bdb stack")
	argparser.add_argument(
		"-s", "--stack", required=True, help="Path to original bdb stack"
	)
	return argparser


def _main_():
	argparser = setup_argparser()

	args = argparser.parse_args()

	path_to_resultsfile = args.resultsfile
	mode = args.mode
	lower_then = args.lower_then
	greater_then = args.greater_then
	path_output = args.output
	path_stack = args.stack
	field = args.field
	field_index_error = 1
	field_index_dr_ratio = 3
	if field == "ERROR":
		field = field_index_error
	elif field == "DR_RATIO":
		field = field_index_dr_ratio

	# Read in error file
	results = np.loadtxt(path_to_resultsfile, delimiter=",")

	# Identify relevant particles
	relevant_selection = [
		a and b
		for a, b in zip(
			results[:, field] <= lower_then, results[:, field] >= greater_then
		)
	]

	# Write stack
	number_of_particles = EMAN2.EMUtil.get_image_count(path_stack)
	local_bdb_stack = EMAN2db.db_open_dict(path_output)
	num_particles_relevant = 0
	for particle_index in range(number_of_particles):

		particle = sp_ctf_refine_io.read_particle(
			path_stack, particle_index, header_only=True
		)
		particle_header = particle.get_attr_dict()

		if relevant_selection[particle_index]:
			pctf = particle_header["ctf"]
			pctf.defocus = results[particle_index, 2]
			particle_header["ctf"] = pctf
			num_particles_relevant = num_particles_relevant + 1

		if mode == "UPDATE":
			local_bdb_stack[particle_index] = particle_header
		elif mode == "EXTRACT" and relevant_selection[particle_index]:
			local_bdb_stack[num_particles_relevant - 1] = particle_header

	EMAN2db.db_close_dict(local_bdb_stack)
	sp_global_def.sxprint("Particles updated/extracted", num_particles_relevant)


if __name__ == "__main__":
	sp_global_def.print_timestamp("Start")
	_main_()
	sp_global_def.print_timestamp("Finish")
