"""
Plotting stuff for CTF Refinement with error assessment in SPHIRE

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
import sys
import os
import matplotlib as mpl

mpl.use("Agg")
from matplotlib import pyplot

from tqdm import tqdm
import sp_ctf_refine_io
import numpy as np


def create_particles_plot_map(
		stack_file_path,
		indices,
		values=None,
		vmin=None,
		vmax=None,
		colmap="viridis",
		title=None,
		colorlabel=None,
		plotpos=(),
):
	"""
	Plot the particle as dots where the color correspondents to the defocus value.
	:param stack_file_path:
	:param indices:
	:param value:
	:return: figure and axis of matplotlib
	"""
	x_coordinates = []
	y_coordinates = []
	value_for_coloring = []
	for value_index, particle_index in enumerate(indices):
		particle = sp_ctf_refine_io.read_particle(
			stack_file_path, particle_index, header_only=True
		)

		particle_coord = particle.get_attr("ptcl_source_coord")
		x_coordinates.append(particle_coord[0])
		y_coordinates.append(particle_coord[1])
		if values is None:
			particle_ctf = particle.get_attr("ctf")
			value_for_coloring.append(particle_ctf.defocus)
		else:
			val_for_col = values[value_index]
			value_for_coloring.append(val_for_col)

	axis = pyplot.subplot(*plotpos)
	scatter = axis.scatter(
		x_coordinates,
		y_coordinates,
		c=value_for_coloring,
		vmin=vmin,
		vmax=vmax,
		marker="o",
		cmap=colmap,
		s=50,
		lw=0,
	)
	pyplot.xlabel("Image dimension x")
	pyplot.ylabel("Image dimension y")
	cbar = pyplot.colorbar(scatter)
	if colorlabel:
		cbar.set_label(colorlabel, rotation=270)
	if title:
		axis.set_title(title, fontsize=7)
	if vmin != None or vmax != None:
		ticks = [t.get_text() for t in cbar.ax.get_yticklabels()]
		ticks[0] = "<" + ticks[0]
		ticks[-1] = ">" + ticks[-1]
		cbar.ax.set_yticklabels(ticks)

	pyplot.setp(axis.get_xticklabels(), rotation=30, horizontalalignment="right")
	return axis


def create_and_save_particle_plots(
		path_output_img,
		stack_file_path,
		refinement_results_per_micrograph,
		min_max_error,
		min_max_ratio,
):
	"""
	Create the plot with defocus map, error map and significance map.
	:param path_output_img: Path for output image file
	:param refinement_results_per_micrograph:
	:param min_max_error:
	:param min_max_ratio:
	:return:
	"""
	mpl.rcParams["figure.dpi"] = 200
	mpl.rcParams.update({"font.size": 7})
	all_errors = []
	with tqdm(total=len(refinement_results_per_micrograph), file=sys.stdout) as pbar:
		for mic_name in refinement_results_per_micrograph:
			particle_indices = refinement_results_per_micrograph[mic_name]["indices"]
			particle_error = refinement_results_per_micrograph[mic_name]["error"]
			all_errors.extend(particle_error)
			particle_defocus = refinement_results_per_micrograph[mic_name]["defocus"]
			particle_drratio = refinement_results_per_micrograph[mic_name]["drratio"]
			old_defocus = refinement_results_per_micrograph[mic_name]["diff"][0] + \
						  refinement_results_per_micrograph[mic_name]["defocus"][0]
			# First plot
			out_img_name = os.path.basename(mic_name) + ".png"
			out_img_name = os.path.join(path_output_img, out_img_name)
			create_particles_plot_map(
				stack_file_path=stack_file_path,
				indices=particle_indices,
				values=particle_defocus,
				title="Defocus map",
				plotpos=(2, 2, 1),
			)

			# Second plot
			create_particles_plot_map(
				stack_file_path=stack_file_path,
				indices=particle_indices,
				values=particle_error,
				# vmin=0,
				# vmax=0.08,
				colmap="coolwarm",
				title="Error map",
				plotpos=(2, 2, 2),
			)

			# Third plot
			create_particles_plot_map(
				stack_file_path=stack_file_path,
				indices=particle_indices,
				values=particle_drratio,
				# vmin=min_max_ratio[0],
				# vmax=min_max_ratio[1],
				colmap="YlOrBr",
				title="Significance map",
				plotpos=(2, 2, 3),
			)

			axis = pyplot.subplot(2, 2, 4)

			from numpy.lib.function_base import _hist_bin_fd
			width = int(np.round(_hist_bin_fd(np.array(particle_defocus))))

			width = max(10, width)
			n, _, _ = axis.hist(particle_defocus, bins=width, facecolor='green', alpha=0.75)
			pyplot.axvline(old_defocus, color='k', linestyle='dashed', linewidth=1)
			pyplot.text(old_defocus, max(n) / 2, "Old defocus", rotation=90,
						verticalalignment='center')
			pyplot.xlabel("Defocus")
			pyplot.ylabel("Frequency")

			pyplot.tight_layout()
			pyplot.savefig(out_img_name)
			pyplot.close()
			pbar.update(1)

		# PLOT CDF OF ERROR:
		out_img_name = os.path.join(
			os.path.dirname(os.path.dirname(path_output_img)), "error_cdf.png"
		)
		pyplot.hist(
			all_errors,
			bins=np.arange(0, np.max(all_errors) - 0.001, 0.001),
			normed=True,
			cumulative=True,
			label="CDF",
			histtype="step",
			color="k",
		)
		pyplot.xlabel("Error")
		pyplot.ylabel("CDF")
		pyplot.grid(b=True, linestyle="--")
		pyplot.tight_layout()
		pyplot.savefig(out_img_name)
		pyplot.close()
