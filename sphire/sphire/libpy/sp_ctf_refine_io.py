
from __future__ import division
from past.utils import old_div

"""
IO stuff for CTF Refinement with error assessment in SPHIRE

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
"""
# pylint: disable=C0330
import EMAN2
import EMAN2db
import numpy
import os
from . import sp_global_def
from . import sp_projection
from . import sp_utilities


def read_meridien_data(meridien_path):
    """
    Reads all necessary paths from the meridien folder.
    :param meridien_path: Path to meridien folder
    :return: Dictoniary with paths
    """
    # TODO: READ MASK AND STACK FROM TRACKER
    files = [
        name
        for name in os.listdir(meridien_path)
        if not os.path.isdir(os.path.join(meridien_path, name))
    ]
    final_params_filename = None
    best_iteration_number = None
    for file in files:
        if file[:13] == "final_params_":
            final_params_filename = file
            best_iteration_number = file[13:-4]

    best_iteration_dir = "main" + best_iteration_number

    final_params_path = os.path.join(meridien_path, final_params_filename)
    first_halfmap_filename = "vol_0_unfil_" + best_iteration_number + ".hdf"
    first_halfmap_path = os.path.join(meridien_path, first_halfmap_filename)
    second_halfmap_filename = "vol_1_unfil_" + best_iteration_number + ".hdf"
    second_halfmap_path = os.path.join(meridien_path, second_halfmap_filename)
    particle_chunk0_filename = "chunk_0_" + best_iteration_number + ".txt"
    particle_chunk1_filename = "chunk_1_" + best_iteration_number + ".txt"
    particle_chunk0_path = os.path.join(
        meridien_path, best_iteration_dir, particle_chunk0_filename
    )
    particle_chunk1_path = os.path.join(
        meridien_path, best_iteration_dir, particle_chunk1_filename
    )

    return {
        "final_params": final_params_path,
        "first_halfmap": first_halfmap_path,
        "second_halfmap": second_halfmap_path,
        "chunk0": particle_chunk0_path,
        "chunk1": particle_chunk1_path,
    }


def write_virtual_bdb_stack(
    output_stack_path, origin_stack_path, refined_ctfs, number_of_particles=None
):
    """
    Write virtual bdb stack to disk
    :param output_stack_path: Output path for bdb database
    :param origin_stack_path: Path to original stack
    :param refined_ctfs: List of the refined CTFs
    :param number_of_particles: Number of particles to write.
    :return: None
    """
    sp_global_def.sxprint("Write results to virtual stack...")

    if number_of_particles is None:
        number_of_particles = EMAN2.EMUtil.get_image_count(origin_stack_path)

    local_bdb_stack = EMAN2db.db_open_dict(output_stack_path)
    particle_headers = sp_utilities.make_v_stack_header(
        origin_stack_path, output_stack_path
    )

    for particle_index in range(number_of_particles):
        particle_header = particle_headers[particle_index]
        particle_header["ctf"] = refined_ctfs[particle_index]
        local_bdb_stack[particle_index] = particle_header

    EMAN2db.db_close_dict(local_bdb_stack)
    sp_global_def.sxprint("Write results to virtual stack done")


def read_meridien_params(path):
    """
    Reads the refinement parameters from meridien.
    :param path: Path to params.txt
    :return: 2D numpy array with refinement parameters.
    """
    return numpy.atleast_2d(numpy.genfromtxt(path))


def read_particle(path, index, header_only=False):
    """
    Reads single particle from stack.
    :param path: Stack from path
    :param index: Index of the particle in the stack
    :param header_only: If true, only the header information will be read
    :return: 2D particle image
    """
    particle_image = EMAN2.EMData()
    particle_image.read_image(path, index, header_only)
    return particle_image


def prepare_volume(volume_path, mask=None, resolution=None, pixel_size=None):
    """
    Prepares the volume for reprojections
    :param volume_path: Path to volume file
    :param mask: Particle mask
    :param resolution: Resolution of the current reconstruction.
    :return:
    """
    vol = sp_utilities.get_im(volume_path)
    if resolution:
        vol = vol.process(
            "filter.lowpass.gauss",
            {"cutoff_freq": old_div(1.0, resolution), "apix": pixel_size},
        )
    if mask:
        vol = vol * mask
    vol = sp_projection.prep_vol(vol, interpolation_method=1)
    return vol


def read_volume(
    path_vol_1, path_vol_2, path_mask=None, resolution=None, pixel_size=None
):
    """
    Reads the 3D density map
    :param path_vol_1: Path to the density map
    :param path_mask: If provided, the volume is masked.
    :return: Volume for density map
    """
    print("Read volumes, masking them and prepare them for reprojections")
    mask_vol = None
    vol2 = None
    if path_mask is not None:
        mask_vol = sp_utilities.get_im(path_mask)

    vol1 = prepare_volume(path_vol_1, mask_vol, resolution, pixel_size)
    if path_vol_2:
        vol2 = prepare_volume(path_vol_2, mask_vol, resolution, pixel_size)

    if mask_vol:
        mask_vol = sp_projection.prep_vol(mask_vol, interpolation_method=1)

    return vol1, vol2, mask_vol


def write_statistics(output_stats_path, mic_stats):
    """
    Write the refinement statistics to file.
    :param output_stats_path: Directory to save the statistics
    :param mic_stats: Micrograph statistics as structured array
    :return: None
    """

    diff_abs_mean_list = []
    diff_abs_std_list = []
    diff_mean_list = []
    diff_std_list = []
    diff_p25_list = []
    diff_p50_list = []
    diff_p75_list = []
    diff_max_list = []
    diff_min_list = []
    mean_error_list = []

    for mic_row in mic_stats:
        diff_abs_mean_list.append(mic_row[2])
        diff_abs_std_list.append(mic_row[4])
        diff_mean_list.append(mic_row[1])
        diff_std_list.append(mic_row[3])
        diff_p25_list.append(mic_row[5])
        diff_p50_list.append(mic_row[6])
        diff_p75_list.append(mic_row[7])
        diff_max_list.append(mic_row[9])
        diff_min_list.append(mic_row[8])
        mean_error_list.append(mic_row[10])

    header = (
        "Statistics about defocus changes (old-new) in micrometer \n "
        "{0: ^35s} {1: ^10s} {2: ^10s} {3: ^10s} {4: ^10s} {5: ^10s} "
        "{6: ^10s} {7: ^10s} {8: ^10s} {9: ^10s} {10: ^10s}".format(
            "MICROGRAPH",
            "MEAN",
            "MEAN_ABS",
            "SD",
            "SD_ABS",
            "25q",
            "50q",
            "75q",
            "MIN",
            "MAX",
            "ERROR",
        )
    )

    footer = (
        "Average over all micrographs \n "
        "{0: ^35s} {1: ^10.5f} {2: ^10.5f} {3: ^10.5f} {4: ^10.5f} {5: ^10.5f} "
        "{6: ^10.5f} {7: ^10.5f} {8: ^10.5f} {9: ^10.5f} {10: ^10.5f}".format(
            "AVERAGE",
            numpy.mean(diff_mean_list),
            numpy.mean(diff_abs_mean_list),
            numpy.mean(diff_std_list),
            numpy.mean(diff_abs_std_list),
            numpy.mean(diff_p25_list),
            numpy.mean(diff_p50_list),
            numpy.mean(diff_p75_list),
            numpy.mean(diff_min_list),
            numpy.mean(diff_max_list),
            numpy.mean(mean_error_list),
        )
    )

    stats_path = os.path.join(output_stats_path, "stats.txt")
    numpy.savetxt(
        stats_path,
        mic_stats,
        fmt=[
            "% 35s",
            "% 10.5f",
            "% 10.5f",
            "% 10.5f",
            "% 10.5f",
            "% 10.5f",
            "% 10.5f",
            "% 10.5f",
            "% 10.5f",
            "% 10.5f",
            "% 10.5f",
        ],
        header=header,
        footer=footer,
    )  # ,fmt='% 10.5f')
