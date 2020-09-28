#! /usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div

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

import argparse
import os
from ..libpy import sp_global_def
from ..libpy import sp_utilities
from ..libpy import sp_morphology
from ..libpy import sp_filter


class NotSmallerZeroAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values < 0:
            parser.error("Minimum value for {0} is 0".format(option_string))
        setattr(namespace, self.dest, values)


def parse_command_line():
    """
	Parse the command line.

	Arguments:
	None:

	Returns:
	Parsed arguments object
	"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("input_volume", type=str, help="Input volume for masking")
    parser.add_argument(
        "output_dir", type=str, help="Output directory containing the output files"
    )
    parser.add_argument(
        "--prefix", type=str, default="sp_mask", help="Prefix for the produced files"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Override existing files with the same output_dir/prefix combination.",
    )

    group_filter = parser.add_argument_group(
        title="Filter",
        description="Filter options. The filter is applied prior to masking.",
    )
    group_filter.add_argument(
        "--low_pass_filter_resolution",
        "--fl",
        type=float,
        default=None,
        action=NotSmallerZeroAction,
        help="Target low-pass filter resolution in angstrom.",
    )
    group_filter.add_argument(
        "--low_pass_filter_falloff",
        "--aa",
        type=float,
        default=0.01,
        action=NotSmallerZeroAction,
        help="Low pass filter falloff in absolute frequencies. Used for filtering the volume.",
    )
    group_filter.add_argument(
        "--pixel_size",
        "--apix",
        type=float,
        default=1.0,
        help="Pixel size of the input volume. Used for filtering the volume.",
    )

    group_mask = parser.add_argument_group(
        title="Mask",
        description="Mask creation related options. If edge_width is 0, a binary mask is produced.",
    )
    group_threshold = group_mask.add_mutually_exclusive_group(required=True)
    group_threshold.add_argument(
        "--threshold",
        "--thr",
        type=float,
        default=None,
        help="Binarization threshold optained for example from chimera.",
    )
    group_threshold.add_argument(
        "--nsigma",
        "--ns",
        type=float,
        default=None,
        action=NotSmallerZeroAction,
        help="Number of times of sigma of the input volume to calculate the binarization threshold.",
    )
    group_threshold.add_argument(
        "--mol_mass",
        "--mm",
        type=float,
        default=None,
        action=NotSmallerZeroAction,
        help="The estimated molecular mass of the target particle in kD.",
    )

    group_mask.add_argument(
        "--ndilation",
        "--nd",
        type=int,
        default=3,
        action=NotSmallerZeroAction,
        help="Number of iterations to dilate the binarized volume. Each iteration adds about two pixels to the mask.",
    )
    group_mask.add_argument(
        "--nerosion",
        "--ne",
        type=int,
        default=0,
        action=NotSmallerZeroAction,
        help="Number of iterations to erode the binarized volume. Each iteration removes about two pixels from the mask.",
    )
    group_mask.add_argument(
        "--edge_width",
        "--ew",
        type=int,
        default=5,
        action=NotSmallerZeroAction,
        help="Width of the soft edge of the mask in pixels",
    )
    group_mask.add_argument(
        "--edge_type",
        "--et",
        type=str,
        default="cosine",
        choices=("cosine", "gaussian"),
        help="Type of the soft edge.",
    )
    group_mask.add_argument(
        "--do_old",
        "--da",
        action="store_true",
        help="Restore the old masking behaviour, which is a bit less smooth.",
    )
    group_mask.add_argument(
        "--allow_disconnected",
        "--ad",
        action="store_true",
        default=False,
        help="Allow disconnected regions in the mask.",
    )
    group_mask.add_argument(
        "--fill_mask",
        "--fill",
        action="store_true",
        default=False,
        help="Fill empty regions in the mask.",
    )

    group_second_mask = parser.add_argument_group(
        title="Second mask",
        description="Second mask creation related options. If edge_width is 0, a binary mask is produced. The second mask will be multiplied with the volume mask. If ",
    )
    group_second_mask_mut = group_second_mask.add_mutually_exclusive_group()
    group_second_mask_mut.add_argument(
        "--second_mask",
        "--sm",
        type=str,
        default=None,
        help="File path of the optional second mask.",
    )
    group_second_mask_mut.add_argument(
        "--second_mask_shape",
        "--sms",
        type=str,
        default=None,
        choices=("cylinder", "sphere", "cube"),
        help="Binarization threshold optained for example from chimera.",
    )

    group_threshold = group_second_mask.add_mutually_exclusive_group()
    group_threshold.add_argument(
        "--s_threshold",
        "--sthr",
        type=float,
        default=None,
        help="Second mask mask_path: Binarization threshold optained for example from chimera.",
    )
    group_threshold.add_argument(
        "--s_nsigma",
        "--sns",
        type=float,
        default=None,
        action=NotSmallerZeroAction,
        help="Second mask mask_path: Number of times of sigma of the input volume to calculate the binarization threshold.",
    )
    group_filter.add_argument(
        "--s_pixel_size",
        "--sapix",
        type=float,
        default=1.0,
        help="Second mask pixel_size: Pixel size of the input volume. Used for filtering the volume.",
    )
    group_threshold.add_argument(
        "--s_mol_mass",
        "--smm",
        type=float,
        default=None,
        action=NotSmallerZeroAction,
        help="Second mask mask_path: The estimated molecular mass of the target particle in kD.",
    )

    group_second_mask.add_argument(
        "--s_radius",
        "--sr",
        type=int,
        default=None,
        action=NotSmallerZeroAction,
        help="Second mask sphere/cylinder: Sphere/cylinder radius.",
    )
    group_second_mask.add_argument(
        "--s_nx",
        "--snx",
        type=int,
        default=None,
        action=NotSmallerZeroAction,
        help="Second mask cube/sphere/cylinder: X Dimension of the mask. The final mask volume will be zero padded to the same dimensions as the input image.",
    )
    group_second_mask.add_argument(
        "--s_ny",
        "--sny",
        type=int,
        default=None,
        action=NotSmallerZeroAction,
        help="Second mask cube/sphere/cylinder: Y Dimension of the mask. The final mask volume will be zero padded to the same dimensions as the input image.",
    )
    group_second_mask.add_argument(
        "--s_nz",
        "--snz",
        type=int,
        default=None,
        action=NotSmallerZeroAction,
        help="Second mask cube/sphere/cylinder: Z Dimension of the mask. The final mask volume will be zero padded to the same dimensions as the input image.",
    )

    group_second_mask.add_argument(
        "--s_ndilation",
        "--snd",
        type=int,
        default=3,
        action=NotSmallerZeroAction,
        help="Second mask: Number of iterations to dilate the binarized volume. Each iteration adds about two pixels to the mask.",
    )
    group_second_mask.add_argument(
        "--s_nerosion",
        "--sne",
        type=int,
        default=0,
        action=NotSmallerZeroAction,
        help="Second mask: Number of iterations to erode the binarized volume. Each iteration removes about two pixels from the mask.",
    )
    group_second_mask.add_argument(
        "--s_edge_width",
        "--sew",
        type=int,
        default=5,
        action=NotSmallerZeroAction,
        help="Second mask: Width of the soft edge of the mask in pixels",
    )
    group_second_mask.add_argument(
        "--s_edge_type",
        "--set",
        type=str,
        default="cosine",
        choices=("cosine", "gaussian"),
        help="Type of the soft edge.",
    )
    group_mask.add_argument(
        "--s_do_old",
        "--sda",
        action="store_true",
        help="Second mask: Restore the old masking behaviour, which is a bit less smooth.",
    )
    group_second_mask.add_argument(
        "--s_allow_disconnected",
        "--sad",
        action="store_true",
        default=False,
        help="Second mask: Allow disconnected regions in the mask.",
    )
    group_mask.add_argument(
        "--s_fill_mask",
        "--sfm",
        action="store_true",
        default=False,
        help="Second mask: Fill empty regions in the mask.",
    )
    group_second_mask.add_argument(
        "--s_invert",
        "--sinv",
        action="store_true",
        default=False,
        help="Second mask: Exclude the second mask region instead of including it.",
    )

    return parser.parse_args()


def sanity_checks(command_args, input_vol):

    if os.path.isfile(
        os.path.join(command_args.output_dir, command_args.prefix + "_mask.hdf")
    ):
        if not command_args.overwrite:

            sp_global_def.ERROR(
                "Output mask already exists! Please provide the overwrite option if you want to overwrite the existing mask."
            )

    if (
        command_args.s_nx is not None
        and command_args.s_ny is None
        and command_args.s_ny is None
        and command_args.second_mask_shape
    ):
        command_args.s_ny = command_args.s_nx
        command_args.s_nz = command_args.s_nx
    elif (
        command_args.s_nx is not None
        and command_args.s_ny is not None
        and command_args.s_ny is not None
        and command_args.second_mask_shape
    ):
        pass
    elif command_args.second_mask_shape:
        sp_global_def.ERROR("You need to specify s_nx only or s_nx and s_ny and s_nz")

    if command_args.second_mask_shape in ("cylinder", "sphere"):
        nx = input_vol.get_xsize()
        ny = input_vol.get_ysize()
        nz = input_vol.get_zsize()
        if command_args.s_radius > old_div(nx, 2) or command_args.s_radius > old_div(
            ny, 2
        ):
            sp_global_def.ERROR(
                "Provided radius is larger than input image dimensions!"
            )

    if command_args.second_mask_shape is not None:
        nx = input_vol.get_xsize()
        ny = input_vol.get_ysize()
        nz = input_vol.get_zsize()
        if command_args.s_nx > nx or command_args.s_ny > ny or command_args.s_nz > nz:
            sp_global_def.ERROR(
                "Provided s_nx, s_ny, s_nz mask dimension is larger than input image dimensions!"
            )


def run():
    """
	Main function.

	Arguments:
	None

	Returns:
	None
	"""

    command_args = parse_command_line()

    # Import volume
    sp_global_def.sxprint("Import volume.")
    input_vol = sp_utilities.get_im(command_args.input_volume)

    # Sanity checks
    sanity_checks(command_args, input_vol)

    try:
        os.makedirs(command_args.output_dir)
    except OSError:
        sp_global_def.sxprint("Output directory already exists. No need to create it.")
    else:
        sp_global_def.sxprint("Created output directory.")
    sp_global_def.write_command(command_args.output_dir)
    output_prefix = os.path.join(command_args.output_dir, command_args.prefix)

    # Filter volume if specified
    if command_args.low_pass_filter_resolution is not None:
        sp_global_def.sxprint(
            "Filter volume to {0}A.".format(command_args.low_pass_filter_resolution)
        )
        input_vol = sp_filter.filt_tanl(
            input_vol,
            old_div(command_args.pixel_size, command_args.low_pass_filter_resolution),
            command_args.low_pass_filter_falloff,
        )
        input_vol.write_image(output_prefix + "_filtered_volume.hdf")
    else:
        sp_global_def.sxprint("Skip filter volume.")

    # Create a mask based on the filtered volume
    sp_global_def.sxprint("Create mask")
    density_threshold = -9999.0
    nsigma = 1.0
    if command_args.mol_mass:
        density_threshold = input_vol.find_3d_threshold(
            command_args.mol_mass, command_args.pixel_size
        )
        sp_global_def.sxprint(
            "Mask molecular mass translated into binary threshold: ", density_threshold
        )
    elif command_args.threshold:
        density_threshold = command_args.threshold
    elif command_args.nsigma:
        nsigma = command_args.nsigma
    else:
        assert False

    if command_args.edge_type == "cosine":
        mode = "C"
    elif command_args.edge_type == "gaussian":
        mode = "G"
    else:
        assert False

    mask_first = sp_morphology.adaptive_mask_scipy(
        input_vol,
        nsigma=nsigma,
        threshold=density_threshold,
        ndilation=command_args.ndilation,
        nerosion=command_args.nerosion,
        edge_width=command_args.edge_width,
        allow_disconnected=command_args.allow_disconnected,
        mode=mode,
        do_approx=command_args.do_old,
        do_fill=command_args.fill_mask,
        do_print=True,
    )

    # Create a second mask based on the filtered volume
    s_mask = None
    s_density_threshold = 1
    s_nsigma = 1.0
    if command_args.second_mask is not None:
        sp_global_def.sxprint("Prepare second mask")
        s_mask = sp_utilities.get_im(command_args.second_mask)
        s_density_threshold = -9999.0
        s_nsigma = 1.0
        if command_args.s_mol_mass:
            s_density_threshold = input_vol.find_3d_threshold(
                command_args.s_mol_mass, command_args.s_pixel_size
            )
            sp_global_def.sxprint(
                "Second mask molecular mass translated into binary threshold: ",
                s_density_threshold,
            )
        elif command_args.s_threshold:
            s_density_threshold = command_args.s_threshold
        elif command_args.s_nsigma:
            s_nsigma = command_args.s_nsigma
        else:
            assert False
    elif command_args.second_mask_shape is not None:
        sp_global_def.sxprint("Prepare second mask")
        nx = mask_first.get_xsize()
        ny = mask_first.get_ysize()
        nz = mask_first.get_zsize()
        if command_args.second_mask_shape == "cube":
            s_nx = command_args.s_nx
            s_ny = command_args.s_ny
            s_nz = command_args.s_nz
            s_mask = sp_utilities.model_blank(s_nx, s_ny, s_nz, 1)
        elif command_args.second_mask_shape == "cylinder":
            s_radius = command_args.s_radius
            s_nx = command_args.s_nx
            s_ny = command_args.s_ny
            s_nz = command_args.s_nz
            try:
                s_mask = sp_utilities.model_cylinder(s_radius, s_nx, s_ny, s_nz)
            except RuntimeError as e:
                sp_global_def.sxprint("An error occured! Please check the error log")
                raise
        elif command_args.second_mask_shape == "sphere":
            s_radius = command_args.s_radius
            s_nx = command_args.s_nx
            s_ny = command_args.s_ny
            s_nz = command_args.s_nz
            try:
                s_mask = sp_utilities.model_circle(s_radius, s_nx, s_ny, s_nz)
            except RuntimeError as e:
                sp_global_def.sxprint("An error occured! Please check the error log")
                raise
        else:
            assert False
        s_mask = sp_utilities.pad(s_mask, nx, ny, nz, 0)

    if s_mask is not None:
        sp_global_def.sxprint("Create second mask")

        if command_args.s_edge_type == "cosine":
            mode = "C"
        elif command_args.s_edge_type == "gaussian":
            mode = "G"
        else:
            assert False

        s_mask = sp_morphology.adaptive_mask_scipy(
            s_mask,
            nsigma=s_nsigma,
            threshold=s_density_threshold,
            ndilation=command_args.s_ndilation,
            nerosion=command_args.s_nerosion,
            edge_width=command_args.s_edge_width,
            allow_disconnected=command_args.s_allow_disconnected,
            mode=mode,
            do_approx=command_args.s_do_old,
            do_fill=command_args.s_fill_mask,
            do_print=True,
        )
        if command_args.s_invert:
            s_mask = 1 - s_mask
        sp_global_def.sxprint("Write outputs.")
        mask_first.write_image(output_prefix + "_mask_first.hdf")
        s_mask.write_image(output_prefix + "_mask_second.hdf")
        masked_combined = mask_first * s_mask
        masked_combined.write_image(output_prefix + "_mask.hdf")
    else:
        sp_global_def.sxprint("Write outputs.")
        mask_first.write_image(output_prefix + "_mask.hdf")

def main():
    sp_global_def.print_timestamp("Start")
    sp_global_def.BATCH = True
    run()
    sp_global_def.BATCH = False
    sp_global_def.sxprint("Done!")
    sp_global_def.print_timestamp("Finish")

if __name__ == "__main__":
    main()
