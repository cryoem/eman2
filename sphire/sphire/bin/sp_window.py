#!/usr/bin/env python

# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2016-2019 (toshio.moriya@kek.jp)
# Author: T. Durmaz 08/29/2014 (tunay.durmaz@uth.tmc.edu)
#
# Copyright (c) 2016-2019 Max Planck Institute of Molecular Physiology
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPHIRE software packages have some GPL dependencies,
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
from __future__ import division
from past.utils import old_div

import EMAN2_cppwrap
import EMAN2db
import EMAN2jsondb
import glob
import inspect
import mpi
import numpy
import optparse
import os
import shutil
from ..libpy import sp_applications
from ..libpy import sp_filter
from ..libpy import sp_fundamentals
from ..libpy import sp_global_def
from ..libpy import sp_morphology
from ..libpy import sp_statistics
from ..libpy import sp_utilities
import sys
import time
from builtins import range




# ========================================================================================
# Define functions for reading coordinates files of different formats.
# One of these will be used in main() through a first-class data type variable of Python
# (like function pointer in C/C++)
# This way, switch statement is unnecessary inside of the coordinates loop.
# ========================================================================================


def read_sphire_coords_file(coords_path):
    # 	coords_list = read_text_row(coords_path)
    coords_list = sp_utilities.read_text_row(coords_path, skip="#")
    return coords_list


def read_eman1_coords_file(coords_path):
    coords_list = sp_utilities.read_text_row(coords_path)
    for i in range(len(coords_list)):
        coords_list[i] = [
            (coords_list[i][0] + old_div(coords_list[i][2], 2)),
            (coords_list[i][1] + old_div(coords_list[i][3], 2)),
        ]
    return coords_list


def read_eman2_coords_file(coords_path):
    coords_list = EMAN2jsondb.js_open_dict(coords_path)["boxes"]
    for i in range(len(coords_list)):
        coords_list[i] = [coords_list[i][0], coords_list[i][1]]
    return coords_list


def read_spider_coords_file(coords_path):
    coords_list = sp_utilities.read_text_row(coords_path)
    for i in range(len(coords_list)):
        coords_list[i] = [coords_list[i][2], coords_list[i][3]]
    return coords_list


def read_cryolo_helical_segmented_coords_file(coords_path):
    """
	Reads segmented helical coordinates as writtin by crYOLO and sxhelixboxer.
	Furthermore it add the filament id, segment id and estimates the angle of the segment at
	each position.

	:param coords_path: Path to the coordinate file
	:return: list with entries [x,y,filamend_id,segment_id,angle_of_segment]
	"""
    try:
        split_indicis = []
        index_first_helix = -1

        with open(coords_path, "r") as csvfile:
            for index, row in enumerate(csvfile.readlines()):
                if row.startswith("#helix"):
                    if index_first_helix == -1:
                        index_first_helix = index
                    else:
                        split_indicis.append(
                            index - index_first_helix - (len(split_indicis) + 1)
                        )
        if index_first_helix == -1:
            return []
        filament_name = "{{0}}_{0:05d}"
        coords_list = []
        coordinates = numpy.atleast_2d(numpy.genfromtxt(coords_path))
        coord_filaments = numpy.split(coordinates, split_indicis)
        for filament_id, filament in enumerate(coord_filaments):
            curr_filaments_coords = []
            if len(filament) <= 1:
                print('Nr. of segments for this filament is {0}! Cannot calculate the angle necessary for subsequent filamentous processing steps. Therefore this filament will be skipped. If you do not want to skip those, please follow the SPA approach and use sphire instead of cryolo_helical_segmented as the --coordinates_format flag and do not provide the --filament_width flag in subsequent processes.'.format(len(filament)))
                continue

            for segment_id, coords in enumerate(filament):
                if segment_id == 0 and len(filament) > 1:
                    angle = estimate_angle(coords, filament[segment_id + 1])
                elif segment_id == (len(filament) - 1):
                    angle = estimate_angle(coords, filament[segment_id - 1])
                else:
                    angle = estimate_angle(
                        filament[segment_id - 1], filament[segment_id + 1]
                    )
                if segment_id == 0 and len(filament) > 1:
                    helical_track_length = 0
                else:
                    helical_track_length += estimate_distance(
                        coords, filament[segment_id - 1]
                    )
                curr_filaments_coords.append(
                    [
                        float(coords[0]),
                        float(coords[1]),
                        filament_name.format(filament_id),
                        segment_id,
                        angle,
                        helical_track_length,
                    ]
                )

            coords_list.extend(curr_filaments_coords)
        return coords_list
    except:
        return []


# ========================================================================================
#  Helper functions
# ========================================================================================


def estimate_distance(coords_a, coords_b):
    """
	Estimates the distance of a line between two  boxes.
	:param coords_a: First coordinate pair
	:param coords_b: Second coordinate pair
	:return: Distance of the line connectiong coords_a and coords_b.
	"""
    delta_x = coords_b[0] - coords_a[0]
    delta_y = coords_b[1] - coords_a[1]
    distance = numpy.sqrt(delta_x ** 2 + delta_y ** 2)

    return distance


def estimate_angle(coords_a, coords_b):
    """
	Estimates the angle given by a line between two  boxes.
	:param coords_a: First coordinate pair
	:param coords_b: Second coordinate pair
	:return: Angle from positive x axis twoards the positiv y axis.
	"""
    delta_x = coords_b[0] - coords_a[0]
    delta_y = coords_b[1] - coords_a[1]
    angle = old_div(-1 * numpy.arctan2(delta_y, delta_x) * 180, numpy.pi)

    return angle


def get_cmd_line():
    """
	Generate command line
	:return: Returns command line
	"""
    cmd_line = ""
    for arg in sys.argv:
        cmd_line += arg + "  "
    cmd_line = "Shell line command: " + cmd_line
    return cmd_line


def get_time_stamp_suffix():
    """
	:return: Returns suffix of current time stamp
	"""

    time_stamp_suffix = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    return time_stamp_suffix


# ----------------------------------------------------------------------------------------
#  Data type checker
# ----------------------------------------------------------------------------------------


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


# ========================================================================================
#  Main function
# ========================================================================================


def run():
    program_name = os.path.basename(sys.argv[0])
    usage = (
        program_name
        + """  input_micrograph_pattern  input_coordinates_pattern  input_ctf_params_source  output_directory  --selection_list=selection_list  --coordinates_format  --box_size=box_size  --skip_invert  --limit_ctf  --astigmatism_error=astigmatism_error  --resample_ratio=resample_ratio  --check_consistency

Window particles from micrographs using the particles coordinates.

All Micrographs Mode - Process all micrographs in a directory:
	Specify path pattern of input micrographs and coordinates files with a wild card (*). 
	Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). 
	The path pattern must be enclosed by single quotes (') or double quotes ("). (Note: sxgui.py automatically adds single quotes (')). 
	The substring at the variable part must be same between a associated pair of input micrograph and coordinates file.
	BDB files can not be selected as input micrographs.
	Next, specify the source of CTF paramters. 
	For cryo data, this should be the file produced by sxcter and normally called partres.txt. 
	For negative staining data, it should be the pixel size [A/Pixels] of input micrographs.
	Finally, specify output directory where all outputs should be saved.
	In this mode, all micrographs matching the path pattern will be processed.

	mpirun  -np  32  sp_window.py  './mic*.hdf'  'info/mic*_info.json'  outdir_cter/partres/partres.txt  particles  --coordinates_format=eman2  --box_size=64

Selected Micrographs Mode - Process all micrographs in a selection list file:
	In addition to input micrographs path pattern, coordinates files path pattern, CTF paramters source, and output directry arguments, 
	specify a name of micrograph selection list text file using --selection_list option.
	In this mode, only micrographs in the selection list which matches the file name part of the pattern (ignoring the directory paths) will be processed.
	If a micrograph name in the selection list does not exists in the directory specified by the micrograph path pattern, processing of the micrograph will be skipped.

	mpirun  -np  32  sp_window.py  './mic*.hdf'  'info/mic*_info.json'  outdir_cter/partres/partres.txt  particles  --selection_list=mic_list.txt  --coordinates_format=eman2  --box_size=64

Single Micrograph Mode - Process a single micrograph:
	In addition to input micrographs path pattern, coordinates files path pattern, CTF paramters source, and output directry arguments, 
	specify a single micrograph name using --selection_list option.
	In this mode, only the specified single micrograph will be processed.
	If this micrograph name does not matches the file name part of the pattern (ignoring the directory paths), the process will exit without processing it.
	If this micrograph name matches the file name part of the pattern but does not exists in the directory which specified by the micrograph path pattern, again the process will exit without processing it.
	Use single processor for this mode. 

	sp_window.py  './mic*.hdf'  'info/mic*_info.json'  outdir_cter/partres/partres.txt  particles  --selection_list=mic0.hdf  --coordinates_format=eman2  --box_size=64

For negative staining data, set the pixel size [A/Pixels] as the source of CTF paramters and use --skip_invert.

	mpirun  -np  32  sp_window.py  './mic*.hdf'  'info/mic*_info.json'  5.2  particles  --coordinates_format=eman2  --box_size=64  --skip_invert

"""
    )
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)
    parser.add_option(
        "--selection_list",
        type="string",
        default=None,
        help="Micrograph selecting list: Specify a name of micrograph selection list text file for Selected Micrographs Mode. The file extension must be '.txt'. Alternatively, the file name of a single micrograph can be specified for Single Micrograph Mode. (default none)",
    )
    parser.add_option(
        "--coordinates_format",
        type="string",
        default="cryolo",
        help="Coordinate file format: Allowed values are 'sphire', 'eman1', 'eman2', or 'spider', 'cryolo', 'cryolo_helical_segmented'. The sphire, eman2, and spider formats use the particle center as coordinates. The eman1, cryolo, and cryolo_helical_segmented' format uses the lower left corner of the box as coordinates. cryolo_helical_segmented is for already segmented helical coordinates. (default cryolo)",
    )
    parser.add_option(
        "--box_size",
        type="int",
        default=256,
        help="Particle box size [Pixels]: The x and y dimensions of square area to be windowed. The box size after resampling is assumed when resample_ratio < 1.0. (default 256)",
    )
    parser.add_option(
        "--skip_invert",
        action="store_true",
        default=False,
        help="Skip invert image contrast: Use this option for negative staining data. By default, the image contrast is inverted for cryo data. (default False)",
    )
    parser.add_option(
        "--limit_ctf",
        action="store_true",
        default=False,
        help="Use CTF limit filter: Frequencies where CTF oscillations can not be properly modeled with the resampled pixel size will be discarded in the images with the appropriate low-pass filter. This has no effects when the CTER partres file is not specified by the CTF paramters source argument. (default False)",
    )
    parser.add_option(
        "--astigmatism_error",
        type="float",
        default=360.0,
        help="Astigmatism error limit [Degrees]: Set astigmatism to zero for all micrographs where the angular error computed by sxcter is larger than the desired value. This has no effects when the CTER partres file is not specified by the CTF paramters source argument. (default 360.0)",
    )
    parser.add_option(
        "--resample_ratio",
        type="float",
        default=1.0,
        help="Ratio between new and original pixel size: Use a value between 0.0 and 1.0 (excluding 0.0). The new pixel size will be automatically recalculated and stored in CTF paramers when resample_ratio < 1.0 is used. (default 1.0)",
    )
    parser.add_option(
        "--check_consistency",
        action="store_true",
        default=False,
        help="Check consistency of dataset: Create a text file containing the list of Micrograph ID entries might have inconsistency among the provided dataset. (i.e. mic_consistency_check_info_TIMESTAMP.txt). (default False)",
    )
    parser.add_option(
        "--filament_width",
        type="int",
        default=-1,
        help="Filament width [Pixels]: Filament width for the creation of the rectangular mask. Default is one third of the box size. (default -1)",
    )

    (options, args) = parser.parse_args(sys.argv[1:])

    # ====================================================================================
    # Prepare processing
    # ====================================================================================
    # ------------------------------------------------------------------------------------
    # Set up MPI related variables
    # ------------------------------------------------------------------------------------
    # Detect if program is running under MPI
    RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ

    main_mpi_proc = 0
    if RUNNING_UNDER_MPI:
        my_mpi_proc_id = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        n_mpi_procs = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
    else:
        my_mpi_proc_id = 0
        n_mpi_procs = 1

        # ------------------------------------------------------------------------------------
        # Set up SPHIRE global definitions
        # ------------------------------------------------------------------------------------
    if sp_global_def.CACHE_DISABLE:
        sp_utilities.disable_bdb_cache()

        # Change the name log file for error message
    original_logfilename = sp_global_def.LOGFILE
    sp_global_def.LOGFILE = (
        os.path.splitext(program_name)[0] + "_" + original_logfilename + ".txt"
    )

    # ------------------------------------------------------------------------------------
    # Print command line
    # ------------------------------------------------------------------------------------
    if my_mpi_proc_id == main_mpi_proc:
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("%s" % get_cmd_line())
        # sxprint(" ")

        # ------------------------------------------------------------------------------------
        # Check error conditions of arguments and options, then prepare variables for arguments
        # ------------------------------------------------------------------------------------

    mic_pattern = None
    coords_pattern = None
    ctf_params_src = None
    root_out_dir = None
    # Not a real while, each "if" statement has the opportunity to use break when errors need to be reported
    error_status = None
    while True:
        # --------------------------------------------------------------------------------
        # Check the number of arguments. If OK, then prepare variables for them
        # --------------------------------------------------------------------------------
        if error_status is None and len(args) != 4:
            error_status = (
                "Please check usage for number of arguments.\n Usage: "
                + usage
                + "\n"
                + "Please run %s -h for help." % (program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        assert len(args) == 4
        mic_pattern = args[0]
        coords_pattern = args[1]
        ctf_params_src = args[2]
        root_out_dir = args[3]

        # --------------------------------------------------------------------------------
        # Check error conditions of arguments
        # --------------------------------------------------------------------------------
        if error_status is None and mic_pattern[: len("bdb:")].lower() == "bdb":
            error_status = (
                "BDB file can not be selected as input micrographs. Please convert the format, and restart the program. Run %s -h for help."
                % (program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        if error_status is None and mic_pattern.find("*") == -1:
            error_status = (
                "Input micrograph file name pattern must contain wild card (*). Please check input_micrograph_pattern argument. Run %s -h for help."
                % (program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        if error_status is None and coords_pattern.find("*") == -1:
            error_status = (
                "Input coordinates file name pattern must contain wild card (*). Please check input_coordinates_pattern argument. Run %s -h for help."
                % (program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        if not is_float(ctf_params_src):
            assert type(ctf_params_src) is str
            # This should be string for CTER partres (CTF parameter) file
            if error_status is None and os.path.exists(ctf_params_src) == False:
                error_status = (
                    "Specified CTER partres file is not found. Please check input_ctf_params_source argument. Run %s -h for help."
                    % (program_name),
                    inspect.getframeinfo(inspect.currentframe()),
                )
                break
        else:
            assert is_float(ctf_params_src)
            if error_status is None and float(ctf_params_src) <= 0.0:
                error_status = (
                    "Specified pixel size is not larger than 0.0. Please check input_ctf_params_source argument. Run %s -h for help."
                    % (program_name),
                    inspect.getframeinfo(inspect.currentframe()),
                )
                break

        if error_status is None and os.path.exists(root_out_dir):
            error_status = (
                "Output directory exists. Please change the name and restart the program.",
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

            # --------------------------------------------------------------------------------
            # Check error conditions of options
            # --------------------------------------------------------------------------------
        if options.selection_list != None:
            if error_status is None and not os.path.exists(options.selection_list):
                error_status = (
                    "File specified by selection_list option does not exists. Please check selection_list option. Run %s -h for help."
                    % (program_name),
                    inspect.getframeinfo(inspect.currentframe()),
                )
                break

        if error_status is None and options.coordinates_format.lower() not in [
            "sphire",
            "eman1",
            "eman2",
            "spider",
            "cryolo",
            "cryolo_helical_segmented",
        ]:
            error_status = (
                "Invalid option value: --coordinates_format=%s. Please run %s -h for help."
                % (options.coordinates_format, program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        if error_status is None and (options.box_size <= 0):
            error_status = (
                "Invalid option value: --box_size=%s. The box size must be an interger larger than zero. Please run %s -h for help."
                % (options.box_size, program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        if error_status is None and (
            options.resample_ratio <= 0.0 or options.resample_ratio > 1.0
        ):
            error_status = (
                "Invalid option value: --resample_ratio=%s. Please run %s -h for help."
                % (options.resample_ratio, program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        break
    sp_utilities.if_error_then_all_processes_exit_program(error_status)
    assert mic_pattern != None
    assert coords_pattern != None
    assert ctf_params_src != None
    assert root_out_dir != None

    # ------------------------------------------------------------------------------------
    # Check warning conditions of options
    # ------------------------------------------------------------------------------------
    if my_mpi_proc_id == main_mpi_proc:
        # This should be string specifying pixel size [A/Pixels]
        assert type(ctf_params_src) is str
        if is_float(ctf_params_src):
            if options.limit_ctf:
                sp_global_def.sxprint(
                    "WARNING!!! --limit_ctf option has no effects since the CTER partres file is not specified with input_ctf_params_source argument..."
                )

            if (
                options.astigmatism_error != 360.0
            ):  # WARN User only when it is obvious that user is trying to use astigmatism_error option
                sp_global_def.sxprint(
                    "WARNING!!! --astigmatism_error option has no effects since the CTER partres file is not specified with input_ctf_params_source argument..."
                )

                # ------------------------------------------------------------------------------------
                # Check the source of CTF parameteres and select the CTF mode
                # (1) Real CTF parameters mode  : Use real CTF parameters stored in a CTER partres (CTF parameter) file for cryo data (use_real_ctf_params is True)
                # (2) Dummy CTF parameters mode : Create dummy CTF parameters for negative staining data (use_real_ctf_params is False)
                # ------------------------------------------------------------------------------------
    use_real_ctf_params = not is_float(ctf_params_src)

    # NOTE: 2017/12/05 Toshio Moriya
    # The following code is to support the old format of CTER partres file. It should be removed near future
    #
    # Define enumerators for indices of parameters of each CTER partres entry in the old format(BEFIRE 2017/12/05).
    # All mpi processes must have access to these indices
    i_enum = -1
    i_enum += 1
    idx_old_cter_def = i_enum  # defocus [um]; index must be same as ctf object format
    i_enum += 1
    idx_old_cter_cs = i_enum  # Cs [mm]; index must be same as ctf object format
    i_enum += 1
    idx_old_cter_vol = i_enum  # voltage[kV]; index must be same as ctf object format
    i_enum += 1
    idx_old_cter_apix = (
        i_enum
    )  # pixel size [A]; index must be same as ctf object format
    i_enum += 1
    idx_old_cter_bfactor = (
        i_enum
    )  # B-factor [A^2]; index must be same as ctf object format
    i_enum += 1
    idx_old_cter_total_ac = (
        i_enum
    )  # total amplitude contrast [%]; index must be same as ctf object format
    i_enum += 1
    idx_old_cter_astig_amp = (
        i_enum
    )  # astigmatism amplitude [um]; index must be same as ctf object format
    i_enum += 1
    idx_old_cter_astig_ang = (
        i_enum
    )  # astigmatism angle [degree]; index must be same as ctf object format
    i_enum += 1
    idx_old_cter_sd_def = i_enum  # std dev of defocus [um]
    i_enum += 1
    idx_old_cter_sd_total_ac = i_enum  # std dev of total amplitude contrast [%]
    i_enum += 1
    idx_old_cter_sd_astig_amp = i_enum  # std dev of ast amp [A]
    i_enum += 1
    idx_old_cter_sd_astig_ang = i_enum  # std dev of ast angle [degree]
    i_enum += 1
    idx_old_cter_cv_def = i_enum  # coefficient of variation of defocus [%]
    i_enum += 1
    idx_old_cter_cv_astig_amp = i_enum  # coefficient of variation of ast amp [%]
    i_enum += 1
    idx_old_cter_spectra_diff = (
        i_enum
    )  # average of differences between with- and without-astig. experimental 1D spectra at extrema
    i_enum += 1
    idx_old_cter_error_def = (
        i_enum
    )  # frequency at which signal drops by 50% due to estimated error of defocus alone [1/A]
    i_enum += 1
    idx_old_cter_error_astig = (
        i_enum
    )  # frequency at which signal drops by 50% due to estimated error of defocus and astigmatism [1/A]
    i_enum += 1
    idx_old_cter_error_ctf = i_enum  # limit frequency by CTF error [1/A]
    i_enum += 1
    idx_old_cter_mic_name = i_enum  # micrograph name
    i_enum += 1
    n_idx_old_cter = i_enum

    # Define enumerators for indices of parameters of each CTER partres entry in the new format (AFTER 2017/12/05).
    # All mpi processes must have access to these indices
    i_enum = -1
    i_enum += 1
    idx_cter_def = i_enum  # defocus [um]; index must be same as ctf object format
    i_enum += 1
    idx_cter_cs = i_enum  # Cs [mm]; index must be same as ctf object format
    i_enum += 1
    idx_cter_vol = i_enum  # voltage[kV]; index must be same as ctf object format
    i_enum += 1
    idx_cter_apix = i_enum  # pixel size [A]; index must be same as ctf object format
    i_enum += 1
    idx_cter_bfactor = i_enum  # B-factor [A^2]; index must be same as ctf object format
    i_enum += 1
    idx_cter_total_ac = (
        i_enum
    )  # total amplitude contrast [%]; index must be same as ctf object format
    i_enum += 1
    idx_cter_astig_amp = (
        i_enum
    )  # astigmatism amplitude [um]; index must be same as ctf object format
    i_enum += 1
    idx_cter_astig_ang = (
        i_enum
    )  # astigmatism angle [degree]; index must be same as ctf object format
    i_enum += 1
    idx_cter_sd_def = i_enum  # std dev of defocus [um]
    i_enum += 1
    idx_cter_sd_total_ac = i_enum  # std dev of total amplitude contrast [%]
    i_enum += 1
    idx_cter_sd_astig_amp = i_enum  # std dev of astigmatism amp [A]
    i_enum += 1
    idx_cter_sd_astig_ang = i_enum  # std dev of astigmatism angle [degree]
    i_enum += 1
    idx_cter_cv_def = i_enum  # coefficient of variation of defocus [%]
    i_enum += 1
    idx_cter_cv_astig_amp = i_enum  # coefficient of variation of astigmatism amp [%]
    i_enum += 1
    idx_cter_error_def = (
        i_enum
    )  # frequency at which signal drops by 50% due to estimated error of defocus alone [1/A]
    i_enum += 1
    idx_cter_error_astig = (
        i_enum
    )  # frequency at which signal drops by 50% due to estimated error of defocus and astigmatism [1/A]
    i_enum += 1
    idx_cter_error_ctf = i_enum  # limit frequency by CTF error [1/A]
    i_enum += 1
    idx_cter_max_freq = (
        i_enum
    )  # visual-impression-based maximum frequency limit [A] (e.g. max frequency of relion; CCC between neighbour zero-crossing pair)
    i_enum += 1
    idx_cter_reserved = (
        i_enum
    )  # reserved spot for maximum frequency limit or error criterion. possibly originated from external program (e.g. CTF figure of merit of RELION)
    i_enum += 1
    idx_cter_const_ac = i_enum  # constant amplitude contrast [%]
    i_enum += 1
    idx_cter_phase_shift = i_enum  # phase shift [degrees]
    i_enum += 1
    idx_cter_mic_name = i_enum  # micrograph name
    i_enum += 1
    n_idx_cter = i_enum

    # ------------------------------------------------------------------------------------
    # Prepare the variables for all sections
    # ------------------------------------------------------------------------------------
    # Micrograph basename pattern (directory path is removed from micrograph path pattern)
    mic_basename_pattern = os.path.basename(mic_pattern)

    # Global entry dictionary (all possible entries from all lists) for all mic id substring
    global_entry_dict = {}  # mic id substring is the key
    subkey_input_mic_path = "Input Micrograph Path"
    subkey_selected_mic_basename = "Selected Micrograph Basename"
    subkey_coords_path = "Input Coordinates File Path"
    subkey_cter_entry = "CTER Partres Entry"

    # List keeps only id substrings of micrographs whose all necessary information are available
    valid_mic_id_substr_list = []

    # ====================================================================================
    # Obtain the list of micrograph id sustrings using a single CPU (i.e. main mpi process)
    # ====================================================================================
    # NOTE: Toshio Moriya 2016/10/24
    # The below is not a real while.
    # It gives if-statements an opportunity to use break when errors need to be reported
    # However, more elegant way is to use 'raise' statement of exception mechanism...
    #
    error_status = None
    while my_mpi_proc_id == main_mpi_proc:
        # --------------------------------------------------------------------------------
        # Prepare variables for this section
        # --------------------------------------------------------------------------------
        # Prefix and suffix of micrograph basename pattern
        # to find the head/tail indices of micrograph id substring
        mic_basename_tokens = mic_basename_pattern.split("*")
        assert len(mic_basename_tokens) == 2
        # Find head index of micrograph id substring
        mic_id_substr_head_idx = len(mic_basename_tokens[0])

        # Prefix and suffix of coordinates file path pattern
        # to find the head/tail indices of coordinates file id substring
        coords_pattern_tokens = coords_pattern.split("*")
        assert len(coords_pattern_tokens) == 2
        # Find head index of coordinates id substring
        coords_id_substr_head_idx = len(coords_pattern_tokens[0])

        # --------------------------------------------------------------------------------
        # Register micrograph id substrings found in the input directory (specified by micrograph path pattern)
        # to the global entry dictionary
        # --------------------------------------------------------------------------------
        # Generate the list of micrograph paths in the input directory
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Checking the input directory...")
        input_mic_path_list = glob.glob(mic_pattern)
        # Check error condition of input micrograph file path list
        sp_global_def.sxprint(
            "Found %d microgarphs in %s."
            % (len(input_mic_path_list), os.path.dirname(mic_pattern))
        )
        if error_status is None and not input_mic_path_list:
            error_status = (
                "No micrograph files are found in the directory specified by micrograph path pattern (%s). Please check input_micrograph_pattern argument. Run %s -h for help."
                % (os.path.dirname(mic_pattern), program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break
        assert input_mic_path_list

        # Register micrograph id substrings to the global entry dictionary
        for input_mic_path in input_mic_path_list:
            # Find tail index of micrograph id substring and extract the substring from the micrograph name
            input_mic_basename = os.path.basename(input_mic_path)
            mic_id_substr_tail_idx = input_mic_basename.rindex(mic_basename_tokens[1])
            mic_id_substr = input_mic_basename[
                mic_id_substr_head_idx:mic_id_substr_tail_idx
            ]
            assert input_mic_path == mic_pattern.replace("*", mic_id_substr)
            if not mic_id_substr in global_entry_dict:
                # sxprint("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
                global_entry_dict[mic_id_substr] = {}
            assert mic_id_substr in global_entry_dict
            global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path
        assert global_entry_dict

        # --------------------------------------------------------------------------------
        # Register micrograph id substrings found in the selection list
        # to the global entry dictionary
        # --------------------------------------------------------------------------------
        # Generate the list of selected micrograph paths in the selection file
        selected_mic_path_list = []
        # Generate micrograph lists according to the execution mode
        if options.selection_list == None:
            sp_global_def.sxprint(" ")
            sp_global_def.sxprint("----- Running with All Micrographs Mode -----")
            # Treat all micrographs in the input directory as selected ones
            selected_mic_path_list = input_mic_path_list
        else:
            assert options.selection_list != None
            if os.path.splitext(options.selection_list)[1] == ".txt":
                sp_global_def.sxprint(" ")
                sp_global_def.sxprint(
                    "----- Running with Selected Micrographs Mode -----"
                )
                sp_global_def.sxprint(" ")
                sp_global_def.sxprint("Checking the selection list...")
                assert os.path.exists(options.selection_list)
                selected_mic_path_list = sp_utilities.read_text_file(
                    options.selection_list
                )

                # Check error condition of micrograph entry lists
                sp_global_def.sxprint(
                    "Found %d microgarph entries in %s."
                    % (len(selected_mic_path_list), options.selection_list)
                )
                if error_status is None and len(selected_mic_path_list) == 0:
                    error_status = (
                        "No micrograph entries are found in the selection list file. Please check selection_list option. Run %s -h for help."
                        % (program_name),
                        inspect.getframeinfo(inspect.currentframe()),
                    )
                    break
                assert len(selected_mic_path_list) > 1
                if error_status is None and not isinstance(
                    selected_mic_path_list[0], str
                ):
                    error_status = (
                        "Invalid format of the selection list file. The first column must contain micrograph paths in string type. Please check selection_list option. Run %s -h for help."
                        % (program_name),
                        inspect.getframeinfo(inspect.currentframe()),
                    )
                    break
            else:
                sp_global_def.sxprint(" ")
                sp_global_def.sxprint("----- Running with Single Micrograph Mode -----")
                sp_global_def.sxprint(" ")
                sp_global_def.sxprint(
                    "Processing a single micrograph: %s..." % (options.selection_list)
                )
                selected_mic_path_list = [options.selection_list]
            assert selected_mic_path_list

            selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
            if selected_mic_directory != "":
                sp_global_def.sxprint(
                    "    NOTE: Program disregards the directory paths in the selection list (%s)."
                    % (selected_mic_directory)
                )

        assert selected_mic_path_list

        # Register micrograph id substrings to the global entry dictionary
        for selected_mic_path in selected_mic_path_list:
            # Find tail index of micrograph id substring and extract the substring from the micrograph name
            selected_mic_basename = os.path.basename(selected_mic_path)
            mic_id_substr_tail_idx = selected_mic_basename.rindex(mic_basename_tokens[1])
            mic_id_substr = selected_mic_basename[
                mic_id_substr_head_idx:mic_id_substr_tail_idx
            ]
            if (
                error_status is None
                and selected_mic_basename
                != mic_basename_pattern.replace("*", mic_id_substr)
            ):
                error_status = (
                    "A micrograph name (%s) in the input directory (%s) does not match with input micrograph basename pattern (%s) (The wild card replacement with '%s' resulted in '%s'). Please correct input micrograph path pattern. Run %s -h for help."
                    % (
                        selected_mic_basename,
                        os.path.dirname(mic_pattern),
                        mic_basename_pattern,
                        mic_id_substr,
                        mic_basename_pattern.replace("*", mic_id_substr),
                        program_name,
                    ),
                    inspect.getframeinfo(inspect.currentframe()),
                )
                break
            if not mic_id_substr in global_entry_dict:
                # sxprint("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list " % (mic_id_substr))
                global_entry_dict[mic_id_substr] = {}
            assert mic_id_substr in global_entry_dict
            global_entry_dict[mic_id_substr][
                subkey_selected_mic_basename
            ] = selected_mic_basename
        assert global_entry_dict

        del selected_mic_path_list  # Do not need this anymore
        del input_mic_path_list  # Do not need this anymore

        # --------------------------------------------------------------------------------
        # Register coordinates id substrings in coordinate path list to the global entry dictionary.
        # coordinates id substring (coords_id_substr) and micrograph id substring (mic_id_substr)
        # should be the same for the associated pair of micrograph and coordnates file.
        # --------------------------------------------------------------------------------
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Checking the coordinates files...")
        coords_path_list = glob.glob(coords_pattern)

        # Check error condition of coordinates file path list
        sp_global_def.sxprint(
            "Found %d coordinates files in %s directory."
            % (len(coords_path_list), os.path.dirname(coords_pattern))
        )

        if error_status is None and not coords_path_list:
            error_status = (
                "No coordinates files are found in the directory specified by coordinates file path pattern (%s). Please check input_coordinates_pattern argument. Run %s -h for help."
                % (os.path.dirname(coords_pattern), program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break
        assert coords_path_list

        for coords_path in coords_path_list:

            # Find tail index of coordinates id substring and extract the substring from the coordinates file path
            coords_id_substr_tail_idx = coords_path.rindex(coords_pattern_tokens[1])
            coords_id_substr = coords_path[
                coords_id_substr_head_idx:coords_id_substr_tail_idx
            ]

            assert coords_path == coords_pattern.replace("*", coords_id_substr)

            if not coords_id_substr in global_entry_dict:
                global_entry_dict[coords_id_substr] = {}

            assert coords_id_substr in global_entry_dict

            global_entry_dict[coords_id_substr][subkey_coords_path] = coords_path

        assert global_entry_dict

        del coords_path_list  # Do not need this anymore

        # --------------------------------------------------------------------------------
        # If necessary, register micrograph id substrings of CTER partres entries to the global entry dictionary
        # --------------------------------------------------------------------------------

        if use_real_ctf_params:
            # This should be string for CTER partres (CTF parameter) file
            sp_global_def.sxprint(" ")
            sp_global_def.sxprint("Checking the CTER partres file...")
            assert os.path.exists(ctf_params_src)
            cter_entry_list = sp_utilities.read_text_row(ctf_params_src)

            # Check error condition of CTER partres entry list
            sp_global_def.sxprint(
                "Found %d CTER partres entries in %s."
                % (len(cter_entry_list), ctf_params_src)
            )

            if error_status is None and not cter_entry_list:
                error_status = (
                    "No CTER partres entries are found in %s. Please check input_ctf_params_source argument. Run %s -h for help."
                    % (ctf_params_src, program_name),
                    inspect.getframeinfo(inspect.currentframe()),
                )
                break

            assert cter_entry_list

            # NOTE: 2017/12/05 Toshio Moriya
            # The following code is to support the old format of CTER partres file. It should be removed near future
            if (
                error_status is None
                and len(cter_entry_list[0]) != n_idx_cter
                and len(cter_entry_list[0]) != n_idx_old_cter
            ):
                error_status = (
                    "The number of columns (%d) has to be %d in %s."
                    % (len(cter_entry_list[0]), n_idx_cter, ctf_params_src),
                    inspect.getframeinfo(inspect.currentframe()),
                )
                break

            assert (
                len(cter_entry_list[0]) == n_idx_cter
                or len(cter_entry_list[0]) == n_idx_old_cter
            )

            # For NEW CTER partres format (AFTER 2017/12/05)
            if len(cter_entry_list[0]) == n_idx_cter:
                cter_mic_directory = os.path.dirname(
                    cter_entry_list[0][idx_cter_mic_name]
                )
                if cter_mic_directory != "":
                    sp_global_def.sxprint(
                        "    NOTE: Program disregards the directory paths in the CTER partres file (%s)."
                        % (cter_mic_directory)
                    )

                for cter_entry in cter_entry_list:
                    # Find tail index of micrograph id substring and extract the substring from the micrograph path of CTER partres entry
                    cter_mic_path = cter_entry[idx_cter_mic_name]
                    cter_mic_basename = os.path.basename(cter_mic_path)
                    mic_id_substr_tail_idx = cter_mic_basename.rindex(
                        mic_basename_tokens[1]
                    )
                    mic_id_substr = cter_mic_basename[
                        mic_id_substr_head_idx:mic_id_substr_tail_idx
                    ]
                    # Between cter_mic_path and mic_path, directory paths might be different but the basenames should be same!
                    if (
                        error_status is None
                        and cter_mic_basename
                        != mic_basename_pattern.replace("*", mic_id_substr)
                    ):
                        error_status = (
                            "A micrograph name (%s) in the CTER partres file (%s) does not match with input micrograph basename pattern (%s) (The wild card replacement with '%s' resulted in '%s'). Please check the CTER partres file and correct input micrograph path pattern. Run %s -h for help."
                            % (
                                cter_mic_basename,
                                ctf_params_src,
                                mic_basename_pattern,
                                mic_id_substr,
                                mic_basename_pattern.replace("*", mic_id_substr),
                                program_name,
                            ),
                            inspect.getframeinfo(inspect.currentframe()),
                        )
                        break

                    if cter_entry[idx_cter_sd_astig_ang] > options.astigmatism_error:
                        sp_global_def.sxprint(
                            "    NOTE: Astigmatism angular SD of %s (%f degree) exceeds specified limit (%f degree). Resetting astigmatism parameters to zeros..."
                            % (
                                cter_mic_basename,
                                cter_entry[idx_cter_sd_astig_ang],
                                options.astigmatism_error,
                            )
                        )
                        cter_entry[idx_cter_astig_amp] = 0.0
                        cter_entry[idx_cter_astig_ang] = 0.0

                    if not mic_id_substr in global_entry_dict:
                        # sxprint("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from cter_entry_list " % (mic_id_substr))
                        global_entry_dict[mic_id_substr] = {}
                    assert mic_id_substr in global_entry_dict
                    global_entry_dict[mic_id_substr][subkey_cter_entry] = cter_entry
                assert global_entry_dict
                # For OLD CTER partres format (BEFORE 2017/12/05)
            else:
                assert len(cter_entry_list[0]) == n_idx_old_cter, "MKR_DEBUG"
                sp_global_def.sxprint(
                    "WARNING!!!: Number of columns is %d in the specified CTER partres file %s. The format might be old because the new one should contain %d colums."
                    % (len(cter_entry_list[0]), ctf_params_src, n_idx_cter)
                )
                sp_global_def.sxprint(
                    "            We will stop supporting the old format in near future. Please consider rerun CTER."
                )

                cter_mic_directory = os.path.dirname(
                    cter_entry_list[0][idx_old_cter_mic_name]
                )
                if cter_mic_directory != "":
                    sp_global_def.sxprint(
                        "    NOTE: Program disregards the directory paths in the CTER partres file (%s)."
                        % (cter_mic_directory)
                    )

                for old_cter_entry in cter_entry_list:
                    # Convert each entry to new format
                    #
                    # assume amplitude amplitude contrast is total amplitude constrast in [%], estimated as variable Volta phase shift. Conver it to [deg].
                    # Also, assuming constant amplitude contrast is zero since there is no information available in the old format.
                    #
                    cter_entry = [None] * n_idx_cter
                    cter_entry[idx_cter_def] = old_cter_entry[idx_old_cter_def]
                    cter_entry[idx_cter_cs] = old_cter_entry[idx_old_cter_cs]
                    cter_entry[idx_cter_vol] = old_cter_entry[idx_old_cter_vol]
                    cter_entry[idx_cter_apix] = old_cter_entry[idx_old_cter_apix]
                    cter_entry[idx_cter_bfactor] = old_cter_entry[idx_old_cter_bfactor]
                    cter_entry[idx_cter_total_ac] = old_cter_entry[
                        idx_old_cter_total_ac
                    ]
                    cter_entry[idx_cter_astig_amp] = old_cter_entry[
                        idx_old_cter_astig_amp
                    ]
                    cter_entry[idx_cter_astig_ang] = old_cter_entry[
                        idx_old_cter_astig_ang
                    ]
                    cter_entry[idx_cter_sd_def] = old_cter_entry[idx_old_cter_sd_def]
                    cter_entry[idx_cter_sd_total_ac] = old_cter_entry[
                        idx_old_cter_sd_total_ac
                    ]
                    cter_entry[idx_cter_sd_astig_amp] = old_cter_entry[
                        idx_old_cter_sd_astig_amp
                    ]
                    cter_entry[idx_cter_sd_astig_ang] = old_cter_entry[
                        idx_old_cter_sd_astig_ang
                    ]
                    cter_entry[idx_cter_cv_def] = old_cter_entry[idx_old_cter_cv_def]
                    cter_entry[idx_cter_cv_astig_amp] = old_cter_entry[
                        idx_old_cter_cv_astig_amp
                    ]
                    cter_entry[idx_cter_error_def] = old_cter_entry[
                        idx_old_cter_error_def
                    ]
                    cter_entry[idx_cter_error_astig] = old_cter_entry[
                        idx_old_cter_error_astig
                    ]
                    cter_entry[idx_cter_error_ctf] = old_cter_entry[
                        idx_old_cter_error_ctf
                    ]
                    cter_entry[idx_cter_max_freq] = old_div(
                        0.5, old_cter_entry[idx_old_cter_apix]
                    )  # Set to Nyquist frequency
                    cter_entry[idx_cter_reserved] = 0.0
                    cter_entry[idx_cter_const_ac] = 0.0
                    cter_entry[idx_cter_phase_shift] = sp_morphology.ampcont2angle(
                        old_cter_entry[idx_old_cter_total_ac]
                    )
                    cter_entry[idx_cter_mic_name] = old_cter_entry[
                        idx_old_cter_mic_name
                    ]
                    assert len(cter_entry) == n_idx_cter

                    # Find tail index of micrograph id substring and extract the substring from the micrograph path of CTER partres entry
                    cter_mic_path = cter_entry[idx_cter_mic_name]
                    cter_mic_basename = os.path.basename(cter_mic_path)
                    mic_id_substr_tail_idx = cter_mic_basename.rindex(
                        mic_basename_tokens[1]
                    )
                    mic_id_substr = cter_mic_basename[
                        mic_id_substr_head_idx:mic_id_substr_tail_idx
                    ]
                    # Between cter_mic_path and mic_path, directory paths might be different but the basenames should be same!
                    if (
                        error_status is None
                        and cter_mic_basename
                        != mic_basename_pattern.replace("*", mic_id_substr)
                    ):
                        error_status = (
                            "A micrograph name (%s) in the CTER partres file (%s) does not match with input micrograph basename pattern (%s) (The wild card replacement with '%s' resulted in '%s'). Please check the CTER partres file and correct input micrograph path pattern. Run %s -h for help."
                            % (
                                cter_mic_basename,
                                ctf_params_src,
                                mic_basename_pattern,
                                mic_id_substr,
                                mic_basename_pattern.replace("*", mic_id_substr),
                                program_name,
                            ),
                            inspect.getframeinfo(inspect.currentframe()),
                        )
                        break

                    if cter_entry[idx_cter_sd_astig_ang] > options.astigmatism_error:
                        sp_global_def.sxprint(
                            "    NOTE: Astigmatism angular SD of %s (%f degree) exceeds specified limit (%f degree). Resetting astigmatism parameters to zeros..."
                            % (
                                cter_mic_basename,
                                cter_entry[idx_cter_sd_astig_ang],
                                options.astigmatism_error,
                            )
                        )

                        cter_entry[idx_cter_astig_amp] = 0.0
                        cter_entry[idx_cter_astig_ang] = 0.0

                    if not mic_id_substr in global_entry_dict:
                        # sxprint("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from cter_entry_list " % (mic_id_substr))
                        global_entry_dict[mic_id_substr] = {}
                    assert mic_id_substr in global_entry_dict
                    global_entry_dict[mic_id_substr][subkey_cter_entry] = cter_entry
                assert global_entry_dict

            del cter_entry_list  # Do not need this anymore

            # --------------------------------------------------------------------------------
            # Clean up variables related to registration to the global entry dictionary
            # --------------------------------------------------------------------------------
        del mic_basename_tokens
        del mic_id_substr_head_idx
        del coords_pattern_tokens
        del coords_id_substr_head_idx

        # --------------------------------------------------------------------------------
        # Create the list containing only valid micrograph id substrings
        # --------------------------------------------------------------------------------
        # Prepare lists to keep track of invalid (rejected) micrographs
        no_input_mic_id_substr_list = []
        no_coords_mic_id_substr_list = []
        no_cter_entry_mic_id_substr_list = []

        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Checking consistency of the provided dataset ...")

        # Loop over substring id list
        for mic_id_substr in global_entry_dict:
            mic_id_entry = global_entry_dict[mic_id_substr]

            warninig_messages = []
            # selected micrograph basename must have been registed always .
            if subkey_selected_mic_basename in mic_id_entry:
                # Check if associated input micrograph exists
                if not subkey_input_mic_path in mic_id_entry:
                    input_mic_path = mic_pattern.replace("*", mic_id_substr)
                    warninig_messages.append(
                        "    associated input micrograph %s." % (input_mic_path)
                    )
                    no_input_mic_id_substr_list.append(mic_id_substr)

                    # Check if associated coordinate file exists
                if not subkey_coords_path in mic_id_entry:
                    coords_path = coords_pattern.replace("*", mic_id_substr)
                    warninig_messages.append(
                        "    associated coordinates file %s." % (coords_path)
                    )
                    no_coords_mic_id_substr_list.append(mic_id_substr)

                if use_real_ctf_params:
                    # Check if associated CTER partres entry exists
                    if not subkey_cter_entry in mic_id_entry:
                        mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
                        warninig_messages.append(
                            "    associated entry with %s in the CTER partres file %s."
                            % (mic_basename, ctf_params_src)
                        )
                        no_cter_entry_mic_id_substr_list.append(mic_id_substr)
                else:
                    assert not use_real_ctf_params
                    # Register dummy CTF parameters
                    # Create dummy CTER partres entry
                    assert float(ctf_params_src) > 0.0
                    dummy_cter_entry = [None] * n_idx_cter
                    assert len(dummy_cter_entry) == n_idx_cter
                    dummy_cter_entry[idx_cter_def] = 0.0
                    dummy_cter_entry[idx_cter_cs] = 0.0
                    dummy_cter_entry[idx_cter_vol] = 300.0
                    dummy_cter_entry[idx_cter_apix] = float(ctf_params_src)
                    dummy_cter_entry[idx_cter_bfactor] = 0.0
                    dummy_cter_entry[idx_cter_total_ac] = 100.0
                    dummy_cter_entry[idx_cter_astig_amp] = 0.0
                    dummy_cter_entry[idx_cter_astig_ang] = 0.0
                    dummy_cter_entry[idx_cter_sd_def] = 0.0
                    dummy_cter_entry[idx_cter_sd_total_ac] = 0.0
                    dummy_cter_entry[idx_cter_sd_astig_amp] = 0.0
                    dummy_cter_entry[idx_cter_sd_astig_ang] = 0.0
                    dummy_cter_entry[idx_cter_cv_def] = 0.0
                    dummy_cter_entry[idx_cter_cv_astig_amp] = 0.0
                    dummy_cter_entry[idx_cter_error_def] = old_div(
                        0.5, dummy_cter_entry[idx_cter_apix]
                    )  # Set to Nyquist frequency
                    dummy_cter_entry[idx_cter_error_astig] = old_div(
                        0.5, dummy_cter_entry[idx_cter_apix]
                    )  # Set to Nyquist frequency
                    dummy_cter_entry[idx_cter_error_ctf] = old_div(
                        0.5, dummy_cter_entry[idx_cter_apix]
                    )  # Set to Nyquist frequency
                    dummy_cter_entry[idx_cter_max_freq] = old_div(
                        0.5, dummy_cter_entry[idx_cter_apix]
                    )  # Set to Nyquist frequency
                    dummy_cter_entry[idx_cter_reserved] = 0.0
                    dummy_cter_entry[idx_cter_const_ac] = 0.0
                    dummy_cter_entry[idx_cter_phase_shift] = 0.0
                    dummy_cter_entry[idx_cter_mic_name] = ""

                    assert not subkey_cter_entry in mic_id_entry
                    global_entry_dict[mic_id_substr][
                        subkey_cter_entry
                    ] = dummy_cter_entry

                if warninig_messages:
                    sp_global_def.sxprint(
                        "WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr)
                    )
                    for warning_message in warninig_messages:
                        sp_global_def.sxprint(warning_message)
                    sp_global_def.sxprint("    Ignores this as an invalid entry.")
                else:
                    valid_mic_id_substr_list.append(mic_id_substr)

        # Check the input dataset consistency and save the result to a text file, if necessary.
        if options.check_consistency:
            # Create output directory
            assert not os.path.exists(root_out_dir)
            os.makedirs(root_out_dir)
            # Open the consistency check file
            mic_consistency_check_info_path = os.path.join(
                root_out_dir,
                "mic_consistency_check_info_%s.txt" % (get_time_stamp_suffix()),
            )
            sp_global_def.sxprint(" ")
            sp_global_def.sxprint(
                "Generating consistency report of the provided dataset in %s..."
                % (mic_consistency_check_info_path)
            )
            mic_consistency_check_info_file = open(mic_consistency_check_info_path, "w")
            mic_consistency_check_info_file.write(
                "# The consistency information about micrograph IDs that might have problmes with consistency among the provided dataset.\n"
            )
            mic_consistency_check_info_file.write("# \n")

            # Loop over substring id list
            for mic_id_substr in global_entry_dict:
                mic_id_entry = global_entry_dict[mic_id_substr]

                consistency_messages = []
                # Check if associated input micrograph path exists
                if not subkey_input_mic_path in mic_id_entry:
                    input_mic_path = mic_pattern.replace("*", mic_id_substr)
                    consistency_messages.append(
                        "    associated input micrograph %s." % (input_mic_path)
                    )

                    # Check if associated selected micrograph basename exists
                if not subkey_selected_mic_basename in mic_id_entry:
                    input_mic_path = mic_pattern.replace("*", mic_id_substr)
                    consistency_messages.append(
                        "    associated selected micrograph %s." % (input_mic_path)
                    )

                    # Check if associated coordinate file exists
                if not subkey_coords_path in mic_id_entry:
                    coords_path = coords_pattern.replace("*", mic_id_substr)
                    consistency_messages.append(
                        "    associated coordinates file %s." % (coords_path)
                    )

                if use_real_ctf_params:
                    # Check if associated CTER partres entry exists
                    if not subkey_cter_entry in mic_id_entry:
                        mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
                        consistency_messages.append(
                            "    associated entry with %s in the CTER partres file %s."
                            % (mic_basename, ctf_params_src)
                        )
                else:
                    assert not use_real_ctf_params
                    # All entry must have dummy cter entry
                    assert subkey_cter_entry in mic_id_entry

                if consistency_messages:
                    mic_consistency_check_info_file.write(
                        "Micrograph ID %s might have problems with consistency among the provided dataset:\n"
                        % (mic_id_substr)
                    )

                    for consistency_message in consistency_messages:
                        mic_consistency_check_info_file.write(consistency_message)
                        mic_consistency_check_info_file.write("\n")

            # Close the consistency check file, if necessary
            mic_consistency_check_info_file.flush()
            mic_consistency_check_info_file.close()

        # Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
        # we need sort the valid_mic_id_substr_list here
        valid_mic_id_substr_list.sort()

        # --------------------------------------------------------------------------------
        # Print out the summary of input consistency
        # --------------------------------------------------------------------------------
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Summary of dataset consistency check...")
        sp_global_def.sxprint(
            "Detected                           : %6d" % (len(global_entry_dict))
        )
        sp_global_def.sxprint(
            "Valid                              : %6d" % (len(valid_mic_id_substr_list))
        )
        sp_global_def.sxprint(
            "Rejected by no coordinates file    : %6d"
            % (len(no_coords_mic_id_substr_list))
        )

        if use_real_ctf_params:
            sp_global_def.sxprint(
                "Rejected by no CTER partres entry  : %6d"
                % (len(no_cter_entry_mic_id_substr_list))
            )

        # --------------------------------------------------------------------------------
        # Clean up variables related to tracking of invalid (rejected) micrographs
        # --------------------------------------------------------------------------------
        del no_input_mic_id_substr_list
        del no_coords_mic_id_substr_list
        del no_cter_entry_mic_id_substr_list

        # --------------------------------------------------------------------------------
        # Check MPI error condition
        # --------------------------------------------------------------------------------
        # if error_status is None and len(valid_mic_id_substr_list) < n_mpi_procs:
        # 	error_status = ( "Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid micrographs that satisfy all criteria to be processed)."
        # 					 % (n_mpi_procs, len(valid_mic_id_substr_list)), getframeinfo(currentframe()) )
        # 	break
        break

    sp_utilities.if_error_then_all_processes_exit_program(
        error_status
    )  # here three should just be thrown an exception

    # ====================================================================================
    # Obtain the list of micrograph id sustrings
    # ====================================================================================
    # --------------------------------------------------------------------------------
    # Prepare variables for this section
    # --------------------------------------------------------------------------------
    # Prepare variables related to options
    box_size = options.box_size
    box_half = old_div(box_size, 2)

    resample_ratio = options.resample_ratio

    # Prepare the function for reading coordinates files with the specified format.
    # This way, the following switch statement is unnecessary inside of the coordinates loop.
    coords_format = options.coordinates_format.lower()
    read_coords_file = None

    if coords_format == "sphire":
        read_coords_file = read_sphire_coords_file
    elif coords_format == "eman1":
        read_coords_file = read_eman1_coords_file
    elif coords_format == "cryolo":
        read_coords_file = read_eman1_coords_file
    elif coords_format == "eman2":
        read_coords_file = read_eman2_coords_file
    elif coords_format == "spider":
        read_coords_file = read_spider_coords_file
    elif coords_format == "cryolo_helical_segmented":
        read_coords_file = read_cryolo_helical_segmented_coords_file
    else:
        assert False  # Unreachable code
    assert read_coords_file != None
    is_filament = False
    if read_coords_file != read_cryolo_helical_segmented_coords_file:
        mask2d = sp_utilities.model_circle(
            old_div(box_size, 2), box_size, box_size
        )  # Create circular 2D mask to Util.infomask of particle images
    else:
        is_filament = True

        # Preapre variables related to CTF limit option
    abs_ctf_limit_histogram = (
        []
    )  # compute the histogram for micrographs cut of by cter_entry limit.

    # Micrograph baseroot pattern (extension are removed from micrograph basename pattern)
    # for substack file names
    mic_baseroot_pattern = os.path.splitext(mic_basename_pattern)[0]

    # Prepare the counters for the global summary of micrographs
    n_mic_process = 0
    n_mic_reject_no_coords_entry = 0
    n_global_coords_detect = 0
    n_global_coords_process = 0
    n_global_coords_reject_out_of_boundary = 0

    # keep a copy of the root output directory where the final bdb will be created
    unsliced_valid_serial_id_list = valid_mic_id_substr_list
    mpi_proc_dir = root_out_dir
    if RUNNING_UNDER_MPI:
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        # All mpi processes should know global entry directory and valid micrograph id substring list
        global_entry_dict = sp_utilities.wrap_mpi_bcast(
            global_entry_dict, main_mpi_proc
        )
        valid_mic_id_substr_list = sp_utilities.wrap_mpi_bcast(
            valid_mic_id_substr_list, main_mpi_proc
        )

        # Slice the list of valid micrograph id substrings for this mpi process
        mic_start, mic_end = sp_applications.MPI_start_end(
            len(valid_mic_id_substr_list), n_mpi_procs, my_mpi_proc_id
        )
        valid_mic_id_substr_list = valid_mic_id_substr_list[mic_start:mic_end]

        # generate subdirectories user root_out_dir, one for each process
        mpi_proc_dir = os.path.join(root_out_dir, "mpi_proc_%03d" % my_mpi_proc_id)

        # Set up progress message and all necessary output directories
    reject_out_of_boundary_dir = os.path.join(root_out_dir, "reject_out_of_boundary")
    if my_mpi_proc_id == main_mpi_proc:
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint(
            "Micrographs processed by main process (including percent of progress):"
        )
        progress_percent_step = old_div(
            len(valid_mic_id_substr_list), 100.0
        )  # the number of micrograms for main node divided by 100

        # Create output directory
        #
        # NOTE: Toshio Moriya 2017/12/12
        # This might not be necessary since particle_img.write_image() will automatically create all directory tree necessary to save the file.
        # However, it is side-effect of the function, so we will explicitly make root output directory here.
        #
        if not os.path.exists(root_out_dir):
            os.makedirs(root_out_dir)
        sp_global_def.write_command(root_out_dir)
        assert not os.path.exists(reject_out_of_boundary_dir), "MRK_DEBUG"
        os.mkdir(reject_out_of_boundary_dir)

    if RUNNING_UNDER_MPI:
        mpi.mpi_barrier(
            mpi.MPI_COMM_WORLD
        )  # all MPI processes should wait until the directory is created by main process
        #
        # NOTE: 2017/12/12 Toshio Moriya
        # To walk-around synchronisation problem between all MPI nodes and a file server,
        # we use exception to assert the existence of directory.
        #
        try:
            os.mkdir(reject_out_of_boundary_dir)
            assert False, "MRK_DEBUG: Unreachable code..."
        except OSError as err:
            pass
    else:
        assert os.path.exists(reject_out_of_boundary_dir), "MRK_DEBUG"

        # ------------------------------------------------------------------------------------
        # Starting main parallel execution
        # ------------------------------------------------------------------------------------
    for mic_id_substr_idx, mic_id_substr in enumerate(valid_mic_id_substr_list):

        # --------------------------------------------------------------------------------
        # Print out progress if necessary
        # --------------------------------------------------------------------------------
        mic_basename = global_entry_dict[mic_id_substr][subkey_selected_mic_basename]
        assert mic_basename == mic_basename_pattern.replace("*", mic_id_substr)
        if my_mpi_proc_id == main_mpi_proc:
            sp_global_def.sxprint(
                "%s ---> % 2.2f%%"
                % (mic_basename, old_div(mic_id_substr_idx, progress_percent_step))
            )

            # --------------------------------------------------------------------------------
            # Read the associated coordinates according to the specified format and
            # make the coordinates the center of particle image if necessary
            # Do this first because error might happen
            # --------------------------------------------------------------------------------
        coords_path = global_entry_dict[mic_id_substr][subkey_coords_path]
        assert os.path.exists(coords_path)
        assert read_coords_file != None
        coords_list = read_coords_file(coords_path)
        if not coords_list:
            sp_global_def.sxprint(
                "For %s, the associate coordinates file %s does not contain any entries. Skipping..."
                % (mic_basename, coords_path)
            )
            n_mic_reject_no_coords_entry += 1
            continue

            # --------------------------------------------------------------------------------
            # Get CTF parameter if necessary
            # Calculate the resampled pixel size and store it to the cter_entry if necessary
            # Do before expensive micrograph processing
            # --------------------------------------------------------------------------------
        cter_entry = global_entry_dict[mic_id_substr][subkey_cter_entry]
        src_pixel_size = cter_entry[idx_cter_apix]
        if resample_ratio < 1.0:
            assert resample_ratio > 0.0
            # store the resampled pixel size to the cter_entry to generate CTF object of this micrograph
            cter_entry[idx_cter_apix] = old_div(src_pixel_size, resample_ratio)
            if my_mpi_proc_id == main_mpi_proc:
                sp_global_def.sxprint(
                    "Resample micrograph to pixel size %6.4f [A/Pixels] from %6.4f [A/Pixels] and window segments from resampled micrograph."
                    % (cter_entry[idx_cter_apix], src_pixel_size)
                )
                # else:
                # 	assert (resample_ratio == 1.0)
                # 	# Do nothing

                # Generate CTF object of this micrograph
                # indexes 0 to 7 (idx_cter_def to idx_cter_astig_ang) must be same in cter format & ctf object format.
                # ctf_obj = generate_ctf(cter_entry)
                #
                # NOTE: 2017/03/07 Toshio Moriya
                # Due to the change of error handling in generate_ctf()
                # the argument have to be a list with length of 6 or 8 now.
                #
        ctf_entry = []
        ctf_entry.append(cter_entry[idx_cter_def])
        ctf_entry.append(cter_entry[idx_cter_cs])
        ctf_entry.append(cter_entry[idx_cter_vol])
        ctf_entry.append(cter_entry[idx_cter_apix])
        ctf_entry.append(cter_entry[idx_cter_bfactor])
        ctf_entry.append(cter_entry[idx_cter_total_ac])
        ctf_entry.append(cter_entry[idx_cter_astig_amp])
        ctf_entry.append(cter_entry[idx_cter_astig_ang])
        assert len(ctf_entry) == 8
        ctf_obj = sp_utilities.generate_ctf(ctf_entry)

        # --------------------------------------------------------------------------------
        # Read micrograph
        # --------------------------------------------------------------------------------
        mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
        assert mic_path == mic_pattern.replace("*", mic_id_substr)
        try:
            mic_img = sp_utilities.get_im(mic_path)
        except:
            sp_global_def.sxprint(
                "Failed to read the associate micrograph %s for %s. The file might be corrupted. Skipping..."
                % (mic_path, mic_basename)
            )
            continue

            # --------------------------------------------------------------------------------
            # Move to the Fourier space processing
            # --------------------------------------------------------------------------------
        sp_fundamentals.fftip(mic_img)  # In-place fft

        # --------------------------------------------------------------------------------
        # If necessary, apply the hyperbolic tangent low-pass Fourier filter based on the (resampled) CTF limit;
        # Cut off frequency components higher than the (resampled) CTF limit.
        # --------------------------------------------------------------------------------
        if options.limit_ctf:
            # Comput absolute frequency of CTF limit (abs_ctf_limit) with the resampled pixel size
            abs_ctf_limit, angstrom_ctf_limit = sp_morphology.ctflimit(
                box_size,
                cter_entry[idx_cter_def],
                cter_entry[idx_cter_cs],
                cter_entry[idx_cter_vol],
                cter_entry[idx_cter_apix],
            )

            # Adjust the CTF limit according to the resampling ratio and box size
            if resample_ratio < 1.0:
                assert resample_ratio > 0.0
                abs_ctf_limit = old_div(resample_ratio * abs_ctf_limit, float(box_size))
            else:
                assert (
                    resample_ratio == 1.0
                )  # -> src_pixel_size == resampled_pixel_size -> src_pixel_size / resampled_pixel_size == 1.0
                abs_ctf_limit = old_div(abs_ctf_limit, float(box_size))

                # If ctf limit is lower than Nyquist frequency, apply the low pass filter with the cutoff at CTF limit frequencye.
            if abs_ctf_limit < 0.5:
                mic_img = sp_filter.filt_tanl(mic_img, abs_ctf_limit, 0.01)
                abs_ctf_limit_histogram.append(abs_ctf_limit)

                # --------------------------------------------------------------------------------
                # Apply the Gaussian high-pass Fourier filter to micrograph based on the (resampled) box size;
                # Cut off frequency components lower than one that the (resampled) box size can express.
                # Then, move back to the real space processing
                # --------------------------------------------------------------------------------
        mic_img = sp_fundamentals.fft(
            sp_filter.filt_gaussh(mic_img, old_div(resample_ratio, box_size))
        )

        # --------------------------------------------------------------------------------
        # Resample micrograph, map coordinates, and window segments from resampled micrograph using new coordinates
        # after resampling by resample_ratio, resampled pixel size = src_pixel_size/resample_ratio
        # --------------------------------------------------------------------------------
        # NOTE: 2015/04/13 Toshio Moriya
        # resample() efficiently takes care of the case resample_ratio = 1.0 but
        # it does not set apix_*. Even though it sets apix_* when resample_ratio < 1.0...
        mic_img = sp_fundamentals.resample(mic_img, resample_ratio)

        # --------------------------------------------------------------------------------
        # If necessary, invert image contrast of this micrograph
        # --------------------------------------------------------------------------------
        if not options.skip_invert:
            mic_stats = EMAN2_cppwrap.Util.infomask(
                mic_img, None, True
            )  # mic_stat[0:mean, 1:SD, 2:min, 3:max]
            EMAN2_cppwrap.Util.mul_scalar(mic_img, -1.0)
            mic_img += 2 * mic_stats[0]

            # --------------------------------------------------------------------------------
            # Generate the output file path of particle stack for this mpi process
            # --------------------------------------------------------------------------------
        mic_baseroot = mic_baseroot_pattern.replace("*", mic_id_substr)
        local_stack_path = "bdb:%s#" % mpi_proc_dir + mic_baseroot + "_ptcls"
        local_mrcs_name = mic_baseroot + "_ptcls.mrcs"
        local_mrcs_path = os.path.join(mpi_proc_dir, local_mrcs_name)

        local_bdb_stack = EMAN2db.db_open_dict(local_stack_path)
        # --------------------------------------------------------------------------------
        # Prepare coordinates loop variables
        # --------------------------------------------------------------------------------
        nx = mic_img.get_xsize()
        ny = mic_img.get_ysize()
        x0 = old_div(nx, 2)
        y0 = old_div(ny, 2)

        coords_reject_out_of_boundary_messages = []

        # Loop through coordinates
        coords_accepted = []
        idx_info = 0
        idx_id = 1
        for coords_id in range(len(coords_list)):
            # Get coordinates
            x = int(coords_list[coords_id][0])
            y = int(coords_list[coords_id][1])

            # Rescale coordinates if necessary
            if resample_ratio < 1.0:
                assert resample_ratio > 0.0
                x = int(x * resample_ratio)
                y = int(y * resample_ratio)
            else:
                assert resample_ratio == 1.0

                # Window a particle at this coordinates
            if (
                (0 <= x - box_half)
                and (x + box_half <= nx)
                and (0 <= y - box_half)
                and (y + box_half <= ny)
            ):
                coords_accepted.append([None, None])
                coords_accepted[-1][idx_info] = [
                    mic_img,
                    box_size,
                    box_size,
                    1,
                    x - x0,
                    y - y0,
                ]
                coords_accepted[-1][idx_id] = coords_id
            else:
                coords_reject_out_of_boundary_messages.append(
                    "coordinates ID = %04d: x = %4d, y = %4d, box_size = %4d "
                    % (coords_id, x, y, box_size)
                )
                continue

        local_particle_id = 0  # can be different from coordinates_id
        if coords_accepted:
            local_mrcs = EMAN2_cppwrap.EMData(box_size, box_size, len(coords_accepted))
            local_mrcs.set_attr(
                "apix_x", 1.0
            )  # particle_img.set_attr("apix_x", resampled_pixel_size)
            local_mrcs.set_attr(
                "apix_y", 1.0
            )  # particle_img.set_attr("apix_y", resampled_pixel_size)
            local_mrcs.set_attr(
                "apix_z", 1.0
            )  # particle_img.set_attr("apix_z", resampled_pixel_size)
            local_mrcs.set_attr(
                "ptcl_source_apix", src_pixel_size
            )  # Store the original pixel size
            for coords_id, entry in enumerate(coords_accepted):

                original_id = entry[idx_id]
                particle_img = EMAN2_cppwrap.Util.window(*entry[idx_info])

                if is_filament:
                    long_axis_radius = old_div(
                        numpy.int(numpy.sqrt(2 * box_size ** 2)), 2
                    )
                    if options.filament_width != -1:
                        short_axis_radius = old_div(options.filament_width, 2)
                    else:
                        short_axis_radius = old_div(box_size, 6)
                    mask2d = sp_utilities.model_rotated_rectangle2D(
                        long_axis_radius,
                        short_axis_radius,
                        box_size,
                        box_size,
                        coords_list[coords_id][4],
                    )

                    # Normalize this particle image
                particle_img = sp_fundamentals.ramp(particle_img)
                particle_stats = EMAN2_cppwrap.Util.infomask(
                    particle_img, mask2d, False
                )  # particle_stats[0:mean, 1:SD, 2:min, 3:max]
                particle_img -= particle_stats[0]
                try:
                    particle_img = old_div(particle_img, particle_stats[1])
                except ZeroDivisionError:
                    sp_global_def.sxprint(
                        "The standard deviation of the particle image associated with %dth coordinates entry in micrograph %s for %s is zero unexpectedly. Skipping..."
                        % (coords_id, mic_path, mic_basename)
                    )
                    continue

                    # Set header entries (attributes) of this particle image
                    #
                    # NOTE: 2015/04/09 Toshio Moriya
                    # ptcl_source_image might be redundant information...
                    # Consider re-organizing header entries...
                    #
                particle_img_dict = particle_img.get_attr_dict()
                particle_img_dict["ptcl_source_image"] = mic_path
                particle_img_dict["ptcl_source_coord"] = [
                    int(coords_list[original_id][0]),
                    int(coords_list[original_id][1]),
                ]
                particle_img_dict["ptcl_source_coord_id"] = coords_id
                particle_img_dict["ptcl_source_box_id"] = original_id
                try:
                    particle_img_dict["filament_id"] = coords_list[original_id][
                        2
                    ].format(mic_path)
                except IndexError:
                    pass
                try:
                    particle_img_dict["segment_id"] = coords_list[original_id][3]
                except IndexError:
                    pass
                try:
                    particle_img_dict["segment_angle"] = coords_list[original_id][4]
                except IndexError:
                    pass
                try:
                    particle_img_dict["filament_track_length"] = coords_list[
                        original_id
                    ][5]
                except IndexError:
                    pass
                particle_img_dict[
                    "data_n"
                ] = (
                    coords_id
                )  # NOTE: Toshio Moriya 2017/11/20: same as ptcl_source_coord_id but the other program uses this header entry key...
                particle_img_dict["data_path"] = "../" + local_mrcs_name
                particle_img_dict["resample_ratio"] = resample_ratio

                particle_img_dict["nx"] = box_size
                particle_img_dict["ny"] = box_size
                particle_img_dict["nz"] = 1
                #
                # NOTE: 2015/04/13 Toshio Moriya
                # Pawel Comment: Micrograph is not supposed to have CTF header info.
                # So, let's assume it does not exist & ignore its presence.
                # assert (not particle_img.has_ctff())
                #
                # NOTE: 2015/04/13 Toshio Moriya
                # Note that resample() "correctly" updates pixel size of CTF header info if it exists
                #
                # NOTE: 2015/04/13 Toshio Moriya
                # apix_* attributes are updated by resample() only when resample_ratio != 1.0
                # Let's make sure header info is consistent by setting apix_* = 1.0
                # regardless of options, so it is not passed down the processing line
                #
                particle_img_dict[
                    "apix_x"
                ] = 1.0  # particle_img.set_attr("apix_x", resampled_pixel_size)
                particle_img_dict[
                    "apix_y"
                ] = 1.0  # particle_img.set_attr("apix_y", resampled_pixel_size)
                particle_img_dict[
                    "apix_z"
                ] = 1.0  # particle_img.set_attr("apix_z", resampled_pixel_size)
                particle_img_dict[
                    "ptcl_source_apix"
                ] = src_pixel_size  # Store the original pixel size
                particle_img_dict["ctf"] = ctf_obj
                particle_img_dict["ctf_applied"] = 0

                # Write the particle image to local stack file
                # sxprint("MRK_DEBUG: local_stack_path, local_particle_id", local_stack_path, local_particle_id)

                local_bdb_stack[local_particle_id] = particle_img_dict
                local_mrcs.insert_clip(particle_img, (0, 0, local_particle_id))
                # particle_img.write_image(local_stack_path, local_particle_id) # NOTE: Insert slice instead of write the image
                local_particle_id += 1

            local_mrcs.write_image(local_mrcs_path)

            # Save the message list of rejected coordinates because of out-of-boundary
            # sxprint("MRK_DEBUG: len(coords_reject_out_of_boundary_messages) := %d" % len(coords_reject_out_of_boundary_messages))
        if coords_reject_out_of_boundary_messages:
            # Open file path to save the message list
            #
            # NOTE: 2017/12/12 Toshio Moriya
            # To walk-around synchronisation problem between all MPI nodes and a file server,
            # we use exception to assert the existence of directory.
            #
            if RUNNING_UNDER_MPI:
                try:
                    os.mkdir(reject_out_of_boundary_dir)
                    assert False, "MRK_DEBUG: Unreachable code..."
                except OSError as err:
                    pass
            else:
                assert os.path.exists(reject_out_of_boundary_dir), "MRK_DEBUG"
            coords_reject_out_of_boundary_path = os.path.join(
                reject_out_of_boundary_dir,
                os.path.splitext(os.path.basename(coords_path))[0]
                + "_reject_out_of_boundary.txt",
            )
            # sxprint("MRK_DEBUG: coords_reject_out_of_boundary_path := %s" % coords_reject_out_of_boundary_path)
            coords_reject_out_of_boundary_file = open(
                coords_reject_out_of_boundary_path, "w"
            )

            for (
                coords_reject_out_of_boundary_message
            ) in coords_reject_out_of_boundary_messages:
                coords_reject_out_of_boundary_file.write(
                    coords_reject_out_of_boundary_message
                )
                coords_reject_out_of_boundary_file.write("\n")

                # Close the consistency check file, if necessary
            coords_reject_out_of_boundary_file.flush()
            coords_reject_out_of_boundary_file.close()

            # Update the counters for the global summary of micrographs
        n_mic_process += 1
        n_global_coords_detect += len(coords_list)
        n_global_coords_process += local_particle_id
        n_global_coords_reject_out_of_boundary += len(
            coords_reject_out_of_boundary_messages
        )

        # Release the data base of local stack from this process
        # so that the subprocess can access to the data base
        EMAN2db.db_close_dict(local_stack_path)

        # ------------------------------------------------------------------------------------
        # Print out CTF limit information
        # ------------------------------------------------------------------------------------
    if options.limit_ctf:
        if RUNNING_UNDER_MPI:
            abs_ctf_limit_histogram = sp_utilities.wrap_mpi_gatherv(
                abs_ctf_limit_histogram, main_mpi_proc
            )

        if my_mpi_proc_id == main_mpi_proc:
            # Print out the summary of CTF limit absolute frequency
            sp_global_def.sxprint(" ")
            sp_global_def.sxprint(
                "Global summary of CTF limit absolute frequency (--limit_ctf)..."
            )
            sp_global_def.sxprint(
                "Percentage of filtered micrographs: %8.2f"
                % (
                    old_div(
                        len(abs_ctf_limit_histogram) * 100.0,
                        len(unsliced_valid_serial_id_list),
                    )
                )
            )
            sp_global_def.sxprint(" ")

            n_bins = 10
            if len(abs_ctf_limit_histogram) >= n_bins:

                cutoff_region, cutoff_counts = sp_statistics.hist_list(
                    abs_ctf_limit_histogram, n_bins
                )
                sp_global_def.sxprint(
                    "Histogram of CTF limit absolute frequency used for the filtering:"
                )
                sp_global_def.sxprint("      CTF limit       counts")
                for bin_id in range(n_bins):
                    sp_global_def.sxprint(
                        " %14.7f     %7d"
                        % (cutoff_region[bin_id], cutoff_counts[bin_id])
                    )
            else:
                sp_global_def.sxprint(
                    "The number of filtered micrographs (%d) is less than the number of bins (%d). No histogram is produced."
                    % (len(abs_ctf_limit_histogram), n_bins)
                )

                # ------------------------------------------------------------------------------------
                # Print out summary of processing
                # ------------------------------------------------------------------------------------
    if RUNNING_UNDER_MPI:
        n_mic_process = mpi.mpi_reduce(
            n_mic_process,
            1,
            mpi.MPI_INT,
            mpi.MPI_SUM,
            main_mpi_proc,
            mpi.MPI_COMM_WORLD,
        )
        n_mic_reject_no_coords_entry = mpi.mpi_reduce(
            n_mic_reject_no_coords_entry,
            1,
            mpi.MPI_INT,
            mpi.MPI_SUM,
            main_mpi_proc,
            mpi.MPI_COMM_WORLD,
        )
        n_global_coords_detect = mpi.mpi_reduce(
            n_global_coords_detect,
            1,
            mpi.MPI_INT,
            mpi.MPI_SUM,
            main_mpi_proc,
            mpi.MPI_COMM_WORLD,
        )
        n_global_coords_process = mpi.mpi_reduce(
            n_global_coords_process,
            1,
            mpi.MPI_INT,
            mpi.MPI_SUM,
            main_mpi_proc,
            mpi.MPI_COMM_WORLD,
        )
        n_global_coords_reject_out_of_boundary = mpi.mpi_reduce(
            n_global_coords_reject_out_of_boundary,
            1,
            mpi.MPI_INT,
            mpi.MPI_SUM,
            main_mpi_proc,
            mpi.MPI_COMM_WORLD,
        )

        # Print out the summary of all micrographs
    if main_mpi_proc == my_mpi_proc_id:
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Summary of micrograph level processing...")
        sp_global_def.sxprint(
            "Valid                              : %6d"
            % (len(unsliced_valid_serial_id_list))
        )
        sp_global_def.sxprint(
            "Processed                          : %6d" % (n_mic_process)
        )
        sp_global_def.sxprint(
            "Rejected by no coordinates entries : %6d" % (n_mic_reject_no_coords_entry)
        )
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Global summary of coordinates level processing...")
        sp_global_def.sxprint(
            "Detected                           : %6d" % (n_global_coords_detect)
        )
        sp_global_def.sxprint(
            "Processed                          : %6d" % (n_global_coords_process)
        )
        sp_global_def.sxprint(
            "Rejected by out of boundary        : %6d"
            % (n_global_coords_reject_out_of_boundary)
        )
        if n_global_coords_reject_out_of_boundary > 0:
            assert os.path.exists(reject_out_of_boundary_dir), "MRK_DEBUG"
            sp_global_def.sxprint(
                "    NOTE: Information of rejected coordinates by out of boundary are saved in %s files."
                % (
                    os.path.join(
                        reject_out_of_boundary_dir,
                        os.path.splitext(os.path.basename(coords_pattern))[0]
                        + "_reject_out_of_boundary.txt",
                    )
                )
            )
        else:
            assert n_global_coords_reject_out_of_boundary == 0, "MRK_DEBUG"
            if os.path.exists(reject_out_of_boundary_dir):
                # This directory must be empty and should not contain any sub-directries
                assert os.listdir(reject_out_of_boundary_dir) == [], "MRK_DEBUG"
                shutil.rmtree(reject_out_of_boundary_dir, ignore_errors=True)

    if RUNNING_UNDER_MPI:
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if main_mpi_proc == my_mpi_proc_id:
        # NOTE: Toshio Moriya 2016/10/27
        # Pawel commented out the command execution below because of an MPI issue.
        # When the program is running, all CPUs are writing to the disk.
        # However, at the end of it, there is no guarantee all files will be closed
        # as seen from the main node.
        # The only way to prevent the crash is to execute e2bdb after the program finished,
        # and even with that one is well advised to wait.
        """Multiline Comment0"""

        if RUNNING_UNDER_MPI:
            e2bdb_command = (
                "e2bdb.py  "
                + root_out_dir
                + "/mpi_proc_*  --makevstack=bdb:"
                + root_out_dir
                + "#data"
            )
            sp_global_def.sxprint(" ")
            sp_global_def.sxprint(
                "Please execute from the command line :  ", e2bdb_command
            )
        else:
            e2bdb_command = (
                "e2bdb.py  "
                + root_out_dir
                + "  --makevstack=bdb:"
                + root_out_dir
                + "#data"
            )
            sp_utilities.cmdexecute(e2bdb_command, printing_on_success=False)

        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("DONE!!!")
        sp_global_def.sxprint(" ")

        # ====================================================================================
        # Clean up
        # ====================================================================================
        # ------------------------------------------------------------------------------------
        # Reset SPHIRE global definitions
        # ------------------------------------------------------------------------------------
    sp_global_def.LOGFILE = original_logfilename

    # ------------------------------------------------------------------------------------
    # Clean up MPI related variables
    # ------------------------------------------------------------------------------------
    if RUNNING_UNDER_MPI:
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    sys.stdout.flush()
    return


# ========================================================================================
# Define main function for command line execution
# ========================================================================================
def main():
    RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
    if RUNNING_UNDER_MPI:
        mpi.mpi_init(0, [])  # On OS X, there is an error if MPI is initialized and not finalized, hence the conditional
    sp_global_def.print_timestamp("Start")
    run()
    sp_global_def.print_timestamp("Finish")

    if RUNNING_UNDER_MPI:
        mpi.mpi_finalize()

if __name__ == "__main__":

    main()

# ========================================================================================
#  END OF FILE
# ========================================================================================
