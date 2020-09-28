#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2018-2019 (toshio.moriya@kek.jp)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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


import EMAN2_cppwrap
import EMAN2db
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
from ..libpy import sp_utilities
import sys
import time
from builtins import range




# ========================================================================================
# General Helper Functions
# ========================================================================================
# ----------------------------------------------------------------------------------------
#  Generate command line
# ----------------------------------------------------------------------------------------


def get_cmd_line():
    cmd_line = ""
    for arg in sys.argv:
        cmd_line += arg + "  "
    cmd_line = "Shell line command: " + cmd_line
    return cmd_line


# ----------------------------------------------------------------------------------------
# Get suffix of current time stamp
# ----------------------------------------------------------------------------------------


def get_time_stamp_suffix():
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


# ----------------------------------------------------------------------------------------
# Resample 2D image
# ----------------------------------------------------------------------------------------


def mrk_resample2d(img, sub_rate, target_size=None):
    """
		mrk_resample2d() is based on resample() in sp_fundamentals.py
		However, the following modifications were made
		(1) Do not automatically adjust some values of header entries below
		- ctf              : apix
		- xform.projection : tx, ty
		- xform.align2d    : tx, ty
		(2) Support only 2D images (not 3D volume)
		(3) Input image has to be square (i.e. nx == ny)
		(4) Keep the dimension of output image same as input image
		resample image based on the value of sub_rate.
		sub_rate < 1.0, subsampling the image.
		sub_rate > 1.0, upsampling the image using new gridding interpolation.
		(??? fit_to_fft will change the ouput image size to an fft_friendly size ???)
	"""

    original_size = img.get_xsize()
    assert original_size > 1
    assert img.get_ysize() == original_size
    assert img.get_zsize() == 1

    resample_img = None
    if sub_rate == 1.0:
        resample_img = img.copy()
    elif sub_rate < 1.0:
        resample_img = sp_fundamentals.subsample(img, sub_rate)
    else:
        assert sub_rate > 1.0
        nn = int(original_size * sub_rate + 0.5)
        resample_img, kb = sp_fundamentals.prepi(
            EMAN2_cppwrap.Util.pad(img, nn, nn, 1, 0, 0, 0, "circumference")
        )
        resample_img = resample_img.rot_scale_conv_new(0.0, 0.0, 0.0, kb, sub_rate)
    assert resample_img is not None

    resampled_size = resample_img.get_xsize()
    assert resample_img.get_xsize() == resampled_size
    assert resample_img.get_ysize() == resampled_size
    assert resample_img.get_zsize() == 1

    if target_size is None:
        # In this case, restore the original dimensions
        target_size = original_size
    assert target_size is not None
    assert target_size > 1

    if resampled_size < target_size:
        resample_img = EMAN2_cppwrap.Util.pad(
            resample_img, target_size, target_size, 1, 0, 0, 0, "circumference"
        )
    elif resampled_size > target_size:
        resample_img = EMAN2_cppwrap.Util.window(
            resample_img, target_size, target_size, 1, 0, 0
        )
    else:
        assert resampled_size == target_size
    assert resample_img.get_xsize() == target_size
    assert resample_img.get_ysize() == target_size
    assert resample_img.get_zsize() == 1

    return resample_img


# ========================================================================================
#  Main function
# ========================================================================================


def run():
    program_name = os.path.basename(sys.argv[0])
    usage = (
        program_name
        + """  input_micrograph_pattern  input_rebox_pattern  output_directory  --selection_list=SELECTION_TEXT_FILE  --box_size=BOX_SIZE  --skip_invert  --mic_resample_ratio=RATIO  --swap_ctf_params=CTER_PARTRES_FILE_PATH  --check_consistency

Rewindow particles from micrographs using the information stored in rebox files.

All Micrographs Mode - Process all micrographs in a directory:
	Specify path pattern of input micrographs and rebox files with a wild card (*). 
	Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). 
	The path pattern must be enclosed by single quotes (') or double quotes ("). (Note: sxgui.py automatically adds single quotes (')). 
	The substring at the variable part must be same between a associated pair of input micrograph and rebox file.
	BDB files can not be selected as input micrographs.
	Finally, specify output directory where all outputs should be saved.
	In this mode, all micrographs matching the path pattern will be processed.

	mpirun  -np  32  sxrewindow.py  './mic*.hdf'  'outdir_rebox/centered_rebox/mic*_centered_rebox.rbx'  outdir_rewindow  --box_size=64

	You can also ignore per-particle CTF information stored in rebox files and use CTF information stored in the CTER partres.txt instead.

	mpirun  -np  32  sxrewindow.py  './mic*.hdf'  'outdir_rebox/centered_rebox/mic*_centered_rebox.rbx'  outdir_rewindow  --box_size=64  --swap_ctf_params='outdir_cter/partres.txt'


Selected Micrographs Mode - Process all micrographs in a selection list file:
	In addition to input micrographs path pattern, rebox files path pattern, and output directry arguments, 
	specify a name of micrograph selection list text file using --selection_list option.
	In this mode, only micrographs in the selection list which matches the file name part of the pattern (ignoring the directory paths) will be processed.
	If a micrograph name in the selection list does not exists in the directory specified by the micrograph path pattern, processing of the micrograph will be skipped.

	mpirun  -np  32  sxrewindow.py  './mic*.hdf'  'outdir_rebox/centered_rebox/mic*_centered_rebox.rbx'  outdir_rewindow  --selection_list=mic_list.txt  --box_size=64

Single Micrograph Mode - Process a single micrograph:
	In addition to input micrographs path pattern, rebox files path pattern, CTF paramters source, and output directry arguments, 
	specify a single micrograph name using --selection_list option.
	In this mode, only the specified single micrograph will be processed.
	If this micrograph name does not matches the file name part of the pattern (ignoring the directory paths), the process will exit without processing it.
	If this micrograph name matches the file name part of the pattern but does not exists in the directory which specified by the micrograph path pattern, again the process will exit without processing it.
	Use single processor for this mode. 

	sxrewindow.py  './mic*.hdf'  'outdir_rebox/centered_rebox/mic*_centered_rebox.rbx'  outdir_rewindow  --selection_list=mic0.hdf  --box_size=64

For negative staining data, use --skip_invert.

	mpirun  -np  32  sxrewindow.py  './mic*.hdf'  'outdir_rebox/centered_rebox/mic*_centered_rebox.rbx'  outdir_rewindow  --box_size=64  --skip_invert

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
        "--box_size",
        type="int",
        default=256,
        help="Particle box size [Pixels]: The x and y dimensions of square area to be windowed. The box size after resampling is assumed when mic_resample_ratio < 1.0. (default 256)",
    )
    parser.add_option(
        "--skip_invert",
        action="store_true",
        default=False,
        help="Skip invert image contrast: Use this option for negative staining data. By default, the image contrast is inverted for cryo data. (default False)",
    )
    parser.add_option(
        "--mic_resample_ratio",
        type="float",
        default=1.0,
        help="Image size reduction factor (<1): Use a value between 0.0 and 1.0 (excluding 0.0). The new pixel size will be automatically recalculated and stored in CTF paramers when mic_resample_ratio < 1.0 is used. (default 1.0)",
    )

    parser.add_option(
        "--swap_ctf_params",
        type="string",
        default=None,
        help="Swap CTF parameters: Swaps CTF parameters by setting the CTF parameters in the specified CTER partres file while ignoring the CTF parameters in the input rebox parameters file. Typically, specify the file produced by sxcter and normally called partres.txt. Alternatively, enter pixel size [A/Pixels] to simulate ideal CTF. By default, the program uses the CTF parameters in the input rebox parameters file. (default None)",
    )
    parser.add_option(
        "--check_consistency",
        action="store_true",
        default=False,
        help="Check consistency of dataset: Create a text file containing the list of Micrograph ID entries might have inconsistency among the provided dataset. (i.e. mic_consistency_check_info_TIMESTAMP.txt). (default False)",
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
    rebox_pattern = None
    root_out_dir = None
    ctf_params_src = options.swap_ctf_params
    # Not a real while, each "if" statement has the opportunity to use break when errors need to be reported
    error_status = None
    while True:
        # --------------------------------------------------------------------------------
        # Check the number of arguments. If OK, then prepare variables for them
        # --------------------------------------------------------------------------------
        if error_status is None and len(args) != 3:
            error_status = (
                "Please check usage for number of arguments.\n Usage: "
                + usage
                + "\n"
                + "Please run %s -h for help." % (program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        assert len(args) == 3
        mic_pattern = args[0]
        rebox_pattern = args[1]
        root_out_dir = args[2]

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

        if error_status is None and rebox_pattern.find("*") == -1:
            error_status = (
                "Input rebox file name pattern must contain wild card (*). Please check input_rebox_pattern argument. Run %s -h for help."
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

        if error_status is None and (options.box_size <= 0):
            error_status = (
                "Invalid option value: --box_size=%s. The box size must be an interger larger than zero. Please run %s -h for help."
                % (options.box_size, program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        if error_status is None and (
            options.mic_resample_ratio <= 0.0 or options.mic_resample_ratio > 1.0
        ):
            error_status = (
                "Invalid option value: --mic_resample_ratio=%s. Please run %s -h for help."
                % (options.mic_resample_ratio, program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        if ctf_params_src is not None:
            if not is_float(ctf_params_src):
                assert type(ctf_params_src) is str
                # This should be string for CTER partres (CTF parameter) file
                if error_status is None and os.path.exists(ctf_params_src) == False:
                    error_status = (
                        "Specified CTER partres file is not found. Please check --swap_ctf_params option. Run %s -h for help."
                        % (program_name),
                        inspect.getframeinfo(inspect.currentframe()),
                    )
                    break
            else:
                assert is_float(ctf_params_src)
                if error_status is None and float(ctf_params_src) <= 0.0:
                    error_status = (
                        "Specified pixel size is not larger than 0.0. Please check --swap_ctf_params option. Run %s -h for help."
                        % (program_name),
                        inspect.getframeinfo(inspect.currentframe()),
                    )
                    break

        break
    sp_utilities.if_error_then_all_processes_exit_program(error_status)
    assert mic_pattern != None
    assert rebox_pattern != None
    assert root_out_dir != None

    # ------------------------------------------------------------------------------------
    # Check warning conditions of options
    # ------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------
    # Define enumerators
    # ------------------------------------------------------------------------------------
    # For indices of rebox parameters in rebox file format
    i_enum = -1
    i_enum += 1
    idx_rebox_mic_coord_id = i_enum  # Original box rebox ID
    i_enum += 1
    idx_rebox_mic_coord_x = i_enum  # Original box center x-coordinate
    i_enum += 1
    idx_rebox_mic_coord_y = i_enum  # Original box center y-coordinate
    i_enum += 1
    idx_rebox_mic_resample_ratio = (
        i_enum
    )  # Resample ratio used previous windowing or rewindowing
    i_enum += 1
    idx_rebox_ctf_defocus = i_enum  # Defocus of CTF
    i_enum += 1
    idx_rebox_ctf_cs = i_enum  # CS of CTF
    i_enum += 1
    idx_rebox_ctf_voltage = i_enum  # Voltage of CTF
    i_enum += 1
    idx_rebox_ctf_apix = i_enum  # Pixel size [A] of CTF
    i_enum += 1
    idx_rebox_ctf_bfactor = i_enum  # B-factor of CTF
    i_enum += 1
    idx_rebox_ctf_ampcont = i_enum  # Amplitude contrast of CTF
    i_enum += 1
    idx_rebox_ctf_dfdiff = i_enum  # Defocus difference or astigmatism amplitude of CTF
    i_enum += 1
    idx_rebox_ctf_dfang = i_enum  # Defocus angle or astigmatism angle of CTF
    i_enum += 1
    idx_rebox_proj_phi = (
        i_enum
    )  # Eulerian angle for 3D reconstruction (azimuthal) of projection parameters (xform.projection)
    i_enum += 1
    idx_rebox_proj_theta = (
        i_enum
    )  # Eulerian angle for 3D reconstruction (tilt) of projection parameters (xform.projection)
    i_enum += 1
    idx_rebox_proj_psi = (
        i_enum
    )  # Eulerian angle for 3D reconstruction (in-plane rotation of projection) of projection parameters (xform.projection)
    i_enum += 1
    idx_rebox_proj_sx = (
        i_enum
    )  # Shift in x direction of projection parameters (xform.projection)
    i_enum += 1
    idx_rebox_proj_sy = (
        i_enum
    )  # Shift in y direction of projection parameters (xform.projection)
    i_enum += 1
    idx_rebox_pp_def_error_accum = i_enum  # Accumulated per-particle defocus error
    i_enum += 1
    idx_rebox_pp_mag_error_accum = (
        i_enum
    )  # Accumulated per-particle magnification error
    i_enum += 1
    idx_rebox_chunk_id = i_enum  # NOT SUPPORTED YET (Toshio Moriya 2018/07/10)
    i_enum += 1
    n_idx_rebox = i_enum

    # For indices of CTF parameters
    i_enum = -1
    i_enum += 1
    idx_ctf_defocus = i_enum  # Defocus of CTF
    i_enum += 1
    idx_ctf_cs = i_enum  # CS of CTF
    i_enum += 1
    idx_ctf_voltage = i_enum  # Voltage of CTF
    i_enum += 1
    idx_ctf_apix = i_enum  # Pixel size [A] of CTF
    i_enum += 1
    idx_ctf_bfactor = i_enum  # B-factor of CTF
    i_enum += 1
    idx_ctf_ampcont = i_enum  # Amplitude contrast of CTF
    i_enum += 1
    idx_ctf_dfdiff = i_enum  # Defocus difference or astigmatism amplitude of CTF
    i_enum += 1
    idx_ctf_dfang = i_enum  # Defocus angle or astigmatism angle of CTF
    i_enum += 1
    n_idx_ctf = i_enum
    assert (
        n_idx_ctf == 8
    )  # Due to the error handling in generate_ctf(), the argument has to be a list with length of 6 or 8.

    # For indices of projection parameters
    i_enum = -1
    i_enum += 1
    idx_proj_phi = (
        i_enum
    )  # Eulerian angle for 3D reconstruction (azimuthal) of projection parameters (xform.projection)
    i_enum += 1
    idx_proj_theta = (
        i_enum
    )  # Eulerian angle for 3D reconstruction (tilt) of projection parameters (xform.projection)
    i_enum += 1
    idx_proj_psi = (
        i_enum
    )  # Eulerian angle for 3D reconstruction (in-plane rotation of projection) of projection parameters (xform.projection)
    i_enum += 1
    idx_proj_sx = (
        i_enum
    )  # Shift in x direction of projection parameters (xform.projection)
    i_enum += 1
    idx_proj_sy = (
        i_enum
    )  # Shift in y direction of projection parameters (xform.projection)
    i_enum += 1
    n_idx_proj = i_enum

    if ctf_params_src is not None:
        # NOTE: 2018/04/25 Toshio Moriya
        # sxrewindow does not support the old format of CTER partres file (BEFIRE 2017/12/05).
        #
        # For indies of CTER partres parameters in the new format (AFTER 2017/12/05).
        # All mpi processes must have access to these indices
        i_enum = -1
        i_enum += 1
        idx_cter_def = i_enum  # defocus [um]; index must be same as ctf object format
        i_enum += 1
        idx_cter_cs = i_enum  # Cs [mm]; index must be same as ctf object format
        i_enum += 1
        idx_cter_vol = i_enum  # voltage[kV]; index must be same as ctf object format
        i_enum += 1
        idx_cter_apix = (
            i_enum
        )  # pixel size [A]; index must be same as ctf object format
        i_enum += 1
        idx_cter_bfactor = (
            i_enum
        )  # B-factor [A^2]; index must be same as ctf object format
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
        idx_cter_cv_astig_amp = (
            i_enum
        )  # coefficient of variation of astigmatism amp [%]
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

        # Check the source of CTF parameteres and select the CTF mode
        # (1) Real CTF parameters mode  : Use real CTF parameters stored in a CTER partres (CTF parameter) file for cryo data (use_real_ctf_params is True)
        # (2) Dummy CTF parameters mode : Create dummy CTF parameters for negative staining data (use_real_ctf_params is False)
        use_real_ctf_params = not is_float(ctf_params_src)

    # ------------------------------------------------------------------------------------
    # Prepare the variables for all sections
    # ------------------------------------------------------------------------------------
    # Micrograph basename pattern (directory path is removed from micrograph path pattern)
    mic_basename_pattern = os.path.basename(mic_pattern)

    # Global entry dictionary (all possible entries from all lists) for all mic id substring
    global_entry_dict = {}  # mic id substring is the key
    subkey_input_mic_path = "Input Micrograph Path"
    subkey_selected_mic_basename = "Selected Micrograph Basename"
    subkey_rebox_path = "Input Rebox File Path"
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

        # Prefix and suffix of rebox file path pattern
        # to find the head/tail indices of rebox file id substring
        rebox_pattern_tokens = rebox_pattern.split("*")
        assert len(rebox_pattern_tokens) == 2
        # Find head index of rebox id substring
        rebox_id_substr_head_idx = len(rebox_pattern_tokens[0])

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
        if error_status is None and len(input_mic_path_list) == 0:
            error_status = (
                "No micrograph files are found in the directory specified by micrograph path pattern (%s). Please check input_micrograph_pattern argument. Run %s -h for help."
                % (os.path.dirname(mic_pattern), program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break
        assert len(input_mic_path_list) > 0

        # Register micrograph id substrings to the global entry dictionary
        for input_mic_path in input_mic_path_list:
            # Find tail index of micrograph id substring and extract the substring from the micrograph name
            input_mic_basename = os.path.basename(input_mic_path)
            mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
            mic_id_substr = input_mic_basename[
                mic_id_substr_head_idx:mic_id_substr_tail_idx
            ]
            assert input_mic_path == mic_pattern.replace("*", mic_id_substr)
            if not mic_id_substr in global_entry_dict:
                # sxprint("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
                global_entry_dict[mic_id_substr] = {}
            assert mic_id_substr in global_entry_dict
            global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path
        assert len(global_entry_dict) > 0

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
                    selected_mic_path_list[0], basestring
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
            assert len(selected_mic_path_list) > 0

            selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
            if selected_mic_directory != "":
                sp_global_def.sxprint(
                    "    NOTE: Program disregards the directory paths in the selection list (%s)."
                    % (selected_mic_directory)
                )

        assert len(selected_mic_path_list) > 0

        # Register micrograph id substrings to the global entry dictionary
        for selected_mic_path in selected_mic_path_list:
            # Find tail index of micrograph id substring and extract the substring from the micrograph name
            selected_mic_basename = os.path.basename(selected_mic_path)
            mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
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
        assert len(global_entry_dict) > 0

        del selected_mic_path_list  # Do not need this anymore
        del input_mic_path_list  # Do not need this anymore

        # --------------------------------------------------------------------------------
        # Register rebox id substrings in rebox path list to the global entry dictionary.
        # rebox id substring (rebox_id_substr) and micrograph id substring (mic_id_substr)
        # should be the same for the associated pair of micrograph and coordnates file.
        # --------------------------------------------------------------------------------
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Checking the rebox files...")
        rebox_path_list = glob.glob(rebox_pattern)

        # Check error condition of rebox file path list
        sp_global_def.sxprint(
            "Found %d rebox files in %s directory."
            % (len(rebox_path_list), os.path.dirname(rebox_pattern))
        )
        if error_status is None and len(rebox_path_list) == 0:
            error_status = (
                "No rebox files are found in the directory specified by rebox file path pattern (%s). Please check input_rebox_pattern argument. Run %s -h for help."
                % (os.path.dirname(rebox_pattern), program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break
        assert len(rebox_path_list) > 0

        for rebox_path in rebox_path_list:
            # Find tail index of rebox id substring and extract the substring from the rebox file path
            rebox_id_substr_tail_idx = rebox_path.index(rebox_pattern_tokens[1])
            rebox_id_substr = rebox_path[
                rebox_id_substr_head_idx:rebox_id_substr_tail_idx
            ]
            assert rebox_path == rebox_pattern.replace("*", rebox_id_substr)
            if not rebox_id_substr in global_entry_dict:
                # sxprint("MRK_DEBUG: Added new rebox_id_substr (%s) to global_entry_dict from rebox_path_list " % (rebox_id_substr))
                global_entry_dict[rebox_id_substr] = {}
            assert rebox_id_substr in global_entry_dict
            global_entry_dict[rebox_id_substr][subkey_rebox_path] = rebox_path
        assert len(global_entry_dict) > 0

        del rebox_path_list  # Do not need this anymore

        # --------------------------------------------------------------------------------
        # If necessary, register micrograph id substrings of CTER partres entries to the global entry dictionary
        # --------------------------------------------------------------------------------
        if ctf_params_src is not None:
            if use_real_ctf_params:
                # This should be string for CTER partres (CTF parameter) file
                sp_global_def.sxprint(" ")
                sp_global_def.sxprint(
                    "The program will use the CTF parameters stored in the CTER partres file specified through -swap_ctf_params option, while ignoring the CTF parameters in the input rebox parameters file..."
                )
                sp_global_def.sxprint(
                    "  It will also ignores all defocus errors extracted from the input rebox paramters file and reset it to zero..."
                )
                sp_global_def.sxprint(" ")
                sp_global_def.sxprint("Checking the CTER partres file...")
                assert os.path.exists(ctf_params_src)
                cter_entry_list = sp_utilities.read_text_row(ctf_params_src)

                # Check error condition of CTER partres entry list
                sp_global_def.sxprint(
                    "Found %d CTER partres entries in %s."
                    % (len(cter_entry_list), ctf_params_src)
                )
                if error_status is None and len(cter_entry_list) == 0:
                    error_status = (
                        "No CTER partres entries are found in %s. Please check input_ctf_params_source argument. Run %s -h for help."
                        % (ctf_params_src, program_name),
                        inspect.getframeinfo(inspect.currentframe()),
                    )
                    break
                assert len(cter_entry_list) > 0

                #
                # NOTE: 2017/12/05 Toshio Moriya
                # The following code is to support the old format of CTER partres file. It should be removed near future
                #
                if error_status is None and len(cter_entry_list[0]) != n_idx_cter:
                    error_status = (
                        "The number of columns (%d) has to be %d in %s."
                        % (len(cter_entry_list[0]), n_idx_cter, ctf_params_src),
                        inspect.getframeinfo(inspect.currentframe()),
                    )
                    break
                assert len(cter_entry_list[0]) == n_idx_cter

                # Support only  NEW CTER partres format (AFTER 2017/12/05)
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
                    mic_id_substr_tail_idx = cter_mic_basename.index(
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

                    # if(cter_entry[idx_cter_sd_astig_ang] > options.astigmatism_error):
                    # 	sxprint("    NOTE: Astigmatism angular SD of %s (%f degree) exceeds specified limit (%f degree). Resetting astigmatism parameters to zeros..." % (cter_mic_basename, cter_entry[idx_cter_sd_astig_ang], options.astigmatism_error))
                    # 	cter_entry[idx_cter_astig_amp] = 0.0
                    # 	cter_entry[idx_cter_astig_ang] = 0.0

                    if not mic_id_substr in global_entry_dict:
                        # sxprint("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from cter_entry_list " % (mic_id_substr))
                        global_entry_dict[mic_id_substr] = {}
                    assert mic_id_substr in global_entry_dict
                    global_entry_dict[mic_id_substr][subkey_cter_entry] = cter_entry
                assert len(global_entry_dict) > 0

                del cter_entry_list  # Do not need this anymore
            else:
                assert not use_real_ctf_params
                # This should be string for CTER partres (CTF parameter) file
                sp_global_def.sxprint(" ")
                sp_global_def.sxprint(
                    "The program will use the simulated ideal CTF with the pixel size %f specified through -swap_ctf_params option, while ignoring the CTF parameters in the input rebox parameters file..."
                    % (float(ctf_params_src))
                )
                sp_global_def.sxprint(
                    "  It will also ignores all defocus errors extracted from the input rebox paramters file and reset it to zero..."
                )

        # --------------------------------------------------------------------------------
        # Clean up variables related to registration to the global entry dictionary
        # --------------------------------------------------------------------------------
        del mic_basename_tokens
        del mic_id_substr_head_idx
        del rebox_pattern_tokens
        del rebox_id_substr_head_idx

        # --------------------------------------------------------------------------------
        # Create the list containing only valid micrograph id substrings
        # --------------------------------------------------------------------------------
        # Prepare lists to keep track of invalid (rejected) micrographs
        no_input_mic_id_substr_list = []
        no_rebox_mic_id_substr_list = []
        no_cter_entry_mic_id_substr_list = []

        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Checking consistency of the provided dataset ...")

        # Loop over substring id list
        for mic_id_substr in global_entry_dict:
            mic_id_entry = global_entry_dict[mic_id_substr]

            warinnig_messages = []
            # selected micrograph basename must have been registed always .
            if subkey_selected_mic_basename in mic_id_entry:
                # Check if associated input micrograph exists
                if not subkey_input_mic_path in mic_id_entry:
                    input_mic_path = mic_pattern.replace("*", mic_id_substr)
                    warinnig_messages.append(
                        "    associated input micrograph %s." % (input_mic_path)
                    )
                    no_input_mic_id_substr_list.append(mic_id_substr)

                # Check if associated rebox file exists
                if not subkey_rebox_path in mic_id_entry:
                    rebox_path = rebox_pattern.replace("*", mic_id_substr)
                    warinnig_messages.append(
                        "    associated rebox file %s." % (rebox_path)
                    )
                    no_rebox_mic_id_substr_list.append(mic_id_substr)

                if ctf_params_src is not None:
                    if use_real_ctf_params:
                        # Check if associated CTER partres entry exists
                        if not subkey_cter_entry in mic_id_entry:
                            mic_basename = mic_basename_pattern.replace(
                                "*", mic_id_substr
                            )
                            warinnig_messages.append(
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

                if len(warinnig_messages) > 0:
                    sp_global_def.sxprint(
                        "WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr)
                    )
                    for warinnig_message in warinnig_messages:
                        sp_global_def.sxprint(warinnig_message)
                    sp_global_def.sxprint("    Ignores this as an invalid entry.")
                else:
                    # sxprint("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
                    valid_mic_id_substr_list.append(mic_id_substr)
            # else:
            # 	assert (not subkey_selected_mic_basename in mic_id_entry)
            # 	# This entry is not in the selection list. Do nothing

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

                # Check if associated rebox file exists
                if not subkey_rebox_path in mic_id_entry:
                    rebox_path = rebox_pattern.replace("*", mic_id_substr)
                    consistency_messages.append(
                        "    associated rebox file %s." % (rebox_path)
                    )

                if ctf_params_src is not None:
                    if use_real_ctf_params:
                        # Check if associated CTER partres entry exists
                        if not subkey_cter_entry in mic_id_entry:
                            mic_basename = mic_basename_pattern.replace(
                                "*", mic_id_substr
                            )
                            consistency_messages.append(
                                "    associated entry with %s in the CTER partres file %s."
                                % (mic_basename, ctf_params_src)
                            )
                    else:
                        assert not use_real_ctf_params
                        # All entry must have dummy cter entry
                        assert subkey_cter_entry in mic_id_entry

                if len(consistency_messages) > 0:
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
        # sxprint("MRK_DEBUG: before sort, valid_mic_id_substr_list := ", valid_mic_id_substr_list)
        valid_mic_id_substr_list.sort()
        # sxprint("MRK_DEBUG: after sort, valid_mic_id_substr_list := ", valid_mic_id_substr_list)

        # --------------------------------------------------------------------------------
        # Print out the summary of input consistency
        # --------------------------------------------------------------------------------
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Summary of dataset consistency check...")
        sp_global_def.sxprint(
            "Detected                  : %6d" % (len(global_entry_dict))
        )
        sp_global_def.sxprint(
            "Valid                     : %6d" % (len(valid_mic_id_substr_list))
        )
        sp_global_def.sxprint(
            "Rejected by no rebox file : %6d" % (len(no_rebox_mic_id_substr_list))
        )
        if ctf_params_src is not None:
            if use_real_ctf_params:
                sp_global_def.sxprint(
                    "Rejected by no CTER partres entry  : %6d"
                    % (len(no_cter_entry_mic_id_substr_list))
                )

        # --------------------------------------------------------------------------------
        # Clean up variables related to tracking of invalid (rejected) micrographs
        # --------------------------------------------------------------------------------
        del no_input_mic_id_substr_list
        del no_rebox_mic_id_substr_list
        del no_cter_entry_mic_id_substr_list

        # --------------------------------------------------------------------------------
        # Check MPI error condition
        # --------------------------------------------------------------------------------
        if error_status is None and len(valid_mic_id_substr_list) < n_mpi_procs:
            error_status = (
                "Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid micrographs that satisfy all criteria to be processed)."
                % (n_mpi_procs, len(valid_mic_id_substr_list)),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        break
    #
    # NOTE: Toshio Moriya 2016/10/24
    # The following function takes care of the case when an if-statement uses break for occurence of an error.
    # However, more elegant way is to use 'exception' statement of exception mechanism...
    #
    sp_utilities.if_error_then_all_processes_exit_program(error_status)

    # ====================================================================================
    # Obtain the list of micrograph id sustrings
    # ====================================================================================
    # --------------------------------------------------------------------------------
    # Prepare variables for this section
    # --------------------------------------------------------------------------------
    # Prepare variables related to options
    box_size = options.box_size
    box_half = old_div(box_size, 2)
    mask2d = sp_utilities.model_circle(
        old_div(box_size, 2), box_size, box_size
    )  # Create circular 2D mask to Util.infomask of particle images
    mic_resample_ratio = options.mic_resample_ratio

    # Micrograph baseroot pattern (extension are removed from micrograph basename pattern)
    # for substack file names
    mic_baseroot_pattern = os.path.splitext(mic_basename_pattern)[0]

    # Prepare the counters for the global summary of micrographs
    n_mic_process = 0
    n_mic_reject_no_rebox_entry = 0
    n_mic_reject_invalid_rebox_format = 0
    n_global_rebox_detect = 0
    n_global_rebox_process = 0
    n_global_rebox_reject_out_of_boundary = 0

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
        # Read the associated rebox according to the specified format and
        # make the rebox the center of particle image if necessary
        # Do this first because error might happen
        # --------------------------------------------------------------------------------
        rebox_path = global_entry_dict[mic_id_substr][subkey_rebox_path]
        assert os.path.exists(rebox_path)
        rebox_list = sp_utilities.read_text_row(rebox_path)

        if len(rebox_list) == 0:
            sp_global_def.sxprint(
                "For %s, the associate rebox file %s does not contain any entries. Skipping..."
                % (mic_basename, rebox_path)
            )
            n_mic_reject_no_rebox_entry += 1
            continue
        assert len(rebox_list) > 0
        if len(rebox_list[0]) != n_idx_rebox:
            sp_global_def.sxprint(
                "For %s, the format of the associate rebox file %s is not valid. It should contain %d columns. Skipping..."
                % (mic_basename, rebox_path, n_idx_rebox)
            )
            n_mic_reject_invalid_rebox_format += 1
            continue

        # --------------------------------------------------------------------------------
        # Get CTF parameter if necessary
        # Calculate the resampled pixel size and store it to the cter_entry if necessary
        # Do before expensive micrograph processing
        # --------------------------------------------------------------------------------
        if ctf_params_src is not None:
            cter_entry = global_entry_dict[mic_id_substr][subkey_cter_entry]

            # Get CTF parameters of this micrograph
            # indexes 0 to 7 (idx_cter_def to idx_cter_astig_ang) must be same in cter format & ctf object format.
            # rebox_ctf_obj = generate_ctf(cter_entry)
            #
            # NOTE: 2017/03/07 Toshio Moriya
            # Due to the change of error handling in generate_ctf()
            # the argument have to be a list with length of 6 or 8 now.
            #
            mic_ctf_entry = []
            mic_ctf_entry.append(cter_entry[idx_cter_def])
            mic_ctf_entry.append(cter_entry[idx_cter_cs])
            mic_ctf_entry.append(cter_entry[idx_cter_vol])
            mic_ctf_entry.append(cter_entry[idx_cter_apix])
            mic_ctf_entry.append(cter_entry[idx_cter_bfactor])
            mic_ctf_entry.append(cter_entry[idx_cter_total_ac])
            mic_ctf_entry.append(cter_entry[idx_cter_astig_amp])
            mic_ctf_entry.append(cter_entry[idx_cter_astig_ang])
            assert len(mic_ctf_entry) == n_idx_ctf

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
        # Apply the Gaussian high-pass Fourier filter to micrograph based on the (resampled) box size;
        # Cut off frequency components lower than one that the (resampled) box size can express.
        # Then, move back to the real space processing
        # --------------------------------------------------------------------------------
        mic_img = sp_fundamentals.fft(
            sp_filter.filt_gaussh(mic_img, old_div(mic_resample_ratio, box_size))
        )

        # --------------------------------------------------------------------------------
        # Resample micrograph, map rebox, and window segments from resampled micrograph using new rebox
        # after resampling by mic_resample_ratio, resampled pixel size = src_pixel_size/mic_resample_ratio
        # --------------------------------------------------------------------------------
        # NOTE: 2015/04/13 Toshio Moriya
        # resample() efficiently takes care of the case mic_resample_ratio = 1.0 but
        # it does not set apix_*. Even though it sets apix_* when mic_resample_ratio < 1.0...
        mic_img = sp_fundamentals.resample(mic_img, mic_resample_ratio)

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

        # --------------------------------------------------------------------------------
        # Prepare rebox loop variables
        # --------------------------------------------------------------------------------
        nx = mic_img.get_xsize()
        ny = mic_img.get_ysize()
        x0 = old_div(nx, 2)
        y0 = old_div(ny, 2)

        local_particle_id = 0  # can be different from coordinates_id
        rebox_reject_out_of_boundary_messages = []

        # Loop through rebox
        for rebox_id in range(len(rebox_list)):
            # Get particle coordinates
            x = int(rebox_list[rebox_id][idx_rebox_mic_coord_x])
            y = int(rebox_list[rebox_id][idx_rebox_mic_coord_y])
            rebox_resample_ratio = float(
                rebox_list[rebox_id][idx_rebox_mic_resample_ratio]
            )  # accumulated resample ratio including per-particle magnification adjustment runs, starting from the original micrograph resample ratio used in the original sxwindow run.

            ctf_entry = [0.0] * n_idx_ctf
            if ctf_params_src is None:
                # Get particle CTF parameters
                ctf_entry[idx_ctf_defocus] = float(
                    rebox_list[rebox_id][idx_rebox_ctf_defocus]
                )
                ctf_entry[idx_ctf_cs] = float(rebox_list[rebox_id][idx_rebox_ctf_cs])
                ctf_entry[idx_ctf_voltage] = float(
                    rebox_list[rebox_id][idx_rebox_ctf_voltage]
                )
                ctf_entry[idx_ctf_apix] = float(
                    rebox_list[rebox_id][idx_rebox_ctf_apix]
                )
                ctf_entry[idx_ctf_bfactor] = float(
                    rebox_list[rebox_id][idx_rebox_ctf_bfactor]
                )
                ctf_entry[idx_ctf_ampcont] = float(
                    rebox_list[rebox_id][idx_rebox_ctf_ampcont]
                )
                ctf_entry[idx_ctf_dfdiff] = float(
                    rebox_list[rebox_id][idx_rebox_ctf_dfdiff]
                )
                ctf_entry[idx_ctf_dfang] = float(
                    rebox_list[rebox_id][idx_rebox_ctf_dfang]
                )
            else:
                ctf_entry = mic_ctf_entry
            assert len(ctf_entry) == n_idx_ctf

            # Get particle projection parameters
            proj_entry = [0.0] * n_idx_proj
            proj_entry[idx_proj_phi] = float(rebox_list[rebox_id][idx_rebox_proj_phi])
            proj_entry[idx_proj_theta] = float(
                rebox_list[rebox_id][idx_rebox_proj_theta]
            )
            proj_entry[idx_proj_psi] = float(rebox_list[rebox_id][idx_rebox_proj_psi])
            proj_entry[idx_proj_sx] = float(rebox_list[rebox_id][idx_rebox_proj_sx])
            proj_entry[idx_proj_sy] = float(rebox_list[rebox_id][idx_rebox_proj_sy])
            assert len(proj_entry) == n_idx_proj

            # Get adjustment parameters
            pp_def_error_accum = float(
                rebox_list[rebox_id][idx_rebox_pp_def_error_accum]
            )
            pp_mag_error_accum = float(
                rebox_list[rebox_id][idx_rebox_pp_mag_error_accum]
            )

            # Get chunk id (Sorted group ID)
            chunk_id = float(rebox_list[rebox_id][idx_rebox_chunk_id])

            if ctf_params_src is not None:
                # In case of swapping CTF parameters, defocus error from rebox file does not mean anything...
                # Therefore, set defocus error to zero!!!
                pp_def_error_accum = 0.0

            # Create CTF object
            ctf_obj = sp_utilities.generate_ctf(ctf_entry)

            # Convert projection shifts to the original scale
            if rebox_resample_ratio < 1.0:
                assert rebox_resample_ratio > 0.0
                if my_mpi_proc_id == main_mpi_proc:
                    sp_global_def.sxprint("MRK_DEBUG: ")
                    sp_global_def.sxprint(
                        "MRK_DEBUG: BEFORE applying rebox_resample_ratio := {}".format(
                            rebox_resample_ratio
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: ctf_obj.apix := {}, proj_entry[idx_proj_sx] := {}, proj_entry[idx_proj_sy] := {}".format(
                            ctf_obj.apix,
                            proj_entry[idx_proj_sx],
                            proj_entry[idx_proj_sy],
                        )
                    )
                ctf_obj.apix *= rebox_resample_ratio
                proj_entry[idx_proj_sx] = old_div(
                    proj_entry[idx_proj_sx], rebox_resample_ratio
                )
                proj_entry[idx_proj_sy] = old_div(
                    proj_entry[idx_proj_sy], rebox_resample_ratio
                )
                if my_mpi_proc_id == main_mpi_proc:
                    sp_global_def.sxprint(
                        "MRK_DEBUG: AFTER applying rebox_resample_ratio := {}".format(
                            rebox_resample_ratio
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: ctf_obj.apix := {}, proj_entry[idx_proj_sx] := {}, proj_entry[idx_proj_sy] := {}".format(
                            ctf_obj.apix,
                            proj_entry[idx_proj_sx],
                            proj_entry[idx_proj_sy],
                        )
                    )
            else:
                assert rebox_resample_ratio == 1.0

            # Adjust CTF and projection shifts by rewindow resample ratio
            src_pixel_size = ctf_obj.apix
            if mic_resample_ratio < 1.0:
                assert mic_resample_ratio > 0.0
                if my_mpi_proc_id == main_mpi_proc:
                    sp_global_def.sxprint("MRK_DEBUG: ")
                    sp_global_def.sxprint(
                        "MRK_DEBUG: BEFORE applying mic_resample_ratio := {}".format(
                            mic_resample_ratio
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: ctf_obj.apix := {}, proj_entry[idx_proj_sx] := {}, proj_entry[idx_proj_sy] := {}".format(
                            ctf_obj.apix,
                            proj_entry[idx_proj_sx],
                            proj_entry[idx_proj_sy],
                        )
                    )
                # store the resampled pixel size to the cter_entry to generate CTF object of this micrograph
                ctf_obj.apix = old_div(src_pixel_size, mic_resample_ratio)
                proj_entry[idx_proj_sx] *= mic_resample_ratio
                proj_entry[idx_proj_sy] *= mic_resample_ratio
                if my_mpi_proc_id == main_mpi_proc:
                    sp_global_def.sxprint(
                        "MRK_DEBUG: AFTER applying mic_resample_ratio := {}".format(
                            mic_resample_ratio
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: ctf_obj.apix := {}, proj_entry[idx_proj_sx] := {}, proj_entry[idx_proj_sy] := {}".format(
                            ctf_obj.apix,
                            proj_entry[idx_proj_sx],
                            proj_entry[idx_proj_sy],
                        )
                    )
                if my_mpi_proc_id == main_mpi_proc:
                    sp_global_def.sxprint(
                        "Resample micrograph to pixel size %6.4f [A/Pixels] from %6.4f [A/Pixels] and window segments from resampled micrograph."
                        % (ctf_obj.apix, src_pixel_size)
                    )
            else:
                assert mic_resample_ratio == 1.0
                # Do nothing

            # Adjust ONLY projection shifts by accumulated per-particle magnification error
            # Since this error essentially brings back to all particles to the same pixel size
            # (i.e. micrograph pixel size or micrograph-level global pixel size),
            # the adjustment of ctf_obj.apix is not necessary after all.
            if pp_mag_error_accum != 1.0:
                assert pp_mag_error_accum > 0.0
                if my_mpi_proc_id == main_mpi_proc:
                    sp_global_def.sxprint("MRK_DEBUG: ")
                    sp_global_def.sxprint(
                        "MRK_DEBUG: BEFORE applying pp_mag_error_accum := {}".format(
                            pp_mag_error_accum
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: ctf_obj.apix := {}, proj_entry[idx_proj_sx] := {}, proj_entry[idx_proj_sy] := {}".format(
                            ctf_obj.apix,
                            proj_entry[idx_proj_sx],
                            proj_entry[idx_proj_sy],
                        )
                    )
                proj_entry[idx_proj_sx] = old_div(
                    proj_entry[idx_proj_sx], pp_mag_error_accum
                )
                proj_entry[idx_proj_sy] = old_div(
                    proj_entry[idx_proj_sy], pp_mag_error_accum
                )
                if my_mpi_proc_id == main_mpi_proc:
                    sp_global_def.sxprint(
                        "MRK_DEBUG: AFTER applying pp_mag_error_accum := {}".format(
                            pp_mag_error_accum
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: ctf_obj.apix := {}, proj_entry[idx_proj_sx] := {}, proj_entry[idx_proj_sy] := {}".format(
                            ctf_obj.apix,
                            proj_entry[idx_proj_sx],
                            proj_entry[idx_proj_sy],
                        )
                    )
            else:
                assert pp_mag_error_accum == 1.0
                # Do nothing

            # Rescale particle coordinates if necessary
            if mic_resample_ratio < 1.0:
                assert mic_resample_ratio > 0.0
                x = int(x * mic_resample_ratio)
                y = int(y * mic_resample_ratio)
            else:
                assert mic_resample_ratio == 1.0

            # Rescale box size if necessary
            resampled_box_size = box_size
            resampled_box_half = box_half
            if pp_mag_error_accum != 1.0:
                ### assert (pp_def_error_accum == 0.0) # For now, do not allow both parameters to be adjusted (Toshio 2018/03/25)
                resampled_box_size = int(
                    numpy.ceil(old_div(box_size, pp_mag_error_accum))
                )
                resampled_box_half = int(old_div(resampled_box_size, 2))
                # sxprint("MRK_DEBUG: pp_mag_error_accum := {}, resampled_box_size := {}, resampled_box_half := {}".format(pp_mag_error_accum, resampled_box_size, resampled_box_half))

            # Window a particle at this rebox
            if (
                (0 <= x - resampled_box_half)
                and (x + resampled_box_half < nx)
                and (0 <= y - resampled_box_half)
                and (y + resampled_box_half < ny)
            ):
                try:
                    particle_img = EMAN2_cppwrap.Util.window(
                        mic_img,
                        resampled_box_size,
                        resampled_box_size,
                        1,
                        x - x0,
                        y - y0,
                    )
                except:
                    sp_global_def.sxprint("MRK_DEBUG: ")
                    sp_global_def.sxprint(
                        "MRK_DEBUG: mic_resample_ratio := {}, box_size := {}, box_half := {}".format(
                            mic_resample_ratio, box_size, box_half
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: pp_mag_error_accum := {}, resampled_box_size := {}, resampled_box_half := {}".format(
                            pp_mag_error_accum, resampled_box_size, resampled_box_half
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: nx := {}, ny := {}".format(nx, ny)
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: x := {}, y := {}, x0 := {}, y0 := {}, x-x0 := {}, y-y0 := {}".format(
                            x, y, x0, y0, x - x0, y - y0
                        )
                    )
                    raise
            else:
                rebox_reject_out_of_boundary_messages.append(
                    "rebox ID = %04d: x = %4d, y = %4d, box_size = %4d "
                    % (rebox_id, x, y, box_size)
                )
                # sxprint("MRK_DEBUG: rebox_reject_out_of_boundary_messages[-1] := %s" % rebox_reject_out_of_boundary_messages[-1])
                continue

            if pp_mag_error_accum != 1.0:
                assert pp_mag_error_accum > 0.0
                particle_img = mrk_resample2d(
                    particle_img, pp_mag_error_accum, box_size
                )
                ### #
                ### # NOTE: Toshio Moriya 2018/05/17
                ### # Since the pp_mag_error_accum will essentially brings back to all particles to the same pixel size
                ### # (i.e. micrograph pixel size or micrograph-level global pixel size),
                ### # the adjustment of ctf_obj.apix is not necessary after all.
                ### #
                ### # Store the resampled pixel size to the cter_entry to generated (local) CTF object of this particle
                ### if my_mpi_proc_id == main_mpi_proc:
                ### 	sxprint("MRK_DEBUG: ")
                ### 	sxprint("MRK_DEBUG: BEFORE applying pp_mag_error_accum := {}".format(pp_mag_error_accum))
                ### 	sxprint("MRK_DEBUG: ctf_obj.apix := {}".format(ctf_obj.apix))
                ### ctf_obj.apix /= pp_mag_error_accum
                ### if my_mpi_proc_id == main_mpi_proc:
                ### 	sxprint("MRK_DEBUG: AFTER applying pp_mag_error_accum := {}".format(pp_mag_error_accum))
                ### 	sxprint("MRK_DEBUG: ctf_obj.apix := {}".format(ctf_obj.apix))
            # else:
            # 	# do nothing
            assert particle_img.get_xsize() == box_size
            assert particle_img.get_ysize() == box_size
            assert particle_img.get_zsize() == 1

            # Adjust defocus if necessary
            if pp_def_error_accum != 0.0:
                if my_mpi_proc_id == main_mpi_proc:
                    sp_global_def.sxprint("MRK_DEBUG: ")
                    sp_global_def.sxprint(
                        "MRK_DEBUG: BEFORE applying pp_def_error_accum := {}".format(
                            pp_def_error_accum
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: ctf_obj.defocus := {}".format(ctf_obj.defocus)
                    )
                ctf_obj.defocus += pp_def_error_accum
                if my_mpi_proc_id == main_mpi_proc:
                    sp_global_def.sxprint(
                        "MRK_DEBUG: AFTER applying pp_def_error_accum := {}".format(
                            pp_def_error_accum
                        )
                    )
                    sp_global_def.sxprint(
                        "MRK_DEBUG: ctf_obj.defocus := {}".format(ctf_obj.defocus)
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
                    "The standard deviation of the particle image associated with %dth rebox entry in micrograph %s for %s is zero unexpectedly. Skipping..."
                    % (rebox_id, mic_path, mic_basename)
                )
                continue

            # Set header entries (attributes) of this particle image
            #
            # NOTE: 2015/04/09 Toshio Moriya
            # ptcl_source_image might be redundant information...
            # Consider re-organizing header entries...
            #
            particle_img.set_attr("ptcl_source_image", mic_path)
            particle_img.set_attr(
                "ptcl_source_coord",
                [
                    int(rebox_list[rebox_id][idx_rebox_mic_coord_x]),
                    int(rebox_list[rebox_id][idx_rebox_mic_coord_y]),
                ],
            )
            particle_img.set_attr("ptcl_source_coord_id", rebox_id)
            particle_img.set_attr(
                "data_n", rebox_id
            )  # NOTE: Toshio Moriya 2017/11/20: same as ptcl_source_coord_id but the other program uses this header entry key...
            ### #
            ### # NOTE: Toshio Moriya 2018/05/17
            ### # Since the pp_mag_error_accum will essentially brings back to all particles to the same pixel size
            ### # (i.e. micrograph pixel size or micrograph-level global pixel size).
            ### # The accumulated_resample_ratio is not necessary after all.
            ### #
            ### particle_img.set_attr("resample_ratio", accumulated_resample_ratio) # NOTE: Toshio Moriya 2018/04/24: Store accumulated resample ratio instead of micrograph resample ratio used this sxrewindow run.
            particle_img.set_attr("resample_ratio", mic_resample_ratio)
            particle_img.set_attr("pp_def_error_accum", pp_def_error_accum)
            particle_img.set_attr("pp_mag_error_accum", pp_mag_error_accum)
            particle_img.set_attr("chunk_id", chunk_id)
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
            # apix_* attributes are updated by resample() only when mic_resample_ratio != 1.0
            # Let's make sure header info is consistent by setting apix_* = 1.0
            # regardless of options, so it is not passed down the processing line
            #
            particle_img.set_attr(
                "apix_x", 1.0
            )  # particle_img.set_attr("apix_x", resampled_pixel_size)
            particle_img.set_attr(
                "apix_y", 1.0
            )  # particle_img.set_attr("apix_y", resampled_pixel_size)
            particle_img.set_attr(
                "apix_z", 1.0
            )  # particle_img.set_attr("apix_z", resampled_pixel_size)
            particle_img.set_attr(
                "ptcl_source_apix", src_pixel_size
            )  # Store the original pixel size
            particle_img.set_attr("ctf", ctf_obj)
            particle_img.set_attr("ctf_applied", 0)

            sp_utilities.set_params_proj(particle_img, proj_entry)

            # Write the particle image to local stack file
            # sxprint("MRK_DEBUG: local_stack_path, local_particle_id", local_stack_path, local_particle_id)
            particle_img.write_image(local_stack_path, local_particle_id)
            local_particle_id += 1

        # Save the message list of rejected rebox because of out-of-boundary
        # sxprint("MRK_DEBUG: len(rebox_reject_out_of_boundary_messages) := %d" % len(rebox_reject_out_of_boundary_messages))
        if len(rebox_reject_out_of_boundary_messages) > 0:
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
                os.path.splitext(os.path.basename(rebox_path))[0]
                + "_reject_out_of_boundary.txt",
            )
            # sxprint("MRK_DEBUG: coords_reject_out_of_boundary_path := %s" % coords_reject_out_of_boundary_path)
            coords_reject_out_of_boundary_file = open(
                coords_reject_out_of_boundary_path, "w"
            )

            for (
                coords_reject_out_of_boundary_message
            ) in rebox_reject_out_of_boundary_messages:
                coords_reject_out_of_boundary_file.write(
                    coords_reject_out_of_boundary_message
                )
                coords_reject_out_of_boundary_file.write("\n")

            # Close the consistency check file, if necessary
            coords_reject_out_of_boundary_file.flush()
            coords_reject_out_of_boundary_file.close()

        # Update the counters for the global summary of micrographs
        n_mic_process += 1
        n_global_rebox_detect += len(rebox_list)
        n_global_rebox_process += local_particle_id
        n_global_rebox_reject_out_of_boundary += len(
            rebox_reject_out_of_boundary_messages
        )

        # Release the data base of local stack from this process
        # so that the subprocess can access to the data base
        EMAN2db.db_close_dict(local_stack_path)

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
        n_mic_reject_no_rebox_entry = mpi.mpi_reduce(
            n_mic_reject_no_rebox_entry,
            1,
            mpi.MPI_INT,
            mpi.MPI_SUM,
            main_mpi_proc,
            mpi.MPI_COMM_WORLD,
        )
        n_mic_reject_invalid_rebox_format = mpi.mpi_reduce(
            n_mic_reject_invalid_rebox_format,
            1,
            mpi.MPI_INT,
            mpi.MPI_SUM,
            main_mpi_proc,
            mpi.MPI_COMM_WORLD,
        )
        n_global_rebox_detect = mpi.mpi_reduce(
            n_global_rebox_detect,
            1,
            mpi.MPI_INT,
            mpi.MPI_SUM,
            main_mpi_proc,
            mpi.MPI_COMM_WORLD,
        )
        n_global_rebox_process = mpi.mpi_reduce(
            n_global_rebox_process,
            1,
            mpi.MPI_INT,
            mpi.MPI_SUM,
            main_mpi_proc,
            mpi.MPI_COMM_WORLD,
        )
        n_global_rebox_reject_out_of_boundary = mpi.mpi_reduce(
            n_global_rebox_reject_out_of_boundary,
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
            "Rejected by no rebox entries       : %6d" % (n_mic_reject_no_rebox_entry)
        )
        sp_global_def.sxprint(
            "Rejected by invalid rebox format   : %6d"
            % (n_mic_reject_invalid_rebox_format)
        )
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Global summary of rebox level processing...")
        sp_global_def.sxprint(
            "Detected                           : %6d" % (n_global_rebox_detect)
        )
        sp_global_def.sxprint(
            "Processed                          : %6d" % (n_global_rebox_process)
        )
        sp_global_def.sxprint(
            "Rejected by out of boundary        : %6d"
            % (n_global_rebox_reject_out_of_boundary)
        )
        if n_global_rebox_reject_out_of_boundary > 0:
            assert os.path.exists(reject_out_of_boundary_dir), "MRK_DEBUG"
            sp_global_def.sxprint(
                "    NOTE: Information of rejected rebox by out of boundary are saved in %s files."
                % (
                    os.path.join(
                        reject_out_of_boundary_dir,
                        os.path.splitext(os.path.basename(rebox_pattern))[0]
                        + "_reject_out_of_boundary.txt",
                    )
                )
            )
        else:
            assert n_global_rebox_reject_out_of_boundary == 0, "MRK_DEBUG"
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

def main():
    mpi.mpi_init(0, [])
    sp_global_def.print_timestamp("Start")
    run()
    sp_global_def.print_timestamp("Finish")
    mpi.mpi_finalize()


# ========================================================================================
# Define main function for command line execution
# ========================================================================================
if __name__ == "__main__":
    main()

# ========================================================================================
#  END OF FILE
# ========================================================================================
