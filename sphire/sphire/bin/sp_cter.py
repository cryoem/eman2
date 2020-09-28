#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Pawel A.Penczek 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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
import inspect
import mpi
import optparse
import os
from ..libpy import sp_global_def
from ..libpy import sp_morphology
from ..libpy import sp_utilities
import sys


sp_global_def.BATCH = True


def run():
    program_name = os.path.basename(sys.argv[0])
    usage = (
        program_name
        + """  input_image_path  output_directory  --selection_list=selection_list  --wn=CTF_WINDOW_SIZE --apix=PIXEL_SIZE  --Cs=CS  --voltage=VOLTAGE  --ac=AMP_CONTRAST  --f_start=FREA_START  --f_stop=FREQ_STOP  --vpp  --kboot=KBOOT  --overlap_x=OVERLAP_X  --overlap_y=OVERLAP_Y  --edge_x=EDGE_X  --edge_y=EDGE_Y  --check_consistency  --stack_mode  --debug_mode

Automated estimation of CTF parameters with error assessment.

All Micrographs Mode - Process all micrographs in a directory: 
	Specify a list of input micrographs using a wild card (*), called here input micrographs path pattern. 
	Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). 
	Running from the command line requires enclosing the string by single quotes (') or double quotes ("). 
	sxgui.py will automatically adds single quotes to the string. 
	BDB files can not be selected as input micrographs. 
	Then, specify output directory where all outputs should be saved. 
	In this mode, all micrographs matching the path pattern will be processed.

	mpirun -np 16 sxcter.py './mic*.hdf' outdir_cter --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0

Selected Micrographs Mode - Process all micrographs in a selection list file:
	In addition to input micrographs path pattern and output directry arguments, 
	specify a name of micrograph selection list text file using --selection_list option 
	(e.g. output of sxgui_unblur.py or sxgui_cter.py). The file extension must be ".txt". 
	In this mode, only micrographs in the selection list which matches the file name part of the pattern (ignoring the directory paths) will be processed. 
	If a micrograph name in the selection list does not exists in the directory specified by the micrograph path pattern, processing of the micrograph will be skipped.

	mpirun -np 16 sxcter.py './mic*.hdf' outdir_cter --selection_list=mic_list.txt --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0

Single Micrograph Mode - Process a single micrograph: 
	In addition to input micrographs path pattern and output directry arguments, 
	specify a single micrograph name using --selection_list option. 
	In this mode, only the specified single micrograph will be processed. 
	If this micrograph name does not matches the file name part of the pattern (ignoring the directory paths), the process will exit without processing it. 
	If this micrograph name matches the file name part of the pattern but does not exists in the directory which specified by the micrograph path pattern, again the process will exit without processing it. 
	Use single processor for this mode.

	sxcter.py './mic*.hdf' outdir_cter --selection_list=mic0.hdf --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0

Stack Mode - Process a particle stack (Not supported by SPHIRE GUI)):: 
	Use --stack_mode option, then specify the path of particle stack file (without wild card "*") and output directory as arguments. 
	This mode ignores --selection_list, --wn --overlap_x, --overlap_y, --edge_x, and --edge_y options. 
	Use single processor for this mode. Not supported by SPHIRE GUI (sxgui.py). 

	sxcter.py bdb:stack outdir_cter --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0 --stack_mode

"""
    )
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)
    parser.add_option(
        "--selection_list",
        type="string",
        default=None,
        help="Micrograph selecting list: Specify path of a micrograph selection list text file for Selected Micrographs Mode. The file extension must be '.txt'. Alternatively, the file name of a single micrograph can be specified for Single Micrograph Mode. (default none)",
    )
    parser.add_option(
        "--wn",
        type="int",
        default=512,
        help="CTF window size [pixels]: The size should be slightly larger than particle box size. This will be ignored in Stack Mode. (default 512)",
    )
    parser.add_option(
        "--apix",
        type="float",
        default=-1.0,
        help="Pixel size [A/Pixels]: The pixel size of input micrograph(s) or images in input particle stack. (default -1.0)",
    )
    parser.add_option(
        "--Cs",
        type="float",
        default=2.0,
        help="Microscope spherical aberration (Cs) [mm]: The spherical aberration (Cs) of microscope used for imaging. (default 2.0)",
    )
    parser.add_option(
        "--voltage",
        type="float",
        default=300.0,
        help="Microscope voltage [kV]: The acceleration voltage of microscope used for imaging. (default 300.0)",
    )
    parser.add_option(
        "--ac",
        type="float",
        default=10.0,
        help="Amplitude contrast [%]: The typical amplitude contrast is in the range of 7% - 14%. The value mainly depends on the thickness of the ice embedding the particles. (default 10.0)",
    )
    parser.add_option(
        "--f_start",
        type="float",
        default=-1.0,
        help="Lowest resolution [A]: Lowest resolution to be considered in the CTF estimation. Determined automatically by default. (default -1.0)",
    )
    parser.add_option(
        "--f_stop",
        type="float",
        default=-1.0,
        help="Highest resolution [A]: Highest resolution to be considered in the CTF estimation. Determined automatically by default. (default -1.0)",
    )
    parser.add_option(
        "--kboot",
        type="int",
        default=16,
        help="Number of CTF estimates per micrograph: Used for error assessment. (default 16)",
    )
    parser.add_option(
        "--overlap_x",
        type="int",
        default=50,
        help="X overlap [%]: Overlap between the windows in the x direction. This will be ignored in Stack Mode. (default 50)",
    )
    parser.add_option(
        "--overlap_y",
        type="int",
        default=50,
        help="Y overlap [%]: Overlap between the windows in the y direction. This will be ignored in Stack Mode. (default 50)",
    )
    parser.add_option(
        "--edge_x",
        type="int",
        default=0,
        help="Edge x [pixels]: Defines the edge of the tiling area in the x direction. Normally it does not need to be modified. This will be ignored in Stack Mode. (default 0)",
    )
    parser.add_option(
        "--edge_y",
        type="int",
        default=0,
        help="Edge y [pixels]: Defines the edge of the tiling area in the y direction. Normally it does not need to be modified. This will be ignored in Stack Mode. (default 0)",
    )
    parser.add_option(
        "--check_consistency",
        action="store_true",
        default=False,
        help="Check consistency of inputs: Create a text file containing the list of inconsistent Micrograph ID entries (i.e. inconsist_mic_list_file.txt). (default False)",
    )
    parser.add_option(
        "--stack_mode",
        action="store_true",
        default=False,
        help="Use stack mode: Use a stack as the input. Please set the file path of a stack as the first argument and output directory for the second argument. This is advanced option. Not supported by sxgui. (default False)",
    )
    parser.add_option(
        "--debug_mode",
        action="store_true",
        default=False,
        help="Enable debug mode: Print out debug information. (default False)",
    )
    parser.add_option(
        "--vpp",
        action="store_true",
        default=False,
        help="Volta Phase Plate - fit smplitude contrast. (default False)",
    )
    parser.add_option(
        "--defocus_min",
        type="float",
        default=0.3,
        help="Minimum defocus search [um] (default 0.3)",
    )
    parser.add_option(
        "--defocus_max",
        type="float",
        default=9.0,
        help="Maximum defocus search [um] (default 9.0)",
    )
    parser.add_option(
        "--defocus_step",
        type="float",
        default=0.1,
        help="Step defocus search [um] (default 0.1)",
    )
    parser.add_option(
        "--phase_min",
        type="float",
        default=5.0,
        help="Minimum phase search [degrees] (default 5.0)",
    )
    parser.add_option(
        "--phase_max",
        type="float",
        default=175.0,
        help="Maximum phase search [degrees] (default 175.0)",
    )
    parser.add_option(
        "--phase_step",
        type="float",
        default=5.0,
        help="Step phase search [degrees] (default 5.0)",
    )
    parser.add_option(
        "--pap",
        action="store_true",
        default=False,
        help="Use power spectrum for fitting. (default False)",
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
        ####mpi.mpi_init( 0, [] )
        my_mpi_proc_id = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        n_mpi_procs = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
        sp_global_def.MPI = True

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
    # Check error conditions of arguments and options, then prepare variables for arguments
    # ------------------------------------------------------------------------------------
    input_image_path = None
    output_directory = None
    # not a real while, an if with the opportunity to use break when errors need to be reported
    error_status = None
    # change input unit
    freq_start = -1.0
    freq_stop = -1.0

    if options.f_start > 0.0:
        if options.f_start <= 0.5:
            sp_global_def.ERROR(
                "f_start should be in Angstrom"
            )  # exclude abs frequencies and spatial frequencies
        else:
            freq_start = old_div(1.0, options.f_start)

    if options.f_stop > 0.0:
        if options.f_stop <= 0.5:
            sp_global_def.ERROR(
                "f_stop should be in Angstrom"
            )  # exclude abs frequencies and spatial frequencies
        else:
            freq_stop = old_div(1.0, options.f_stop)

    while True:
        # --------------------------------------------------------------------------------
        # Check the number of arguments. If OK, then prepare variables for them
        # --------------------------------------------------------------------------------
        if len(args) != 2:
            error_status = (
                "Please check usage for number of arguments.\n Usage: "
                + usage
                + "\n"
                + "Please run %s -h for help." % (program_name),
                inspect.getframeinfo(inspect.currentframe()),
            )
            break

        # NOTE: 2015/11/27 Toshio Moriya
        # Require single quotes (') or double quotes (") when input micrograph pattern is give for input_image_path
        #  so that sys.argv does not automatically expand wild card and create a list of file names
        #
        input_image_path = args[0]
        output_directory = args[1]

        # --------------------------------------------------------------------------------
        # NOTE: 2016/03/17 Toshio Moriya
        # cter_mrk() will take care of all the error conditions
        # --------------------------------------------------------------------------------

        break
    sp_utilities.if_error_then_all_processes_exit_program(error_status)
    #  Toshio, please see how to make it informative
    assert input_image_path != None, " directory  missing  input_image_path"
    assert output_directory != None, " directory  missing  output_directory"

    if options.vpp == False:
        wrong_params = False
        vpp_options = [
            "--defocus_min",
            "--defocus_max",
            "--defocus_step",
            "--phase_min",
            "--phase_max",
            "--phase_step",
        ]
        for command_token in sys.argv:
            for vppo in vpp_options:
                if str.find(command_token, vppo) > -1:
                    wrong_params = True
                if wrong_params:
                    break
            if wrong_params:
                break
        if wrong_params:
            sp_global_def.ERROR(
                "Some options are valid only for Volta Phase Plate command  s"
                % command_token,
                myid=my_mpi_proc_id,
            )

    if my_mpi_proc_id == main_mpi_proc:
        command_line = ""
        for command_token in sys.argv:
            command_line += command_token + "  "
        sp_global_def.sxprint(" ")
        sp_global_def.sxprint("Shell line command:")
        sp_global_def.sxprint(command_line)

    if options.vpp:
        vpp_options = [
            options.defocus_min,
            options.defocus_max,
            options.defocus_step,
            options.phase_min,
            options.phase_max,
            options.phase_step,
        ]
        result = sp_morphology.cter_vpp(
            input_image_path,
            output_directory,
            options.selection_list,
            options.wn,
            options.apix,
            options.Cs,
            options.voltage,
            options.ac,
            freq_start,
            freq_stop,
            options.kboot,
            options.overlap_x,
            options.overlap_y,
            options.edge_x,
            options.edge_y,
            options.check_consistency,
            options.stack_mode,
            options.debug_mode,
            program_name,
            vpp_options,
            RUNNING_UNDER_MPI,
            main_mpi_proc,
            my_mpi_proc_id,
            n_mpi_procs,
        )
    elif options.pap:
        result = sp_morphology.cter_pap(
            input_image_path,
            output_directory,
            options.selection_list,
            options.wn,
            options.apix,
            options.Cs,
            options.voltage,
            options.ac,
            freq_start,
            freq_stop,
            options.kboot,
            options.overlap_x,
            options.overlap_y,
            options.edge_x,
            options.edge_y,
            options.check_consistency,
            options.stack_mode,
            options.debug_mode,
            program_name,
            RUNNING_UNDER_MPI,
            main_mpi_proc,
            my_mpi_proc_id,
            n_mpi_procs,
        )
    else:
        result = sp_morphology.cter_mrk(
            input_image_path,
            output_directory,
            options.selection_list,
            options.wn,
            options.apix,
            options.Cs,
            options.voltage,
            options.ac,
            freq_start,
            freq_stop,
            options.kboot,
            options.overlap_x,
            options.overlap_y,
            options.edge_x,
            options.edge_y,
            options.check_consistency,
            options.stack_mode,
            options.debug_mode,
            program_name,
            RUNNING_UNDER_MPI,
            main_mpi_proc,
            my_mpi_proc_id,
            n_mpi_procs,
        )

    if RUNNING_UNDER_MPI:
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if main_mpi_proc == my_mpi_proc_id:
        if options.debug_mode:
            sp_global_def.sxprint("Returned value from cter_mrk() := ", result)
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
    RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
    if RUNNING_UNDER_MPI:
        mpi.mpi_init(
            0, []
        )  # On OS X, there is an error if MPI is initialized and not finalized, hence the conditional
    sp_global_def.print_timestamp("Start")
    run()
    sp_global_def.print_timestamp("Finish")
    if RUNNING_UNDER_MPI:
        mpi.mpi_finalize()

if __name__ == "__main__":
    main()
