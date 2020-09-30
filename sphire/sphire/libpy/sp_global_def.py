from __future__ import print_function
from __future__ import division
from past.utils import old_div

"""
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
Author: Pawel A.Penczek, 2006-09-09 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology

Copyright (c) 2000-2006 The University of Texas - Houston Medical School This 
software is issued under a joint BSD/GNU license. You may use the source code 
in this file under either license. However, note that the complete EMAN2 and 
SPARX software packages have some GPL dependencies, so you are responsible for 
compliance with the licenses of these packages if you opt to use BSD licensing.

The warranty disclaimer below holds in either instance.

This complete copyright notice must be included in any revised version of the
source code. Additional authorship citations may be added, but existing author
citations must be preserved.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details. You 
should have received a copy of the GNU General Public License along with this 
program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, 
Suite 330, Boston, MA  02111-1307 USA
"""

# --------------------------------------------------------------------[ header ]
import EMAN2_cppwrap
import EMAN2_meta
import inspect
import mpi
import random
import sys
import time
from . import sp_utilities as util
import os
import sphire

def get_timestamp(file_format=False):
    """
	Utility function to get a properly formatted timestamp.

	Args:
		file_format (bool): If true, timestamp will not include ':' characters
			for a more OS-friendly string that can be used in less risky file
			names [default: False ]
	"""
    if file_format:
        return time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    else:
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


# ------------------------------------------------------------[ util functions ]


def print_timestamp(tag=""):
    """
	Utility function to print a generic time stamp plus an optional tag.

	   NOTE: Using the tag "Start" will synchronize the SXPRINT_LOG path across
	all mpi processes (if mpi is being used). Throughout SPHIRE, this tag is
	being used before calling main() and after mpi_init() has been called.

	Args:
		tag (string): optional string that can be added to the time stamp to
			provide more information [default: ""]

	Example:
		print_timestamp( "Start" )
		[Start] : 2019-02-07 11:29:37
	"""

    # are we using mpi?
    try:
        mpi_rank = int(os.environ["OMPI_COMM_WORLD_RANK"])

        # printing the "Start"-tag will sync up the log file names so that all processes use the same logfile
        global SXPRINT_LOG, SXPRINT_LOG_SYNC
        if not SXPRINT_LOG_SYNC and "start" in tag.lower():
            SXPRINT_LOG = util.send_string_to_all(SXPRINT_LOG, 0)
            SXPRINT_LOG_SYNC = True

    # if there is no such thing as OMPI_COMM_WORLD_RANK, then we're not using mpi
    except KeyError:
        mpi_rank = 0

    if mpi_rank == 0:
        if tag != "":
            sxprint("[" + tag + "] : ", end="", print_timestamp=False)
        sxprint(get_timestamp(), print_timestamp=False)


def write_command(output_folder=None):
    try:
        mpi_rank = int(os.environ["OMPI_COMM_WORLD_RANK"])

    # if there is no such thing as OMPI_COMM_WORLD_RANK, then we're not using mpi
    except KeyError:
        mpi_rank = 0

    if mpi_rank == 0:
        command = " ".join(sys.argv) + "\n"
        global SXPRINT_OUTDIR_PATH, PRE_LOG
        if output_folder and PRE_LOG is not None:

            SXPRINT_OUTDIR_PATH = output_folder
            with open(os.path.join(SXPRINT_OUTDIR_PATH, SXPRINT_OUTDIR_FILENAME), "a+") as write:
                write.write('\n')
                write.write('\n')
                write.write('NEW COMMAND: {}'.format(command))
                write.write('\n')
                write.write('\n')
                write.write('\n'.join(PRE_LOG))
                write.write('\n')
            PRE_LOG = None

            with open(
                os.path.join(output_folder, "command.txt"), "a+"
            ) as the_command:
                the_command.write(command)

        global SXPRINT_LOG_EXISTS, SXPRINT_CMD_SKIP
        if not SXPRINT_CMD_SKIP:
            if not SXPRINT_LOG_EXISTS:
                try:
                    os.makedirs(SXPRINT_LOG_PATH)
                except OSError:
                    pass
                SXPRINT_LOG_EXISTS = True

            try:
                with open(SXPRINT_CMD, "a+") as f:
                    f.write(command)
            except IOError:
                print(
                    "Permission denied! Could not write to {0}.".format(
                        SXPRINT_LOG_PATH
                    )
                )
                SXPRINT_CMD_SKIP = True

        sxprint(command)


def sxprint(*args, **kwargs):
    """
	Generic print function that includes time stamps and caller id. Everything
	that is printed is also logged to file <SXPRINT_LOG> (can be disabled by
	setting SXPRINT_LOG to "").

	Args:
		*args: Variable number of arguments

		**kwargs: Dictionary containing (separate) variable number of (keyword)
			arguments

		filename (string): If a file name is provided the message is also
			written to file. If a file of the same name already exists the new
			message is appended to the end of the file. If no file of the given
			name exists it is created.

	Example:
		 sxprint( "This is " + "a %s" % "test" + ".", filename="out.log" )
		2019-02-07 13:36:50 <module> => This is a test.
	"""
    end = kwargs.get("end", "\n")
    print_timestamp = kwargs.get("print_timestamp", True)

    # prepend timestamp
    t = get_timestamp()
    f = sys._getframe(1).f_code.co_name
    if print_timestamp:
        m = t + " " + f + " => " + "  ".join(map(str, args))
    else:
        m = "  ".join(map(str, args))

    # print message to stdout
    print(m, end=end)
    sys.stdout.flush()

    if SXPRINT_OUTDIR_PATH is None:
        global PRE_LOG
        PRE_LOG.append(m)
    else:
        assert PRE_LOG is None
        with open(os.path.join(SXPRINT_OUTDIR_PATH, SXPRINT_OUTDIR_FILENAME), "a+") as write:
            write.write(m + end)

    # print message to SPHIRE execution log
    global SXPRINT_LOG_EXISTS, SXPRINT_LOG_SKIP
    if not SXPRINT_LOG_SKIP:
        if not SXPRINT_LOG_EXISTS:
            try:
                os.makedirs(SXPRINT_LOG_PATH)
            except OSError:
                pass
            SXPRINT_LOG_EXISTS = True

        try:
            with open(SXPRINT_LOG, "a+") as f:
                f.write(m + end)
        except IOError:
            print("Permission denied! Could not write to {0}.".format(SXPRINT_LOG_PATH))
            SXPRINT_LOG_SKIP = True

    # return printed message
    return m


def ERROR(message, where="", action=1, myid=0):
    """
	Utility function for consistent error throwing across sparx functions.

	Args:
		where (string): Location of error. Note that this will be determined automatically!
			(sxpipe.py is one exception)
		message (string): Error message
		action (0/1): Choose (1) error and abort, or (0) warning and continue [default: 1]
		myid (integer): mpi rank; used to only print error on main process (myid == 0)
	"""
    global BATCH
    global MPI

    if myid == 0:

        file = sys._getframe(1).f_code.co_filename  # NOTE: for this, inspect.stack can/
        func = sys._getframe(1).f_code.co_name  # should be used but the exact use
        line = sys._getframe(1).f_lineno  # differs from Python 2 to Python 3

        if action:
            sxprint(
                "ERROR reported by function '"
                + func
                + "' in file '"
                + file
                + "', line "
                + str(line)
                + ": "
            )
        else:
            sxprint(
                "WARNING reported by function '"
                + func
                + "' in file '"
                + file
                + "', line "
                + str(line)
                + ": "
            )

        sxprint(message)

    if action == 1 and BATCH:

        if MPI:
            mpi.mpi_finalize()
            MPI = False
            BATCH = False
            if myid == 0:
                sxprint("EXIT")
            sys.exit(1)

        else:
            BATCH = False
            if myid == 0:
                sxprint("EXIT")
            sys.exit(1)


# import


# set global random seed
rand_seed = EMAN2_cppwrap.Util.get_randnum_seed()
random.seed(rand_seed)

rand_seed = EMAN2_cppwrap.Util.get_randnum_seed()
EMAN2_cppwrap.Util.set_randnum_seed(rand_seed)

# ___________________________________________ User settings: change as necessary

# system-wide parameters
interpolation_method_2D = (
    "quadratic"
)  # set 2-D interpolation method ("linear", "quadratic", "gridding")
Eulerian_Angles = "SPIDER"  # set Euler angle convention ("SPIDER", "EMAN", "IMAGIC"):

"""Multiline Comment0"""
BATCH = False

"""Multiline Comment1"""
MPI = False

"""Multiline Comment2"""
CACHE_DISABLE = False

# ________________________________________ System settings: please do not change
#test_sphire_v_string = sphire.__init__.version

SPARXVERSION = "SPHIRE v"+sphire.__version__+" (GitHub: " + EMAN2_meta.DATESTAMP + ")"
SPARX_MPI_TAG_UNIVERSAL = 123456
SPARX_DOCUMENTATION_WEBSITE = "http://sparx-em.org/sparxwiki/"

# -------------------------------------------------------------------[ logging ]


# classic, user-exposed logfile (non-uniform use throughout sparx modules)
LOGFILE = "logfile_" + get_timestamp(file_format=True)
LOGFILE_HANDLE = 0
IS_LOGFILE_OPEN = False

# sxprint log (sxprint logging can be disabled by setting this to "")
SXPRINT_OUTDIR_PATH = None
SXPRINT_OUTDIR_FILENAME = 'logfile.txt'
PRE_LOG = []

SXPRINT_LOG_PATH = "SPHIRE_LOG_HISTORY"
SXPRINT_LOG_EXISTS = os.path.exists(SXPRINT_LOG_PATH)
SXPRINT_LOG_SKIP = False
SXPRINT_CMD_SKIP = False
try:
    init_func = inspect.re.search(
        "([^/]*).py", sys._getframe(len(inspect.stack()) - 1).f_code.co_filename
    ).group(1)
except AttributeError:
    init_func = "none"

SXPRINT_LOG = os.path.join(
    SXPRINT_LOG_PATH, get_timestamp(file_format=True) + "_" + init_func + ".log"
)
SXPRINT_LOG_SYNC = (
    False
)  # denotes whether SXPRINT_LOG has been synchronized across mpi processes
SXPRINT_CMD = os.path.join(SXPRINT_LOG_PATH, "commands.txt")
