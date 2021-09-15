#!/usr/bin/env python

import EMAN2_cppwrap
import EMAN2db

import glob
import argparse
from sphire.libpy import sp_global_def
import os
import sys
import numpy as np
import time
import mpi

from sphire.libpy import sp_utilities
try:
	from pyStarDB import sp_pystardb as star
except:
	print("sp_star module requires a pyStarDB package")

def main():
    program_name = os.path.basename(sys.argv[0])
    parser = argparse.ArgumentParser(description="Similar to e2bdb.py for creating star file stacks")
    parser.add_argument(
        "--version", action="version", version=sp_global_def.SPARXVERSION
    )


    parser.add_argument(
        "root_dir",
        type=str,
        nargs='+',
        help="The root directory where all the data is present"
    )

    parser.add_argument(
        "--makevstack",
        type=str,
        help="Creates a 'virtual' star stack with its own metadata, "
             "but the data taken from the (filtered) list of stacks"
    )

    parser.add_argument(
        "--list",
        type=str,
        help="Specify the name of a file with a list of images to use in creation"
             " of virtual stacks. Please see source for details.",
    )

    parser.add_argument(
        "--exlist",
        type=str,
        help="Specify the name of a file with a list of images to not use in creation"
             " of virtual stacks. Please see source for details.",
    )

    parser.add_argument(
        "--step",
        type=str,
        default="0,1",
        help="Specify <init>,<step>[,<max>]. Processes only a subset of the input data."
             " For example, 0,2 would process only the even numbered particles")

    options = parser.parse_args()
    args = sys.argv[1:]


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

    if options.root_dir != None:
       root_dir = options.root_dir
    else:
        root_dir = os.getcwd()


    start = time.time()
    newout = []
    [
        newout.append(entry)  if entry.endswith('.star') and not os.path.isdir(entry)
        else newout.extend(sorted([entry for entry in glob.glob(os.path.join(entry, '*.star'))
                                 if os.path.isfile(entry)]))
        for entry in args if entry.endswith('.star') or os.path.isdir(entry) or '*' in entry
    ]

    if options.makevstack is None:
        options.makevstack = sys.argv[-1]

    if options.makevstack and options.list:
        indicies = list(np.loadtxt(options.list, dtype = int))
        vstack = star.StarFile.add_star(newout, inc_list=indicies )
        vstack.write_star_file(options.makevstack)

    elif options.makevstack and options.exlist:
        indicies = list(np.loadtxt(options.exlist, dtype = int))
        vstack = star.StarFile.add_star(newout, exc_list=indicies)
        vstack.write_star_file(options.makevstack)

    elif options.makevstack and options.step:
        try:
            newstep = list(map(int, options.step.split(',')))
        except ValueError:
            print("Step values need to be integers")
            sys.exit(1)

        if len(newstep) < 2 or len(newstep) > 3:
            print("Invalid --step specification, requires atleast two values")
            sys.exit(1)


        vstack = star.StarFile.add_star(newout, step=newstep)
        vstack.write_star_file(options.makevstack)

    else:
        vstack = star.StarFile.add_star(newout)
        vstack.write_star_file(options.makevstack)

    print("Time it took to create vstack", time.time() - start)






# ========================================================================================
# Define main function for command line execution
# ========================================================================================

if __name__ == "__main__":

    RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
    if RUNNING_UNDER_MPI:
        mpi.mpi_init(0, [])  # On OS X, there is an error if MPI is initialized and not finalized, hence the conditional
    sp_global_def.print_timestamp("Start")
    main()
    sp_global_def.print_timestamp("Finish")

    if RUNNING_UNDER_MPI:
        mpi.mpi_finalize()