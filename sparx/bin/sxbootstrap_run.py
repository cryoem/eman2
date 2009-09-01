#!/usr/bin/env python

import global_def
from global_def import *
from optparse import OptionParser
import os
import sys
import string

arglist = []
for arg in sys.argv:
	arglist.append( arg )

progname = os.path.basename( arglist[0] )
usage = progname + " prj_stack buf_prefix outdir nvol --CTF --snr=signal_noise_ratio --sym=symmetry -verbose=(0|1) --MPI"
parser = OptionParser(usage, version=SPARXVERSION)

parser.add_option("--CTF", action="store_true", default=False, help="  Consider CTF correction during the alignment" )
parser.add_option("--snr", type="float", default=1.0, help="Signal-to-Noise Ratio of input projection" )
parser.add_option("--sym", type="string", default="c1", help="Symmetry of input projection" )
parser.add_option("--verbose", type="int", default=0, help="verbosity" )
parser.add_option("--MPI", action="store_true", default=False,     help="  whether using MPI version ")

(options,args)=parser.parse_args(arglist[1:])

if len(args) != 4 :
	print usage
	sys.exit(-1)

prj_stack = args[0]
buf_prefix= args[1]
outdir = args[2]
nvol = string.atoi(args[3])

from applications import bootstrap_run

if global_def.CACHE_DISABLE:
	from utilities import disable_bdb_cache
	disable_bdb_cache()

if options.MPI:
	from mpi import mpi_init
	sys.argv = mpi_init(len(sys.argv), sys.argv)

global_def.BATCH = True

bootstrap_run(prj_stack, buf_prefix, outdir, nvol, options.CTF, options.snr, options.sym, options.verbose, options.MPI)
global_def.BATCH = False

