#! /usr/bin/env python
import global_def
from global_def import *

from random import seed
from string import atoi,replace,atof
from mpi import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
from projection import cmlines_voronoi_optall
import sys

if global_def.CACHE_DISABLE:
	from utilities import disable_bdb_cache
	disable_bdb_cache()

if MPI:
	sys.argv = mpi_init( len(sys.argv), sys.argv )
	
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	myid  = mpi_comm_rank(MPI_COMM_WORLD)

	if len(sys.argv) != 6 :
    		print "Usage: sxcmlines_optall.py stack seed outdir prog delta"
		sys.exit(-1)

	prj_stack = sys.argv[1]

	all_seed  = atoi( sys.argv[2] ) + myid
	out_dir   = "%s%4d" % (sys.argv[3],myid)
	out_dir   = replace( out_dir, ' ', '0' )

	progress  = "%s%4d" % (sys.argv[4],myid)
	progress  = replace( progress, ' ', '0' )

	delta = atof( sys.argv[5] )

	seed(all_seed)

	fproj = open( progress, 'w' )

	from sys import stdout
	cmlines_voronoi_optall(prj_stack, delta, out_dir, fproj )
else:
	argv = sys.argv
	if len(argv) !=6 :
		print "Usage: sxcmlines_optall.py stack seed outdir prog delta"
		sys.exit(-1)
	prj_stack = sys.argv[1]
	all_seed  = atoi(sys.argv[2])
	out_dir   = sys.argv[3]
	progress  = sys.argv[4]
	delta = atof(sys.argv[5])
	seed(all_seed)
	fproj = open(progress,'w')
	cmlines_voronoi_optall(prj_stack, delta, out_dir, fproj)
	
	


