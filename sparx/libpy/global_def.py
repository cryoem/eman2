
from __future__ import print_function
"""
Author: Pawel A.Penczek, 2006-09-09 (Pawel.A.Penczek@uth.tmc.edu)
Author: Fabian Schoenfeld, 2019-02-07 (fabian.schoenfeld@mpi-dortmund.mpg.de)

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

#--------------------------------------------------------------------[ header ]

# import
import sys
import time

import mpi

from EMAN2  import Util, EMData, EMUtil, Transform
from random import seed

# system-wide parameters
interpolation_method_2D = "quadratic"  # set 2-D interpolation method ("linear", "quadratic", "gridding")
Eulerian_Angles         = "SPIDER"     # set Euler angle convention ("SPIDER", "EMAN", "IMAGIC"):

# set global random seed
rand_seed = Util.get_randnum_seed()
seed(rand_seed)

rand_seed = Util.get_randnum_seed()
Util.set_randnum_seed(rand_seed)

#___________________________________________ User settings: change as necessary

"""
BATCH flag should generally be set to False, which indicates that the output
should be both displayed on the screen and written to the log file.
However, the user may change it to True (either here or in other programs) so
that the output is only written to the log file.
"""
BATCH = False

"""
NOTICE: beginning from version 0.70, we will no longer use MPI as a global 
variable. Instead, the user would add mpi as a parameter for command line, 
example: $ mpirun -np 10 sxali2d_c.py  ...  --mpi
NOTE: Version 0.70 is pretty old by now; it's unclear, however, whether this 
can be removed safely [2019-0207,Fabian]
"""
MPI = False

"""
variable for disabling bdb cache use, For running sparx on clusters, set it to
True to disable cache
"""
CACHE_DISABLE = False

# global logfile setup
global LOGFILE
LOGFILE = "logfile" + time.strftime("_%Y-%m-%d_%H-%M-%S", time.localtime())
LOGFILE_HANDLE  = 0
IS_LOGFILE_OPEN = False

#________________________________________ System settings: please do not change

from EMAN2_meta import DATESTAMP

global SPARXVERSION
SPARXVERSION = "SPARX v4.0" + ' (GITHUB: ' + DATESTAMP +')'

global SPARX_MPI_TAG_UNIVERSAL
SPARX_MPI_TAG_UNIVERSAL = 123456

global SPARX_DOCUMENTATION_WEBSITE
SPARX_DOCUMENTATION_WEBSITE = "http://sparx-em.org/sparxwiki/"

#------------------------------------------------------------[ util functions ]

def print_timestamp( tag="" ):
	"""
	Utility function to print a generic time stamp.

	Args:
		tag (string): optional string that can be added to the time stamp to
			provide more information [default: ""]

	Example:
		>>>  print_timestamp( "Start" )
		[Start] : 2019-02-07 11:29:37
	"""

	if tag != "": 
		print( "["+tag+"] : ", end="" )

	print( time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) )

def ERROR( message, where, action=1, myid=0 ):
	"""
	Utility function for consistent error throwing across sparx functions.

	Args:
		where (string): Location of error (e.g. "sxsummovie.main")
		message (string): Error message
		action (0/1): Choose (1) error and abort, or (0) warning and continue [default: 1]
		myid (integer): mpi rank; used to only print error on main process (myid == 0)
	"""
	global BATCH
	global MPI
	
	if myid == 0:

		if action: 
			print( "\n  *****  ERROR in: %s" % where )
		else:      
			print( "\n  *****  WARNING in: %s" % where )
		print("  *****  %s\n" % message)

	if action == 1 and BATCH:
		
		if  MPI:
			mpi.mpi_finalize()
			MPI   = False
			BATCH = False
			print_timestamp( "ABORT" )
			sys.exit(1)

		else:
			BATCH = False
			print_timestamp( "ABORT" )
			sys.exit(1)
