#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
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

''' variables governing system performance - can be changed by the user'''
# 2-D interpolation method:
#    "linear", "quadratic", "gridding"

interpolation_method_2D = "quadratic"

# Eulerian angles:
#    SPIDER, EMAN, IMAGIC

Eulerian_Angles = "SPIDER"

# NOTICE: beginning from version 0.70, we will no longer use MPI as a global variable
# Instead, the user would add mpi as a parameter for command line, example
# mpirun -np 10 sxali2d_c.py  ...  --mpi

# We read the global seed here. If the user wish to repeat the random results twice,
# he/she should first set the rand_seed to a fixed number and then run the program twice.
from   EMAN2   import Util, EMData, EMUtil, Transform
from   e2version import CVSDATESTAMP
from   random  import seed

rand_seed = Util.get_randnum_seed()
seed(rand_seed)

rand_seed = Util.get_randnum_seed()
Util.set_randnum_seed(rand_seed)

# BATCH flag should generally be set to False, which indicates that the output should be both displayed on the screen and written to the log file.
# However, the user may change it to True (either here or in other programs) so that the output is only written to the log file.
BATCH = False

# variable for disabling bdb cache use, For running sparx on clusters, set it to True to disable cache,
CACHE_DISABLE = False


global LOGFILE 
LOGFILE = "logfile"
from time import localtime, strftime
# timestring = strftime("_%d_%b_%Y_%H_%M_%S", localtime())
timestring = strftime("_%Y_%m_%d_%H_%M_%S", localtime())
LOGFILE = LOGFILE+timestring
LOGFILE_HANDLE = 0
IS_LOGFILE_OPEN = False
'''   SYSTEM FUNCTIONS - please do not change the text below '''
global SPARXVERSION
SPARXVERSION = "SPARX v3.0" + ' (CVS' + CVSDATESTAMP[6:-2] +')'

global SPARX_DOCUMENTATION_WEBSITE
SPARX_DOCUMENTATION_WEBSITE = "http://sparx-em.org/sparxwiki/"


def ERROR(message, where, action = 1, myid = 0):
	"""
		General error function for sparx system
		where:   function name
		message: error message
		action: 1 - fatal error, exit; 0 - non-fatal, print a warning
	"""
	if myid == 0:
		if action: print  "\n  *****  ERROR in: %s"%(where)
		else:      print  "\n  *****  WARNING in: %s"%(where)
		print "  *****  %s"%message
		print ""
	if action and BATCH:
		from sys import exit
		exit()
