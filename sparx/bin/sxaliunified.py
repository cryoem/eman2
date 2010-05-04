#!/usr/bin/env python
#
# Author: Chao Yang, 09/09/2006 (CYang@lbl.gov)
#         Lawrence Berkeley National Laboratory, Berkeley, CA
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


import os
import global_def
from global_def import *
from optparse import OptionParser
import sys
def main():
	if(MPI):
		from mpi import mpi_init
   		sys.argv = mpi_init( len(sys.argv), sys.argv )
        arglist = []
        for arg in sys.argv:
        	arglist.append( arg )
	progname = os.path.basename(arglist[0])
	usage = progname + " stack initial_vol --param=initial_angles_shifts_filename --maxit=max_iter --CTF --sym=c1 "
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--param",            type="string", default= "-1",               help="  initial angles and shifts contained in an ASCII file")	
	parser.add_option("--output",           type="string", default= "junk",             help="  the name of the output volume file")	
	parser.add_option("--ncpus",            type="string", default= "1",                help="  number of CPUS to be used")
	parser.add_option("--maxit",            type="string", default= 10,                 help="  maximum number of iterations allowed (set to 10 by default) ")
	parser.add_option("--CTF",              action="store_true", default=False,         help="  Consider CTF correction during the alignment (not implemented yet)")
	parser.add_option("--sym",              type="string", default= "c1",               help="  symmetry of the structure (not implemented yet)")
	(options, args) = parser.parse_args(arglist[1:])

	exefile = "/home/pawel/EMAN2/src/eman2/libEM/sparx/temp/rununified"

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

	if len(args) < 2 or len(args) >7 :
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
        elif (os.path.isfile(exefile)): 
        	xcommand = "mpirun -np %s %s -data=%s -model=%s -param=%s -out=%s -maxit=%s\n" % (options.ncpus,exefile,args[0],args[1],options.param,options.output,options.maxit)
		print "running " +  xcommand + " ... "
		global_def.BATCH = True
		os.system(xcommand)
		global_def.BATCH = False
	else:
		print "the binary executable 'rununified' is not available"
		print "cd into $SPXROOT/eman2/eman2/libEM/sparx/temp, and"
		print "type 'make rununified' first go generate this file"
if __name__ == "__main__":
	main()
