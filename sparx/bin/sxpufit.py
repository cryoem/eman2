#!/usr/bin/python

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



from global_def import *
import global_def
from applications import pufit
from numpy import * 
from scipy import *
from scipy.optimize import fmin
from optparse import OptionParser
import sys
# Levenberg Marquardt least square optimization 
from scipy.optimize import leastsq # least square minimization
# Nelder-Mead simplex algorithm
from scipy.optimize import fmin # Nelder-Mead simplex algorithm
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " runmodel PDBTMP.spi roo.model --ps=PS (in Angstrom) --sgm1=half_width of 1/10 Angstrom peak --sgm2=half_width of 1/5 Angstrom peak "

	parser = OptionParser(usage,version=SPARXVERSION)
  
	parser.add_option("--ps", type="float", default=1.0, help="  pixel size in 1/angstrom ")
	
	parser.add_option("--sgm1", type="float",default=.05, help=" sigma of 1/10 Angstrom peak")
	
	parser.add_option("--sgm2", type="float",default=.2, help=" sigma of 1/5 Angstrom peak")
	
	(options, args) = parser.parse_args()
    
	if len(args) != 3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	runlist=[]
	runlist.insert(0,args[1])
	runlist.insert(1,args[0])
	os.system("echo CALL SPIDER TO GET THE DENSITY MAP")
	os.system(runlist[1])  
	pufit(args[1],args[2],options.sgm1,options.sgm2,options.ps)

if __name__ == "__main__":
	        main()



