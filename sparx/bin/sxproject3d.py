#! /usr/bin/env python

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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#


import os
import global_def
from   global_def import *
from   optparse import OptionParser
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " volume stack  <maskfile> --delta=angular_step --method=S --phiEqpsi=Minus --symmetry=c1"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--delta",    type="float", default=2, help="  angular step ")
	parser.add_option("--phiEqpsi", type="string",  default="Minus", help="  if Minus, psi is set to minus phi (default), if Zero, set to zero ")
	parser.add_option("--method",   type="string", default="S", help="  method of quasi-uniformly distributing Eulerian angles S (default) or P")
	parser.add_option("--symmetry", type="string", default="c1", help="  symmetry group")
	parser.add_option("--angles",   type="string", default=None, help="  List of angles (phi, theta, psi)")
	
	# extensions to generate noise, ctf and to use MPI
	parser.add_option("--noise",type="float",default=None,help="add Gaussian noise with relative SNR of N")
	parser.add_option("--CTF",type="string",default=None,help="list of defocus values")
	parser.add_option("--MPI",action="store_true",default=False,help="use MPI")

	(options, args) = parser.parse_args()
	if(len(args) < 2 or len(args) > 3):
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if len(args) == 2:
			mask = None
		else:
			mask = args[2]
			
		from   applications import project3d
		global_def.BATCH = True
		project3d(args[0], args[1], mask, options.delta, options.method, options.phiEqpsi, options.symmetry, options.angles, listctfs=options.CTF,noise=options.noise,MPI=options.MPI)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
