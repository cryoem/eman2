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
from global_def import *
from optparse import OptionParser
import sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " roodir --mtd=s --oup=w --ps=1. --volt=120 --cs=2 --wgh=.1 --rndf=100 --dz_max=50000. --f_l=30 --f_h=8 --nr1=5 --nr2=5 --prf=roo --f=a --skip=# --micdir=mics --pnt=n"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--mtd",    type="string", default="f",   help=" options for algorithm, mtd=s will use slow while accurate algorithm ")
	parser.add_option("--oup",    type="string", default="w",   help=" 1. defocus output to text file (2) and return as list(3)")
	parser.add_option("--ps",     type="float",  default=1,     help=" pixel size of the 1D-power-spectrum-calculated micrograph ")
	parser.add_option("--volt",   type="float",  default=120,   help=" accelerating voltage ")
	parser.add_option("--cs",     type="float",  default=2,     help=" spherical abberation constant of the microscope")
	parser.add_option("--wgh",    type="float",  default=.1,    help=" amplitude constrast ratio")
	parser.add_option("--rndf",   type="float",  default=100,   help=" round off value of the estimated defocus ")
	parser.add_option("--dz_max", type="float",  default=50000, help=" The the maximum value of the estimated defocus ")
	parser.add_option("--f_l",    type="float",  default=30,    help=" low frequency cut_off  ")
	parser.add_option("--f_h",    type="float",  default=8,     help=" high frequency cut_off ")		
	parser.add_option("--nr1",    type="float",  default=5,     help=" highest polynomial rank for envelope fitting  ")
	parser.add_option("--nr2",    type="float",  default=5,     help=" highest polynomial rank for baseline noise fitting  ")
	parser.add_option("--prf",    type="string", default="roo", help=" prefix of 1D power spectrum text file  ")
	parser.add_option("--f",      type="string", default="s",   help=" format of 1D power spectrum text file  ")
	parser.add_option("--skip",   type="string", default=";",   help=" commented lines in 1D power spectrum text file ")
	parser.add_option("--micdir", type="string", default="no",  help=" write CTF parameters back into micrograph header")
	parser.add_option("--pnt",    type="string", default="no",  help=" print defocus guessing process if the flag is p ")				
	(options, args) = parser.parse_args()
    	if len(args) != 1:
        	print "usage: " + usage
        	print "Please run '" + progname + " -h' for detailed options"
	else: 	
		from applications import defocus_calc

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		global_def.BATCH = True
		defocus_calc(args[0],options.mtd, options.oup, options.ps, options.volt, options.cs, options.wgh, options.rndf, options.dz_max, options.f_l, options.f_h, options.nr1, options.nr2, options.prf, options.f, options.skip, options.micdir, options.pnt)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
