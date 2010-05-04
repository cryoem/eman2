#! /usr/bin/env python
#
# Author: Zhong Huang, 09/09/2006 (zhong.huang@uth.tmc.edu)
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
import os
import global_def
from global_def import *
from optparse import OptionParser
import sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " mics power --the=theta,--bnd=band --yo=Y_overlap --xd=X_edge --yd=Y_edge --r=radius --prm=prefix of mics --w=256"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--the",  type="float",default=0.0, help=" tilt angle ")
	parser.add_option("--bnd",  type="float",default=-1, help=" number of bands in the micrographs along tilt axis ")
	parser.add_option("--yo", type="float",default=50.,  help="  Y overlap ratio 0-100 ")
	parser.add_option("--xd", type="float",default=500., help="  X edge 0,100,500,...")
	parser.add_option("--yd", type="float",default=500., help="  Y edge 0,100,500,...")
	parser.add_option("--r",  type="float",default=10.,  help="  mask radius for display, set as 10 " )
	parser.add_option("--prm",type="string",default="mic_", help="  mask radius for display, set as 10 " )
	parser.add_option("--w",  type="float",default=256,  help=" image size " )	
	(options, args) = parser.parse_args()
	if len(args) != 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	else:
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
		from applications  import  pw2sp_tilt
		global_def.BATCH = True
		pw2sp_tilt(args[0],args[1],options.the,options.bnd,options.yo,options.xd,options.yd,options.r,options.prm,options.w)
		global_def.BATCH = False

if __name__ == "__main__":
	        main()
