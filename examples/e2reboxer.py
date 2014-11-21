#!/usr/bin/env python

# Author: Steven Ludtke, 06/16/14 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
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
import sys
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]
This program will re-extract boxed particles from micrographs based on the box locations in the info/*.json files.
This is normally used to change the particle box-size.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--boxsize", type=int, help="Box size for particle extraction, may be larger than the standard size for the project to compensate for motion",default=-1)
	parser.add_argument("--invert",action="store_true",help="Extracted images are multiplied by -1 ",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if options.boxsize<10 : 
		print "Please specify a box size"
		sys.exit(1)

	box=options.boxsize

	try: os.mkdir("particles")
	except: pass

	logid=E2init(sys.argv)

	micros=[base_name(i) for i in os.listdir("info") if "_info.json" in i]
	
	for m in micros:
		if not os.path.exists("micrographs/{}.hdf".format(m)):
			if options.verbose>0 : print "No micrograph for: ",m
			continue
		
		db=js_open_dict(info_name(m))
		if not db.has_key("boxes") :
			db.close()
			if options.verbose>0 : print "No boxes in: ",m
			continue
		
		microfsp="micrographs/{}.hdf".format(m)
		micro=EMData(microfsp,0)
		if options.verbose>0 : print "{} has {} boxes".format(microfsp,len(db["boxes"]))
		
		for b in db["boxes"]:
			ptcl=fm.get_clip(Region(b[0]-box/2,b[1]-box/2,box,box))
#			ptcl.process_inplace("mask.zeroedgefill")
			if options.invert : ptcl.mult(-1.0)
			ptcl["ptcl_source_coord"]=b[:2]
			ptcl["ptcl_source_image"]=microfsp
			ptcl.write_image(outfsp,i+ib*n)

		db.close()
	
		
	E2end(logid)



if __name__ == "__main__":
	main()
