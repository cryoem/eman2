#!/usr/bin/env python
#
# Author: Steve Ludtke 02/13/20 (sludtke@bcm.edu)
# Copyright (c) 2014- Baylor College of Medicine
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


from builtins import range
from EMAN2 import *
from math import *
import os
import sys
import time
import traceback
from EMAN2star import StarFile


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <relion STAR file> <output set name>

This program will take a STAR file representing a subset of data from a Relion project and convert it into an EMAN2
set, IF the EMAN2 project contains the same micrographs and exact same boxed out particles as the Relion project
which produced the STAR file. e2refinetoeman.py can be used to convert an entire Relion project into an EMAN2
project. This program can then be used in this established project to copy information about particle subsets
selected in Relion into EMAN2.
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#options associated with e2refine.py
	parser.add_pos_argument(name="star_file",help="Select STAR file", default="", guitype='filebox', browser="EMParticlesEditTable(withmodal=True,multiselect=False)",  row=6, col=0,rowspan=1, colspan=3)
	parser.add_pos_argument(name="output_set_name",help="Name of output set", default="from_relion", guitype='str', row=7, col=0,rowspan=1, colspan=1)
	parser.add_argument("--origname",action="store_true",help="Adds the original STAR name as a comment on each image",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.verbose>0 : print("Parsing STAR file")
	star=StarFile(args[0])

	logid=E2init(sys.argv,options.ppid)
	
	if "rlnMicrographName" in star : difkey="rlnMicrographName"
	else: difkey="rlnDefocusU"
	
	if '.lst' in args[1]:
		if not "sets/" in args[1]: args[1]=f"sets/{args[1]}"
	else:
		args[1]=f"sets/{args[1]}.lst"
	
	try: os.unlink(args[1])
	except: pass
	out=LSXFile(args[1])
	
	orig=None
	for i in range(len(star["rlnImageName"])):
		rlnname=star["rlnImageName"][i]
		name=rlnname.split("@")[1]
		imgnum=int(rlnname.split("@")[0])-1

		emanname=os.path.split(name)[1]
		emanname=os.path.join("particles",os.path.splitext(emanname)[0]+".hdf")

		if not os.path.exists(emanname) :
			badname=emanname
			# this won't always work right, but it's worth a try
			emanname=os.path.split(name)[1]
			emanname=os.path.join("particles",name.split("/")[-2]+"-"+os.path.splitext(emanname)[0]+".hdf")
			if not os.path.exists(emanname) :
				raise Exception(f"Cannot find:{emanname} or {badname} (from {rlnname})")

		if options.origname: orig="# "+rlnname
		out.write(-1,imgnum,emanname,orig)

	E2end(logid)

if __name__ == "__main__":
    main()
