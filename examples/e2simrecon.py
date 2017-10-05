#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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


from EMAN2 import *
from math import *
import os
import sys
import traceback
import random
import numpy as np

def run(command):
	"Mostly here for debugging, allows you to control how commands are executed (os.system is normal)"

	print "{}: {}".format(time.ctime(time.time()),command)

	ret=launch_childprocess(command)

	# We put the exit here since this is what we'd do in every case anyway. Saves replication of error detection code above.
	if ret !=0 :
		print "Error running: ",command
		sys.exit(1)

	return


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <volume in> <volume out>

Simulates the effects of a 3D reconstruction by including noise and rotational uncertainty """

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--apix", "-A", type=float, help="A/voxel", default=1.0)
	parser.add_argument("--res", "-R", type=float, help="Resolution in A. This is not a filter, this is a target resolution to achieve a 0.143 FSC with noise present",default=2.8)
	parser.add_argument("--anglesigma", type=float, help="Uncertainty in orientation determination in degrees")
	parser.add_argument("--sym", dest = "sym", default="c1",help = "Specify symmetry - choices are: c<n>, d<n>, tet, oct, icos.")
	parser.add_argument("--threads", default=4 ,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")

	vol=EMData(args[0],0,True)
	nz=vol["nz"]
	if nz==1 :
		print "Input must be a volume"
		sys.exit(1)
		
	da=360.0/nz
	
	com="e2project3d.py {} --outfile simproj.hdf -f --orientgen eman:delta={}:inc_mirror=0:random_phi=1 --sym {} --parallel thread:{}".format(args[0],da,options.sym,options.threads)
	run(com)
	
	n=EMUtil.get_image_count("simproj.hdf")
	for i in xrange(n):
		h=EMData(args[0],i,True)
		xf=h["xform.projection"]
		k=xf.get_rotation("eman")
		k["alt"]+=random.gauss(0,options.anglesigma)
		k["az"]+=random.gauss(0,options.anglesigma)
		k["phi"]+=random.gauss(0,options.anglesigma)
		xf.set_rotation(k)
		h["xform.projection"]=xf
		h.write_image(args[0],i,IMAGE_UNKNOWN,True)
	
