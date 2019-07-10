#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

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


from past.utils import old_div
from builtins import range
from EMAN2 import *
from math import *
import os
import sys
import traceback
import random
import numpy as np

def run(command):
	"Mostly here for debugging, allows you to control how commands are executed (os.system is normal)"

	print("{}: {}".format(time.ctime(time.time()),command))

	ret=launch_childprocess(command)

	# We put the exit here since this is what we'd do in every case anyway. Saves replication of error detection code above.
	if ret !=0 :
		print("Error running: ",command)
		sys.exit(1)

	return


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <volume in> <volume out>

Simulates the effects of a 3D reconstruction by including noise and rotational uncertainty """

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--apix", "-A", type=float, help="A/voxel", default=0.0)
	parser.add_argument("--fsc", type=str, help="A text file containing a FSC curve to simulate",default=None)
	parser.add_argument("--anglesigma", type=float, help="Uncertainty in orientation determination in degrees")
	parser.add_argument("--sym", dest = "sym", default="c1",help = "Specify symmetry - choices are: c<n>, d<n>, tet, oct, icos.")
	parser.add_argument("--threads", default=4 ,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")

	if options.fsc==None:
		print("ERROR: must specify FSC curve")
		sys.exit(1)
	fsc=XYData()
	fsc.read_file(options.fsc)

	vol=EMData(args[0],0,True)
	if options.apix<=0 : options.apix=vol["apix_x"]

	nz=vol["nz"]
	if nz==1 :
		print("Input must be a volume")
		sys.exit(1)

	logger=E2init(sys.argv,options.ppid)


#	da=2*360.0/nz

	com="e2project3d.py {} --outfile simproj.hdf -f --orientgen rand:phitoo=1:n={} --sym {} --parallel thread:{}".format(args[0],nz*3,options.sym,options.threads)
	run(com)

	n=EMUtil.get_image_count("simproj.hdf")
	for i in range(n):
		h=EMData("simproj.hdf",i,True)
		xf=h["xform.projection"]
		k=xf.get_rotation("eman")
		k["alt"]+=random.gauss(0,options.anglesigma)
		k["az"]+=random.gauss(0,options.anglesigma)
		k["phi"]+=random.gauss(0,options.anglesigma)
		xf.set_rotation(k)
		h["xform.projection"]=xf
		h["ptcl_repr"]=1
		h.write_image("simproj.hdf",i,IMAGE_UNKNOWN,True)

	com="e2make3dpar.py --input simproj.hdf --sym {} --output nonoise.hdf --keep 1 --pad -1 --mode gauss_var --threads {}".format(options.sym,options.threads)
	run(com)

	recon=EMData("nonoise.hdf")
	reconf=recon.do_fft()

	ds=old_div(1.0,(options.apix*nz))
	pspec=reconf.calc_radial_dist(int(nz*1.8),0,1,1)
#	pspec=reconf.calc_radial_dist(len(fsc),fsc.get_x(0)/ds,(fsc.get_x(1)-fsc.get_x(0))/ds,1)
#	print pspec
	fscl=[min(max(fsc.get_yatx(i*ds),.001),0.9999) for i in range(len(pspec))]
#	print fscl
	noise=[old_div(sqrt(pspec[i]*(1.0-fscl[i])/(fscl[i])),(4000.0)) for i in range(len(pspec))]
	print(noise)

	noisemap=EMData(nz,nz,nz)
	noisemap.process_inplace("testimage.noise.gauss")
	noisemap.process_inplace("filter.radialtable",{"table":noise})

	reconc=recon.copy()
	reconc.add(noisemap)
	reconc.write_image("{}_even.hdf".format(args[1].rsplit(".",1)[0]),0)

	noisemap=EMData(nz,nz,nz)
	noisemap.process_inplace("testimage.noise.gauss")
	noisemap.process_inplace("filter.radialtable",{"table":noise})

	reconc=recon.copy()
	reconc.add(noisemap)
	reconc.write_image("{}_odd.hdf".format(args[1].rsplit(".",1)[0]),0)


	E2end(logger)


if __name__=="__main__":
	main()
