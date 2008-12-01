#!/usr/bin/env python

#
# Author: Steven Ludtke, 11/30/2008 (ludtke@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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
from optparse import OptionParser
import random
from math import *
import os
import sys
from e2simmx import cmponetomany

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] 
	Initial model generator"""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--input", dest="input", default=None,type="string", help="The name of the image containing the particle data")
	parser.add_option("--iter", type = "int", default=8, help = "The total number of refinement iterations to perform")
	parser.add_option("--tries", type="int", default=1, help="The number of different initial models to generate in search of a good one")
	parser.add_option("--sym", dest = "sym", help = "Specify symmetry - choices are: c<n>, d<n>, h<n>, tet, oct, icos",default="c1")
	parser.add_option("--verbose","-v", dest="verbose", default=False, action="store_true",help="Toggle verbose mode - prints extra infromation to the command line while executing")

	(options, args) = parser.parse_args()


	ptcls=EMData.read_images(options.input)
	if not ptcls or len(ptcls)==0 : parser.error("Bad input file")
	boxsize=ptcls[0].get_xsize()

	# angles to use for refinement
	sym_object = parsesym(options.sym)
	orts = sym_object.gen_orientations("eman",{"delta":7.5})

	logid=E2init(sys.argv)

	# We make one new reconstruction for each loop of t 
	for t in range(options.tries):
		threed=[make_random_map(boxsize)]		# initial model
		apply_sym(threed[0],options.sym)		# with the correct symmetry
		
		# This is the refinement loop
		for it in range(options.iter):
			projs=[threed[j].project("standard",ort) for ort in orts]		# projections
			
			for i in range(len(ptcls)):
				sim=cmponetomany(projs,ptcls[i],align="rotate_translate_flip",alicmp=("frc",{}))
				n=sim.index(min(sim))
				ptcls[i]["xform.projection"]=orts[n]	# best orientation set in the original particle
				
			
			
	E2end(logid)


def make_random_map(boxsize):
	"""This will make a map consisting of random noise, low-pass filtered and center-weighted for use
	as a random starting model in initial model generation. Note that the mask is eliptical and has random aspect."""
	
	ret=EMData(boxsize,boxsize,boxsize)
	ret.process_inplace("testimage.noise.gauss",{"mean":0.02,"sigma":1.0})
	ret.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})
#	ret.process_inplace("mask.gaussian",{"inner_radius":boxsize/3.0,"outer_radius":boxsize/12.0})
	ret.process_inplace("mask.gaussian.nonuniform",{"radius_x":boxsize/random.uniform(2.0,5.0),"radius_y":boxsize/random.uniform(2.0,5.0),"radius_z":boxsize/random.uniform(2.0,5.0)})
	
	return ret
	
def apply_sym(data,sym):
	"""applies a symmetry to a 3-D volume in-place"""
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(sym)
	ref=data.copy()
	for i in range(1,nsym):
		dc=ref.copy()
		dc.transform(xf.get_sym(sym,i))
		data.add(dc)
	data.mult(1.0/nsym)	
	

if __name__ == "__main__":
    main()

