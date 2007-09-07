#!/usr/bin/env python

#
# Author: David Woolford, 09/07/2007 (woolford@bcm.edu)
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
from math import *
import os
import sys


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input particles> <sim mx> <output>
	Produces class averages """
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--iter", type="int", help="The number of iterations to perform, default is 1", default=1)
	parser.add_option("--hard", type="float", help="The quality metric threshold. Default is off (0)", default=0.0)
	
	(options, args) = parser.parse_args()
	
	if len(args)<3 : parser.error("Input and output files required")

	# check to see if the image exists
	for i in range(0,2):
		if not os.path.exists(args[i]):
			parser.error("File %s does not exist" %args[i])
	
		
	total_particle_images = EMUtil.get_image_count(args[0])
	total_sim_images =  EMUtil.get_image_count(args[1])
	(xsize, ysize ) = gimme_image_2dimensions(args[1]);
	
	(rxsize, rysize ) = gimme_image_2dimensions(args[0]);
	
	if ( ysize != total_particle_images ):
		parser.error("Error, you naughty naughty EMAN2 user")
	
	print "Total images are %d simmx dims are %d %d and contains %d images" %(total_particle_images, xsize, ysize, total_sim_images)
	
	simmx = EMData()
	# just get the images scores
	simmx.read_image(args[1], 0)
	alix = EMData()
	alix.read_image(args[1],1)
	aliy = EMData()
	aliy.read_image(args[1],2)
	aliaz = EMData()
	aliaz.read_image(args[1],3)
	
	classification = []
	
	for i in xrange(0,ysize):
		max_score = simmx.get_value_at(0,i)
		max_idx = 0
		for j in xrange(1,xsize):
			new_score = simmx.get_value_at(j,i)
			if ( new_score < max_score ):
				max_score = new_score
				max_idx = j
		print max_idx
		print max_score
		classification.append(max_idx)
	
	for idx in xrange(0,xsize):
		average = EMData()
		average.set_size(rxsize,rysize,1)
		average.to_zero()
		ptcl_repr = 0.0
		for i in xrange(0,total_particle_images):
			if ( classification[i] == idx ):
				t3d = Transform3D(EULER_EMAN,aliaz.get_value_at(idx),0,0)
				t3d.set_posttrans(alix.get_value_at(idx,i), aliy.get_value_at(idx,i))
				image = EMData()
				image.read_image(args[0], i)
				image.rotate_translate(t3d)
				average.add(image)
				ptcl_repr += 1.0
		
		if ( ptcl_repr != 0 ):
			average.mult(1.0/ptcl_repr)
			
		average.set_attr("ptcl_repr", ptcl_repr)
		average.write_image("classes.img",-1)
	
	#particles = EMData(args[0]);
	#simmx = EMData(args[1]);
	
if __name__ == "__main__":
    main()