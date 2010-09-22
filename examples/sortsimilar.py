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

# sortsimilar.py  08/09/2005  Steven Ludtke
# This will sort projection/image pairs based on their mutual similarity
# more similar pairs will be placed in the output file first

from EMAN2 import *
import sys

cmps=[]

n=EMUtil.get_image_count(sys.argv[1])

# read pairs of images and calculate similarity
# store in a list containing similarity and image number
for i in range(0,n,2):
	l=EMData.read_images(sys.argv[1],(i,i+1))
#	l[1].process_inplace("normalize.edgemean")
#	l[1].process_inplace("filter.lowpass.gauss",{"cutoff_abs":.08})
#	cmps.append((l[0].cmp("phase",l[1],{})+l[0].cmp("optvariance",l[1],{"radweight":1})/l[0].get_xsize(),i))
	cmps.append((l[0].cmp("optvariance",l[1],{"matchamp":1})/l[0].get_xsize(),i))


#	cmps.append((l[0].cmp("dot",l[1],{"normalize":1}),
#		l[0].cmp("frc",l[1],{}),
#		l[0].cmp("optvariance",l[1],{}),
#		l[0].cmp("phase",l[1],{}),
#		l[0].cmp("quadmindot",l[1],{"normalize":1}),
#		l[0].cmp("SqEuclidean",l[1],{})))

# sort the list in ascending order. Smaller similarity values are better
cmps.sort()

# write the images to the output file
for i in cmps:
	l=EMData.read_images(sys.argv[1],(i[1],i[1]+1))
	l[0].write_image(sys.argv[2],-1)
	l[1].write_image(sys.argv[2],-1)
