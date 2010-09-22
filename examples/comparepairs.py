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

# comparepairs.py  08/09/2005  Steven Ludtke

from EMAN2 import *
import sys
from math import *

cmps=[]

n=EMUtil.get_image_count(sys.argv[1])

# read pairs of images and calculate a variety of
# similarity metrics between the pairs
for i in range(0,n,2):
	l=EMData.read_images(sys.argv[1],(i,i+1))
	l[1].process_inplace("normalize.edgemean")
	l[1].process_inplace("filter.lowpass.gauss",{"cutoff_abs":.1})

	cmps.append((
		l[0].cmp("phase",l[1],{}),
		sqrt(l[0].cmp("optvariance",l[1],{"radweight":1})/l[0].get_xsize()),
		sqrt(l[0].cmp("optvariance",l[1],{})),
		l[0].cmp("dot",l[1],{"normalize":1}),
		l[0].cmp("frc",l[1],{}),
		l[0].cmp("sqeuclidean",l[1],{})))

# calculate the mean for each measure
sm=[0,0,0,0,0,0]
for i in cmps:
	for j,k in enumerate(i):
		sm[j]+=k

sm=[j/len(cmps) for j in sm]

# now print out each similarity measure for each pair
for i,j in enumerate(cmps):
	print "%d.\t"%i,
	for k,l in enumerate(j):
#		print "%1.4f\t"%-(l-sm[k]),
		print "%1.4f\t"%l,
	print ""
