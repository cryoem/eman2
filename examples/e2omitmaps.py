#!/usr/bin/env python

#
# Author: Steve Ludtke 05/23/2016 (sludtke@bcm.edu)
# Copyright (c) 2013- Baylor College of Medicine
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
from sys import argv,exit

if len(argv)<4 :
	print "e2omitmap.py <map> <nseg> <mass>"
	sys.exit(1)

model=EMData(argv[1],0)
nseg=int(argv[2])
mass=float(argv[3])

model.process_inplace("normalize.bymass",{"thr":1,"mass":mass})
seg=model.process("segment.kmeans",{"ampweight":1,"nseg":nseg,"thr":0.7})
for i in range(nseg):
	seg2=seg.process("threshold.binaryrange",{"low":i-0.1,"high":i+0.1})	# by subtracting 1, we don't remove anything from the first map
	seg2.process_inplace("math.linear",{"scale":-1.0,"shift":1.0})
	model2=model*seg2
	model2.write_image("omit_{:02d}.hdf".format(i+1),0)
