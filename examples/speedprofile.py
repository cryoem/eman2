#!/usr/bin/env python

#
# Author: Steven Ludtke, 01/06/2010 (sludtke@bcm.edu)
# Copyright (c) 2000-2010 Baylor College of Medicine
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

# this program will run speedtests for all possible box sizes over a range
# and measure the relative program performance, it takes a LOOONG time to run

from EMAN2 import *
import sys
import time
import random

NTT=100

def init(SIZE=96,NTT=100):
	data=[]
	for i in xrange(NTT):
		data.append(test_image(0,size=(SIZE,SIZE)))
		data[i].transform(Transform({"type":"2d","alpha":random.uniform(0,360.0),"tx":random.uniform(-5.0,5.0),"ty":random.uniform(-5.0,5.0)}))
		data[i].add(test_image(1,size=(SIZE,SIZE)))
		data[i].process_inplace('normalize.circlemean')
		data[i].process_inplace('mask.sharp', {'outer_radius':data[i].get_xsize()/2})

	return data

def catime(SIZE=96,NTT=100):
	data=init(SIZE,NTT)
	ref=test_image(0,size=(SIZE,SIZE))

	start=time.time()
	for i in xrange(NTT):
		x=data[i].align("rotate_translate_flip",ref,{"maxshift":6.0},"dot",{"normalize":0})
		x=x.align("refine",ref,{"maxshift":6.0},"dot",{"normalize":0})
		y=x.cmp("phase",ref)

	return (time.time()-start)/NTT

print "establishing baseline"
base=catime(SIZE=32,NTT=10000)

print "testing"
out=file("profile.txt","w")
for i in xrange(32,513):
	t=catime(i,16000/i)
	print	"%d\t%1.2f\t%1.3f"%(i,t/base,t)
	out.write("%d\t%1.3f\n"%(i,t/base))
	out.flush()

