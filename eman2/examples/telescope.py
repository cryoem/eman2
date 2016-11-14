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
from emimage import *
import time

def live():
	b=EMData()
	b.read_image("/dev/video0",0)
	c=EMImage(b)
	for i in range(10):
		b.read_image("/dev/video0",0)
		c.setdata(b)
		time.sleep(1)

def rawavg(noframes):
	b=EMData()
	b.read_image("/dev/video0",0)
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
		b+=a
		if i%10==0 : print i
	return b

def aliavg(noframes):
	b=EMData()
	b.read_image("/dev/video0",0)
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
#		ba=a.align("translational",b,{},"optvariance",{"matchfilt":1})
		ba=a.align("translational",b,{},"dot",{})
		b+=ba
		print i
	return b


def rawframes(noframes,outfile=None):
	ret=[]
	b=EMData()
	b.read_image("/dev/video0",0)
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
		ret.append(a)
		b+=a
	if outfile: write_set(ret,outfile)
	return (ret,b)
		
def rawframesdisk(noframes,outfile):
	for i in range(noframes):
		a=EMData()
		a.read_image("/dev/video0",0)
		a.write_image(outfile,-1)
		if i%10==0: print i

def avgframesdisk(noframes,noavg,outfile):
	for i in range(noframes):
	  a=EMData()
	  a.read_image("/dev/video0",0)
	  for j in range(noavg-1):
		b=EMData()
		b.read_image("/dev/video0",0)
		a+=b
	  a.write_image(outfile,-1)
	  if i%10==0: print i

def write_set(lst,outfile):
	for i in lst:
		i.write_image(outfile,-1)

a=rawframes(10)[0][9]
display(a)
