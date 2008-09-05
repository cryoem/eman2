#!/bin/env python

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

from EMAN2  import *
from sparx  import *

    

vol = EMData()
vol.read_image("../model001.tcp")
#info(e)
mask=model_circle(32,75,75)

volft,kb=prep_vol(vol)

print  "spin"
#   phi theta psi
cp=[12.1,54.7,77.7,0.4,-0.3]
start=[11.2,55.2,77.3,0,0]
proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali_spin(proj,mask,vol,start)

print  "Euler"

proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali(proj,mask,vol,start)

print  "xyz"


proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali_xyz(proj,mask,vol,start)

print " SMALL THETA "

print  "spin"

cp=[12.1,4.7,77.7,0.4,-0.3]
start=[2.2,5.2,77.3,0,0]
proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali_spin(proj,mask,vol,start)

print  "Euler"
proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali(proj,mask,vol,start)


print  "xyz"

proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali_xyz(proj,mask,vol,start)

print "  THETA 90"

print  "spin"

cp=[12.1,90.1,77.7,0.4,-0.3]
start=[2.2,93,77.3,0,0]
proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali_spin(proj,mask,vol,start)

print  "Euler"
proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali(proj,mask,vol,start)


print  "xyz"

proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali_xyz(proj,mask,vol,start)

print  " END "
cp=[12.1,54.7,77.7,0.4,-0.3]
proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali(proj,mask,vol,[10.2,56.2,76.3,0,0])


cp=[12.1,54.7,77.7,0.4,-0.3]
proj=prgs(volft, kb, cp)
proj.write_image("proj.spi")
print  ref_3d_ali_xyz(proj,mask,vol,[10.2,56.2,76.3,0,0])
