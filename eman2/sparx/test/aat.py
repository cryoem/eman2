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

from sys import exit

#vol=get_image("model001.tcp")
vol=Util.window(get_image("/mnt/shared/pawel/VAR/VOL.hdf"),95,95,95,0,0,0)
mask = model_circle(45,95,95)
#vol=get_image("/mnt/shared/pawel/VAR/VOL.hdf")
#mask = model_circle(45,96,96)
#M=75
#msk=model_circle(30,M,M,1)
#vol=model_circle(0.6,M,M,M)
#  params:  phi, theta, psi, sx, sy
volft,kb=prep_vol(vol)



params = [0.0,0.0,0.0,0.0,0.0]
proj=prgs(volft,kb,params)
info(proj)
drop_image(proj,"xxx1.hdf")



prl = project(vol,params,1000)
info(prl)
drop_image(proj,"xxx2.hdf")

print  ccc(proj,prl)
f = fsc(proj,prl,1.,"zero.txt")

d = im_diff(proj,prl)
p = periodogram(d[0])*mask
drop_image(p,"xxx5.hdf")


params = [19.0,77.0,85.0,0.0,0.0]
proj=prgs(volft,kb,params)
info(proj)
drop_image(proj,"xxx3.hdf")



prl = project(vol,params,1000)
info(prl)
drop_image(proj,"xxx4.hdf")



print  ccc(proj,prl)
f = fsc(proj,prl,1.,"weird.txt")




d = im_diff(proj,prl)
p = periodogram(d[0])*mask
drop_image(p,"xxx6.hdf")



exit()



params = [19.0,77.0,0.0,0.0,0.0]
proj0=prgs(volft,kb,params)
drop_image(proj0,"xxx0.spi")
a1 = 33.
x1 = 2.0
y1 = -3.0
o = rtshg(proj,a1,x1,y1)
drop_image(o,"xxx2.spi")
a2 = -10.
x2 = 3.0
y2 = 2.0
t = rtshg(o,a2,x2,y2)
drop_image(t,"xxx3.spi")

ia2,ix2,iy2,s = inverse_transform2(a2,x2,y2)


ia1,ix1,iy1,s = inverse_transform2(a1,x1,y1)
io = rtshg(o,ia1,ix1,iy1)

print ccc(io,proj,msk)

# inverse transfromation in one step
ta,tx,ty,s=compose_transform2(ia2,ix2,iy2,1.0,ia1,ix1,iy1,1.0)
print  ta,tx,ty

q = rtshg(t,ta,tx,ty)
drop_image(q,"xxx4.spi")
print ccc(q,proj,msk)

# inverse trasformation in one step with additional rotation by psi added

tta,ttx,tty,s=compose_transform2(ta,tx,ty,1.0,-85.,0,0,1.0)
print  tta,ttx,tty

p = rtshg(t,tta,ttx,tty)
drop_image(p,"xxx5.spi")
print ccc(p,proj0,msk)

# work until here


# Now invert and transfer all these parameters to the projection of the volume....
ia,ix,iy,s=inverse_transform2(tta,ttx,tty)
print  ia,ix,iy
params = [19.0,77.0,ia,ix,iy]
projt=prgs(volft,kb,params)
drop_image(projt,"xxx7.spi")
print ccc(t,projt,msk)
