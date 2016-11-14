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

    


vol=getImage("../model001.tcp")
#info(e)
vol = vol.FourInterpol(320,320,320)

dtheta=12
mask3D=[]#(32,75,75,75)
mask=model_circle(135,75,75,1)#(32,75,75,1)
start=0
end=3
first_r =1
last_r = 135#32
rstep=1
xrng=2
yrng=2
step=1
newp = proj_ali(vol, mask3D, "grj/bprj{*****}.tcp", start, end, "grj/prj{*****}.tcp", first_r,last_r,rstep,xrng,yrng,step,dtheta)
print  newp
#newp=[0.0,0.0,0.0,0.0,0.0]
#print  newp
#volft,kb=prep_vol(vol)
#proj=prgs(volft, kb, newp)
#print  ref_3d_ali(proj,mask,vol,newp)
refp = ref_ali(vol, mask3D, mask, "grj/bprj{*****}.tcp", start, end, "grj/prj{*****}.tcp", dtheta, newp)
print  refp
