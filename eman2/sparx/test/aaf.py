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


# this program demonstrates nearly perfect rotational alignment if rings two times oversampled

first_ring=1
last_ring=25
rstep=1



e = EMData()
e.read_image("tf2d0001.tfc")
mask=model_circle(27,64,64)
stat=Util.infomask(e,mask)
ee=(e-stat[0])/stat[1]  #*mask
mask=model_circle(25,64,64)
sx=0
sy=0

line = EMData()

for i in range(0,30+1):
  for ix in range(0,20+1):
    for iy in range(0,10+1):
      sx=ix/20.0
      sy=iy/10.0
      a=i/10.0
      o = rtshg(ee,a,sx,sy)


      [angs,sxs,sys,mirror,peak]=ormq(o,ee,first_ring,last_ring,1,2,2,1,"F")
      #[angsi,sxsi,sysi,mirror,peak]=ormqi(o,ee,first_ring,last_ring,1,2,2,1,"F")
      [angsl,sxsl,sysl,mirror,peak]=ormql(o,ee,first_ring,last_ring,1,2,2,1,"F")
      if(angs>180):
        angs=angs-360
      #if(angsi>180):
      #  angsi=angsi-360
      if(angsl>180):
        angsl=angsl-360
	

      ia=inverse_transform2(a,sx,sy)
      #print  ia[0],"  ",ia[1],"  ",ia[2],"  ",angs,"  ",sxs,"  ",sys,"  ",angsi,"  ",sxsi,"  ",sysi,"  ",angsl,"  ",sxsl,"  ",sysl
      #print  ia[0],"  ",ia[1],"  ",ia[2],"  ",angs,"  ",sxs,"  ",sys,"  ",angsl,"  ",sxsl,"  ",sysl
      q1=rtshg(o,ia[0],ia[1],ia[2])
      c1=ccc(ee,q1,mask)
      d1=q1.cmp("SqEuclidean", ee, {"mask":mask})
      q1=rtshg(o,angs,sxs,sys)
      c2=ccc(ee,q1,mask)
      d2=q1.cmp("SqEuclidean", ee, {"mask":mask})
      #q1=rtshg(o,angsi,sxsi,sysi)
      #c3=ccc(ee,q1,mask)
      #d3=q1.cmp("SqEuclidean", ee, {"mask":mask})
      q1=rtshg(o,angsl,sxsl,sysl)
      c4=ccc(ee,q1,mask)
      d4=q1.cmp("SqEuclidean", ee, {"mask":mask})
      #print  c1,d1,c2,d2,c3,d3,c4,d4
      ib=compose_transform2(a,sx,sy,1.0,angs,sxs,sys,1.0)
      derrs=2*last_ring*abs(sin(ib[0]*pi/180./2.))+sqrt(ib[1]*ib[1]+ib[2]*ib[2])
      ib=compose_transform2(a,sx,sy,1.0,angsl,sxsl,sysl,1.0)
      derrl=2*last_ring*abs(sin(ib[0]*pi/180./2.))+sqrt(ib[1]*ib[1]+ib[2]*ib[2])
      #print  ib
      print  a,sx,sy,c1,d1,c2,d2,derrs,c4,d4,derrl
