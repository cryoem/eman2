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
last_ring=125
rstep=1



#e = EMData()
#e.read_image("tf2d0001.tfc")
#mask=model_circle(27,64,64)
#stat=Util.infomask(e,mask)
#ee=(e-stat[0])/stat[1]  #*mask
sx=0
sy=0
ee=model_gauss_noise(1.0,256,256)

for i in range(0,400+1):
      a=i/100.0
      o = rtshg(ee,a,sx,sy)


      #[angs,sxs,sys,mirror,peak]=ormq(ee,o,first_ring,last_ring,1,2,2,1,"F")
      #[angsi,sxsi,sysi,mirror,peak]=ormqi(ee,o,first_ring,last_ring,1,2,2,1,"F")
      #[angsl,sxsl,sysl,mirror,peak]=ormqlo(ee,o,first_ring,last_ring,1,2,2,1,"F")
      [angs,sxs,sys,mirror,peak]=ormq(ee,o,first_ring,last_ring,1,0,0,1,"F")
      [angsi,sxsi,sysi,mirror,peak]=ormqi(ee,o,first_ring,last_ring,1,0,0,1,"F")
      [angsl,sxsl,sysl,mirror,peak]=ormql(ee,o,first_ring,last_ring,1,0,0,1,"F")
      print  a,"  ",sx,"  ",sy,"  ",angs,"  ",sxs,"  ",sys,"  ",angsi,"  ",sxsi,"  ",sysi,"  ",angsl,"  ",sxsl,"  ",sysl
