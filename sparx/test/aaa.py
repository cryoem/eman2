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
import sparx.libpy

# For convenience, import all functions into the top-level namespace. 
for nm in dir(sparx.libpy):
    if nm.startswith("__"): continue
    exec "from %s import *" % nm
	    
    
    
out = file("tccc_angle", "w")



e = get_image("tst1.spi")#("tf2d0001.tfc")
mask=model_circle(32,75,75)#(27,64,64)
stat=Util.infomask(e,mask)
ee=(e-stat[0])/stat[1]*mask
drop_image(ee,'rtg1.spi')
#e.set_size(64,64,1)
#e.to_zero()
#e[0,0]=64*64/4
#print  Util.infomask(ee,mask)
#rops_textfile(ee,"epw")
mas=model_circle(30,75,75)#(25,64,64)
info(ee,mas)
inorm = ee.cmp("dot",ee,{"mask":mas,"negative":0})
for i in range(0,0+1):
	  o = ee.rot_scale_trans2D(i*pi/180.0)
	  #print  Util.infomask(o,mask)
	  o = o*mask
	  drop_image(o,'rst2.spi')
	  #rops_textfile(o,"opw4")
	  u = o.rot_scale_trans2D(-i*pi/180.0)
	  u=u*mask
	  stat=Util.infomask(u,mask)
	  print  " STAT T",stat[0],"  ",stat[1],"  ",stat[2],"  ",stat[3]
	  #u=(u-stat[0])/stat[1]
	  drop_image(u,'rst3.spi')
	  fsc(ee,u,1,"fsct")
	  ct= ccc(ee,u,mas)
          d=u-ee
	  info(d,mas)
	  drop_image(d,'rst4.spi')
	  a=Util.im_diff(u,ee,mas)
	  d=a["imdiff"]
	  info(d,mas)
	  tnorm = d.cmp("dot",d,{"mask":mas,"negative":0})
	  drop_image(d,'rst5.spi')

	  o = rtshg(ee,i,0.,0.)
	  #print  Util.infomask(o,mask)
	  o = o*mask
	  drop_image(o,'rtg2.spi')
	  #rops_textfile(o,"opw4")
	  u = rtshg(o,-i,0.,0.)
	  u=u*mask
	  stat=Util.infomask(u,mask)
	  print  " STAT G",stat[0],"  ",stat[1],"  ",stat[2],"  ",stat[3]
	  #u=(u-stat[0])/stat[1]
	  drop_image(u,'rtg3.spi')
	  fsc(ee,u,1,"fscg")
	  cg= ccc(ee,u,mas)
          d=u-ee
	  info(d,mas)
	  drop_image(d,'rtg4.spi')
	  a=Util.im_diff(u,ee,mas)
	  d=a["imdiff"]
	  info(d,mas)
	  gnorm = d.cmp("dot",d,{"mask":mas,"negative":0})
	  drop_image(d,'rtg5.spi')

	  s=ee.copy()
	  s.rotate_translate(-i,0,0,0,0,0,0,0,0)
	  drop_image(s,'sss2.spi')
	  s.rotate_translate(i,0,0,0,0,0,0,0,0)
	  stat=Util.infomask(s,mask)
	  print  " STAT L",stat[0],"  ",stat[1],"  ",stat[2],"  ",stat[3]
          #info(s)
	  #s=(s-stat[0])/stat[1]
	  drop_image(s,'sss3.spi')
	  fsc(ee,s,1,"fscs")
	  cs= ccc(ee,s,mas)
          d=s-ee
	  drop_image(d,'sss4.spi')
	  info(d,mas)
	  a=Util.im_diff(s,ee,mas)
	  d=a["imdiff"]
	  snorm = d.cmp("dot",d,{"mask":mas,"negative":0})
	  info(d,mas)
	  drop_image(d,'sss5.spi')

          print  i, cs, cg, ct, inorm, snorm, tnorm, gnorm
	  if(snorm==0.0):
	    snorm=1.0
	  if(tnorm==0.0):
	    tnorm=1.0
	  out.write("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n" % (i, cs, ct, cg, inorm, snorm, tnorm, gnorm, inorm/snorm, inorm/tnorm, inorm/gnorm))

