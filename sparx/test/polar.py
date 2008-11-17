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

from EMAN2 import *
from sparx  import *

numr=Numrinit(5,5,1,'f')
print  numr
wr=ringwe(numr,'f')
#print  wr

#nring = len(numr)/3
#for i in range(0,nring):
#  print  i+1,numr[2+i*3],wr[i]

e = EMData()
e.read_image("test.spi")

#info(e)


f = EMData()
f.read_image("test30.spi")

#info(f)


circ1=Util.Polar2Dm(e,33.,33.,numr,'f')
#info(circ1)
#circ1.print_image()
#circ=Util.Polar2D(e,numr,'f')
#info(circ)
#circ.print_image()

Util.Frngs(circ1, numr)
#circ.print_image()

circ2=Util.Polar2Dm(f,33.,33.,numr,'f')
#info(circ2)

Util.Frngs(circ2, numr)
#circ2.print_image()

Applyws(circ2, numr,wr)
#circ2.print_image()

neg=0
qn,tot,neg=Util.Crosrng_e(circ2,circ1,numr,neg)

print  qn,tot,ang_n(tot,'f', numr[len(numr)-1])

retvals = Util.Crosrng_ms(circ2,circ1,numr)
qn = retvals["qn"]
tot = retvals["tot"]
qm = retvals["qm"]
tmt = retvals["tmt"]
print  qn,tot,ang_n(tot,'f', numr[len(numr)-1])
print  qm,tmt,ang_n(tmt,'f', numr[len(numr)-1])

retvals = Util.Crosrng_msr(circ2,circ1,numr)
qn = retvals["qn"]
tot = retvals["tot"]
qm = retvals["qm"]
tmt = retvals["tmt"]
print  qn,tot,ang_n(tot,'f', numr[len(numr)-1])
print  qm,tmt,ang_n(tmt,'f', numr[len(numr)-1])

line = Util.Crosrng_msg(circ2,circ1,numr)
print_col(line)
M=lin1.get_xsize()
# no pading
npad=1
N=M*npad
# support of the window
K=6
alpha=1.75
r=M/2
v=K/2.0/N
kb = Util.KaiserBessel(alpha, K, r, K/(2.*N), N)
params = {"filter_type" : Processor.fourier_filter_types.KAISER_SINH_INVERSE,"alpha" : alpha, "K":K,"r":r,"v":v,"N":N}
q=Processor.EMFourierFilter(line,params)
#  remember that x is counted to zero and is divided by two in get_pixel_conv !!
x=2.2
q.get_pixel_conv(x,0,0,kb)
