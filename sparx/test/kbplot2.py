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

# Get input values
#alpha = input("Alpha: ")
#K = input("K: ")
#N = input("N: ")
M = 512 # image size
alpha = 1.75
K = 6
N = M*2  # npad*image size
r=M/2
v=K/2.0/N
# set up Kaiser Bessel windows
kb = Util.KaiserBessel(alpha, K, r, v, N)

# Write out Kaiser Bessel I0 window
dx = 0.66007E-03
x = 0.0  #-(K/2+1)
outI0 = file('I0_'+str(alpha)+'.out', 'w')
while x < K/2:
    outI0.write("%g\t%g\n"%(x, kb.i0win(x/N)))
    x += dx
outI0.close()

# Write out Kaisser Bessel sinh window
outsinh = file('sinh_'+str(alpha)+'.out', 'w')
for nu in range (M/2):
    outsinh.write("%g\t%g\n"%(nu, kb.sinhwin(nu)))
outsinh.close()


# Create and write out discretized I0 window
e = EMData()
e.set_size(N, 1, 1)
outI0img = file('I0img_'+str(alpha)+'.out', 'w')
for ix in range(N):
    x = ix
    if (ix > N/2):
        x = ix - N
    val = kb.i0win(x)
    e.set_value_at(ix, 0, 0, val)
    outI0img.write("%g\t%g\n" % (x, val))
outI0img.close()
# Compute and write out FFT'd I0 window
efft = fft(e)
dnu = 1./N
outI0fft = file('I0fft_'+str(alpha)+'.out', 'w')
for inu in range(N/2):
    nu = inu*dnu
    val2 = efft.get_value_at(2*inu, 0, 0)**2 + \
           efft.get_value_at(2*inu+1, 0, 0)**2
    val = sqrt(val2)
    outI0fft.write("%g\t%g\n"%(nu, val/efft.get_value_at(0,0,0)))
outI0fft.close()

# Divide I0 window by sinh window in transform space and write out results
edeconv = filt_kaisersinhinv(e, alpha)
efilt = file('I0filt_'+str(alpha)+'.out', 'w')
for ix in range(N):
    x = ix - N//2
    efilt.write("%g\t%g\n"%(x, edeconv.get_value_at(ix, 0, 0)/edeconv.get_value_at(0,0,0)))
efilt.close()

# Compute and write out discrete sinh window
newft = EMData()
Nft = N+1
odd = True
if N%2 == 0:
    Nft += 1
    odd = False
newft.set_size(Nft, 1, 1)
newft.set_complex(True)
if odd:
    newft.set_fftodd(True)
dnu = 1./N
ftdisc = file("sinh_disc_"+str(alpha)+".out", 'w')
for inu in range(N//2):
    nu = inu*dnu
    val = kb.sinhwin(nu)
    newft.set_value_at(2*inu, 0, 0, val)
    newft.set_value_at(2*inu+1, 0, 0, 0.)
    ftdisc.write("%g\t%g\n"%(nu, val))
ftdisc.close()
# Compute and write out back-FFT'd sinh window
backft = fft(newft)
ftback = file("I0_backft_"+str(alpha)+".out", 'w')
for ix in range(N):
    x = ix
    if ix > N//2:
        x = ix - N
    val = backft.get_value_at(ix,0,0)
    ftback.write("%g\t%g\n"%(x, val/backft.get_value_at(0,0,0)))
ftback.close()
