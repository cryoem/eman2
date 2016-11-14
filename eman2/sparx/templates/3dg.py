#!/bin/env python

#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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


#  TEST:  generate projections of a model, compute reconstructions and fsc curve.

vol = EMData()
vol.read_image("../test/model001.tcp")
info(vol)
nx = vol.get_xsize()
delta = 10.0

angles = even_angles(delta,0.,90.,0,359.99,"S")

#angles=even_angles(delta_theta,0.,180.,0,359.99/2,"P")

volft,kb = prep_vol(vol)

stack_data = "data.hdf"
import os
os.system("rm -f  "+stack_data)
#print  angles
nangles = len(angles)
#dropSpiderDoc("angles.txt",angles)
ppp = []
s2x=0
s2y=0
from random import random,randint
for i in xrange(nangles):
	s2x = 4.0*randint(-1,1)
	s2y = 4.0*randint(-1,1)
	ppp.append([angles[i][0], angles[i][1], angles[i][2], s2x, s2y])
	projo = prgs(volft, kb, [angles[i][0], angles[i][1], angles[i][2], -s2x, -s2y])
	#apply CTF
	defocus = randint(20,40)*1000.0
	proj = filt_ctf(projo,defocus,2.0,300,2.5,0.1)
	proj += model_gauss_noise(10.0, nx, nx)
	# Set all parameters for the new 2D image
	# three angles and two shifts to zero
	proj.set_attr_dict({'phi':angles[i][0], 'theta':angles[i][1], 'psi':angles[i][2], 's2x':s2x, 's2y':s2y})
	# CTF parameters, if defocus zero, they are undetermined
	proj.set_attr_dict({'defocus':defocus, 'amp_contrast':0.1, 'voltage':300, 'Cs':2.0, 'Pixel_size':2.5, 'B_factor':0.})
	# flags describing the status of the image (1 = true, 0 = false)
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# proj.set_attr_dict({'active':1, 'ctf_applied':0})
	proj.set_attr_dict({'ctf_applied':0})
	proj.write_image(stack_data, i)
del volft
dropSpiderDoc("params.txt",ppp)
del ppp
from sys import exit
snr = 2.0
list_p = range(0,nangles,2)
vol1   = recons3d_4nn_ctf(stack_data, list_p, snr)

list_p = range(1,nangles,2)
vol2   = recons3d_4nn_ctf(stack_data, list_p, snr)

mask3d = model_circle(nx//2-5,nx,nx,nx)

Util.mul_img(vol1, mask3d)
Util.mul_img(vol2, mask3d)
del mask3d
fsc(vol1,vol2,0.5,"tdt.txt")
del vol1, vol2
volt = recons3d_4nn_ctf(stack_data, range(nangles), snr)
dropImage(volt, "volt.spi", "s")   # to be displayed in chimera
