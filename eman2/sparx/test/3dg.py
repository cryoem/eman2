#!/usr/bin/env python

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

def cut_pie(a, alpha = 45.0):
	f = fft(a)
	nx = f.get_xsize()
	ny = f.get_ysize()
	from math import atan2, pi
	zer = f.get_value_at(0,0)
	for j in xrange(ny):
		if(j <= ny//2): jy = ny
		else:           jy = j - ny
		for i in xrange(0,nx,2):
			ix = i/2
			phi = atan2(abs(jy), ix)/pi*180
			if(phi <45.0):
				f.set_value_at(i, j, 0.0)
				f.set_value_at(i+1, j, 0.0)
	f.set_value_at(0,0,zer)
	return  fft(f)
			

#image = EMData()
#image.read_image("tf2d0001.tfc")
vol=EMData()
vol.read_image("model001.tcp")
#vol=vols.FourInterpol(256,256,256)
info(vol)


cp=[0,70,0,0.,0.]
volft,kb=prep_vol(vol)
proj=prgs(volft, kb, cp)
nx = proj.get_xsize()
ny = proj.get_ysize()
info(proj)
proj.write_image("proj.hdf")
cp1 = cut_pie(proj)
cp2 = cut_pie(rot_shift2D(proj,-35.))
cp3 = rot_shift2D(cut_pie(rot_shift2D(proj,-35.)), 35.0)
cp1.write_image("pcut.hdf")
cp2.write_image("pcut2.hdf")
cp3.write_image("pcut3.hdf")


from sys import exit
#exit()

st = Util.infomask(proj, None, True)
q = 2.0
t = []
n = 100
for i in xrange(n):
	t.append( proj + filt_gaussl(model_gauss_noise(q*st[1], nx, ny), 0.2))

rssnr, rsumsq, rvar, ssnr, sumsq, var = ssnr2d(t)
pwi = rops_table(proj)
pwn = rops_table(filt_gaussl(model_gauss_noise(q*st[1], nx, ny), 0.2))
fr = []
ps = rot_avg_table(ssnr)
print  len(ps),len(pwn)
ps.extend([0.0])
print  len(ps),len(pwi)
for i in xrange(len(ps)-1):
	fr.append(n*pwi[i]/pwn[i])
	#print  i,ps[i]
	ps[i] = max(0.0, ps[i])

write_text_file([rssnr, rsumsq, rvar, pwi, pwn, fr, ps], "sss.txt")
ssnr.write_image("ssnr.hdf")
sumsq.write_image("sumsq.hdf")
var.write_image("var.hdf")


st = Util.infomask(cp1, None, True)
q = 2.0
t = []
n = 100
for i in xrange(n):
	if i%2: t.append( cp1 + filt_gaussl(model_gauss_noise(q*st[1], nx, ny), 0.2))
	else: t.append( cp2 + filt_gaussl(model_gauss_noise(q*st[1], nx, ny), 0.2))

rssnr, rsumsq, rvar, ssnr, sumsq, var = ssnr2d(t)
pwi = rops_table((cp1+cp2)/2)
pwn = rops_table(filt_gaussl(model_gauss_noise(q*st[1], nx, ny), 0.2))
fr = []
ps = rot_avg_table(ssnr)
ps.extend([0.0])
for i in xrange(len(pwi)):
	fr.append(n*pwi[i]/pwn[i])
	ps[i] = max(0.0, ps[i])

write_text_file([rssnr, rsumsq, rvar, pwi, pwn, fr, ps], "2sss.txt")
ssnr.write_image("2ssnr.hdf")
sumsq.write_image("2sumsq.hdf")
var.write_image("2var.hdf")



st = Util.infomask(cp1, None, True)
q = 2.0
t = []
n = 100
for i in xrange(n):
	if i%2: t.append( cp1 + filt_gaussl(model_gauss_noise(q*st[1], nx, ny), 0.2))
	else: t.append( cp3 + filt_gaussl(model_gauss_noise(q*st[1], nx, ny), 0.2))

rssnr, rsumsq, rvar, ssnr, sumsq, var = ssnr2d(t)
pwi = rops_table((cp1+cp3)/2)
pwn = rops_table(filt_gaussl(model_gauss_noise(q*st[1], nx, ny), 0.2))
fr = []
ps = rot_avg_table(ssnr)
ps.extend([0.0])
for i in xrange(len(pwi)):
	fr.append(n*pwi[i]/pwn[i])
	ps[i] = max(0.0, ps[i])

write_text_file([rssnr, rsumsq, rvar, pwi, pwn, fr, ps], "3sss.txt")
ssnr.write_image("3ssnr.hdf")
sumsq.write_image("3sumsq.hdf")
var.write_image("3var.hdf")

