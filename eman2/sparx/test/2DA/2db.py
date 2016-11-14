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

#  WORKING VERSION!

vol = EMData()
vol.read_image("../model001.tcp")
#info(e)

delta_theta=15

angles=even_angles(delta_theta,0.,180.)
#angles = [[72,27,0],[52,60,45],[180+12,180-12,180+125]]

volft,kb=prep_vol(vol)

stack_data="data.hdf"

#print  angles
nangles = len(angles)
dropSpiderDoc("angles",angles)
ppp = []
seed()
for i in xrange(nangles):
	sx=3.0*random()-2  #int(5.0*random()+0.99999)-3
	sy=3.0*random()-2
	ppp.append([angles[i][0],angles[i][1],angles[i][2],sx,sy])
	proj=prgs(volft, kb, [angles[i][0],angles[i][1],angles[i][2],-sx,-sy])
	# Set all parameters for the new 2D image
	# three angles and two shifts
	sx=0
	sy=0
	proj.set_attr_dict({'phi':angles[i][0], 'theta':angles[i][1], 'psi':angles[i][2], 'sx':sx, 'sy':sy, 'mirror': 0})
	# CTF parameters, if defocus zero, they are undetermined
	proj.set_attr_dict({'defocus':0.0, 'amp_contrast':0.1, 'voltage':200, 'Cs':2.0, 'pixel':2.2})
	# flags describing the status of the image (1 = true, 0 = false)
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# proj.set_attr_dict({'active':1, 'ctf_applied':0})
	proj.set_attr_dict({'ctf_applied':0})
	proj.write_image(stack_data, i)
#dropSpiderDoc("params",ppp)
#del ppp

#  BEGIN

stack_data="data.hdf"
nima = EMUtil.get_image_count(stack_data)
temp = EMData()
temp.read_image(stack_data,0)
# prepare KB interpolants
nx=temp.get_xsize()
kb = kbt(nx)

# read images and prepare them for gridding
data = []
for im in xrange(nima):
	if(im>0):
		temp = EMData()
		temp.read_image(stack_data,im)
	phi = temp.get_attr('phi')
	theta = temp.get_attr('theta')
	psi = temp.get_attr('psi')
	sx =  temp.get_attr('sx')
	sy =  temp.get_attr('sy')
	mn =  temp.get_attr('mirror')
	data.append(prepg(temp,kb))
	data[im].set_attr_dict({'phi':phi, 'theta':theta, 'psi':psi, 'sx':sx, 'sy':sy, 'mirror': mn})

mask = model_circle(34,nx,nx)
# calculate total average using current alignment parameters
tave,tvar = ave_var_series_g(data,kb)
dropImage(tave,"a1.spi")
dropImage(tvar,"a2.spi")
a0 = tave.cmp("dot", tave, {"negative":0,"mask":mask})
print  "initial ",a0
# do the alignment
# IMAGES ARE SQUARES!
#  center is in SPIDER convention
cnx = int(nx/2)+1
cny = cnx

first_ring = 1
last_ring = 35
rstep = 1
xrng  = 1
yrng  = 1
step  = 1
mode = "F"
#precalculate rings
numr=Numrinit(first_ring,last_ring,rstep,mode)
wr=ringwe(numr,mode)
for iter in xrange(20):
	print " ITERATION #",iter
	again = False
	for im in xrange(nima):
		# subtract current image from the average
		psi = data[im].get_attr('psi')
		sx =  data[im].get_attr('sx')
		sy =  data[im].get_attr('sy')
		mirror =  data[im].get_attr('mirror')
		temp = rtshgkb(data[im], psi, sx, sy, kb)
		if  mirror: temp = temp.process("mirror",{"axis":'x'})
		#  Subtract current image from the average
		refim = tave - temp/nima
		cimage=Util.Polar2Dmi(prepg(refim,kb),cnx,cny,numr,mode,kb)
		Util.Frngs(cimage, numr)
		Applyws(cimage, numr, wr)
		# align current image to the reference minus this image
		[angt,sxst,syst,mirrort,peakt]=ormqip(prepg(temp,kb),cimage,first_ring,last_ring,rstep,xrng,yrng,step,mode,numr,kb,cnx,cny,nx)
		# combine parameters and set them to the header
		[psin,sxn,syn,mn]=combine_params2(psi,sx,sy,mirror,angt,sxst,syst,mirrort)
		# apply params to the image
		temp = rtshgkb(data[im], psin, sxn, syn, kb)
		if  mn: temp = temp.process("mirror",{"axis":'x'})
		#check whether the criterion actually increased
		# add current image to the average
		temp = refim + temp/nima
		# calculate the criterion
		a1 = temp.cmp("dot", temp, {"negative":0,"mask":mask})
		if(a1>a0):
			print  im,"  ",a1,"  ",mirror,"  ",mirrort,"  ",psi,"  ",angt,"  ",sxst,"  ",syst
			# replace the average by the improved average and set the new parameters to the image, otherwise, do nothing
			tave = temp.copy()
			data[im].set_attr_dict({'psi':psin, 'sx':sxn, 'sy':syn,'mirror': mn})
			a0 = a1
			again = True
	if(again):
		# calculate total average using current alignment parameters
		av1,av2 = ave_oe_series_g(data,kb)
		frsc = fsc(av1,av2,1.0,"drm%02d"%iter)
		# write the current average
		dropImage(tave,"aqm%04d.spi"%iter)
	else:
		break

for im in xrange(nima):
	temp.read_image(stack_data,im)
	psi = data[im].get_attr('psi')
	sx =  data[im].get_attr('sx')
	sy =  data[im].get_attr('sy')
	mn =  data[im].get_attr('mirror')
	temp.set_attr_dict({'psi':psi, 'sx':sx, 'sy':sy,'mirror': mn})
	temp.write_image(stack_data,im)
