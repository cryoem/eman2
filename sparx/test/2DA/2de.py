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

#  AMOEBA

#  BEGIN

stack_data="data2.hdf"
nima = EMUtil.get_image_count(stack_data)
temp = EMData()
temp.read_image(stack_data,0)
# prepare KB interpolants
nx=temp.get_xsize()
kb = kbt(nx)

# read images and prepare them for gridding
data = []
params = []
for im in xrange(nima):
	if(im>0):
		temp = EMData()
		temp.read_image(stack_data,im)
	phi = temp.get_attr('phi')
	theta = temp.get_attr('theta')
	psi = temp.get_attr('psi')
	params.append(psi)
	sx =  temp.get_attr('sx')
	params.append(sx)
	sy =  temp.get_attr('sy')
	params.append(sy)
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
# load stuff for amoeba
stuff = []
stuff.insert(0,kb)
stuff.insert(1,mask)
stuff.insert(2,nima)
#stuff.insert(3,tave)  # current average
#stuff.insert(4,data)  # current image in the gridding format
weights = [1.75]*3  # weights define initial bracketing, so one would have to figure how to set them correctly
for iter in xrange(20):
	print " ITERATION #",iter,a0
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
		stuff.append(refim)  # curent ave-1
		stuff.append(data[im])  # curent image
		# perform amoeba alignment
		params = [psi, sx, sy]
		outparams =  amoeba(params, weights, crit2d, 1.e-4, 1.e-4, 500, stuff)
		del stuff[3]
		del stuff[3]
		# set parameters to the header
		data[im].set_attr_dict({'psi':outparams[0][0], 'sx':outparams[0][1], 'sy':outparams[0][2],'mirror': mirror})
		# update the average
		temp = rtshgkb(data[im], outparams[0][0], outparams[0][1], outparams[0][2], kb)
		if  mirror: temp = temp.process("mirror",{"axis":'x'})
		#  Subtract current image from the average
		tave = refim + temp/nima
		print  im,tave.cmp("dot", tave, {"negative":0,"mask":mask}),params,outparams[0],outparams[2]

	# calculate total average using current alignment parameters
	av1,av2 = ave_oe_series_g(data,kb)
	frsc = fsc(av1,av2,1.0,"dra%02d"%iter)
	tave,tvar = ave_var_series_g(data,kb)
	a0 = tave.cmp("dot", tave, {"negative":0,"mask":mask})
	# write the current average
	dropImage(tave,"aam%04d.spi"%iter)

for im in xrange(nima):
	temp.read_image(stack_data,im)
	psi = data[im].get_attr('psi')
	sx =  data[im].get_attr('sx')
	sy =  data[im].get_attr('sy')
	mn =  data[im].get_attr('mirror')
	temp.set_attr_dict({'psi':psi, 'sx':sx, 'sy':sy,'mirror': mn})
	temp.write_image(stack_data,im)
