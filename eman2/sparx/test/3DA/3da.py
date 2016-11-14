#!/bin/env python

#
# Author: Pawel A. Penczek, 04/10/2003 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas
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

from random import randint
import os



print  "start"
vol=getImage("../model001.tcp")
#info(vol)
#vol = vol.FourInterpol(320,320,320)
volft,kb=prep_vol(vol)

#angles=even_angles(15.,0.,180.,0.,359.9,'S')#,'Minus')
#angles=even_angles(1.,20.,25.,0.,359.9,'S')#,'Minus')
angles=even_angles(12.,0.,180.0,0.,359.9,'S')#,'Minus')

print  angles
nangles = len(angles)

stack = "proj.hdf"
stack2 = "proj6.hdf"
if os.path.exists(stack):  os.system('rm  '+stack)
if os.path.exists(stack2):  os.system('rm  '+stack2)
sx=0
sy=0

ast = 720/nangles
apsi=-ast
for i in xrange(nangles):
	apsi += ast
	angles[i][2] = apsi
	print  angles[i]
	proj=prgs(volft, kb, [angles[i][0], angles[i][1], angles[i][2], 0., 0.])
	proj.set_attr_dict({'phi':angles[i][0], 'theta':angles[i][1], 'psi':0, 'alpha':0, 'sx':0, 'sy':0, 'mirror':0})
	proj.write_image(stack2, i)

	sx = randint(-4,4)
	sy = randint(-3,3)
	alpha = randint(0,359)
	ass,sxn,syn,sc = compose_transform2(0, -sx, -sy, 1, -alpha,0,0,1)

	'''
	# without mirror
	ps = rtshg(proj, ass, sxn, syn)
	#  alpha, sx, sy are TRUE alignment params, as required to put projection 'in place'
	print  i, angles[i][0], angles[i][1], angles[i][2], sx, sy
	ast, sxt, syt, sc = compose_transform2(0, sx, sy, 1, -alpha,0,0,1)
	print  i, angles[i][0], angles[i][1], ast,sxt,syt
	#ps = rtshg(ps, alpha, sx, sy)
	# Set all parameters for the new 2D image
	# three angles and two shifts
	ps.set_attr_dict({'phi':angles[i][0], 'theta':angles[i][1], 'psi':ast, 'sx':sxt, 'sy':syt, 'mirror':0.})
	#ps.set_attr_dict({'phi':angles[i][0], 'theta':angles[i][1], 'psi':0, 'sx':0, 'sy':0, 'mirror':0.})
	# CTF parameters, if defocus zero, they are undetermined
	ps.set_attr_dict({'defocus':0.0, 'amp_contrast':0.1, 'voltage':200, 'Cs':2.0, 'pixel':2.2})
	# flags describing the status of the image (1 = true, 0 = false)
	ps.set_attr_dict({'ctf_applied':0})
	ps.write_image(stack, i)
	#  without mirror ends here
	'''
	#  try MIRROR
	ps = proj.process("mirror",{"axis":'x'})
	ps = rtshg(ps, ass, sxn, syn)

	angles[i][0] = (540+angles[i][0])%360.0         # phi
	angles[i][1] = 180-angles[i][1]         # theta
	angles[i][2] = (540-angles[i][2])%360.0         #  psi
	#proj=prgs(volft, kb, [angles[i][0], angles[i][1], angles[i][2], 0., 0.])
	#print  ccc(proj,ps)
	#pt = fshift(proj, sx, sy)
	#ps = rtshg(pt, -angles[i][2], 0.0, 0.0)
	print  i, angles[i][0], angles[i][1], angles[i][2], -sx, -sy
	ast, sxt, syt, sc = compose_transform2(0, sx, sy, 1, -alpha,0,0,1)
	print  i, angles[i][0], angles[i][1], ast,sxt,syt
	#ps = rtshg(pt, alpha, sx, sy)
	
	#ps = rtshg(pt, alpha, sx, sy)
	#ps = ps.process("mirror",{"axis":'x'})
	#ps.set_attr_dict({'phi':angles[i][0], 'theta':angles[i][1], 'psi':angles[i][2], 'sx':0., 'sy':0., 'mirror':0.})
	ps.set_attr_dict({'phi':angles[i][0], 'theta':angles[i][1], 'psi':angles[i][2]+ast, 'alpha':0, 'sx':sxt, 'sy':syt, 'mirror':0.})
	ps.set_attr_dict({'defocus':0.0, 'amp_contrast':0.1, 'voltage':200, 'Cs':2.0, 'pixel':2.2})
	# horatio active_refactoring Jy51i1EwmLD4tWZ9_00000_1
	# ps.set_attr_dict({'active':1, 'ctf_applied':0})
	ps.set_attr_dict({'ctf_applied':0})
	ps.write_image(stack2, i)

from sys import exit
exit()

# b=EMData().read_images("proj.hdf")
#stack = "proj.hdf"
m = model_circle(30,75,75,75)
list_proj=range(len(angles))

#v = recons3d_4nn(stack, list_proj, "c1")
#print ccc(v,vol,m)
# 0.943
#dropImage(v,"vv.spi")


v = recons3d_4nn(stack2, list_proj, "c1")
print ccc(v,vol,m)

#dropImage(v,"vt.spi")
#  0.919
'''

v = recons3d_wbp(stack, list_proj, "exact", 75)

dropImage(v,"aa.spi")
v = recons3d_wbp(stack, list_proj, "general")

#v = recons3d_4nn(stack, list_proj, "c1")

dropImage(v,"bb.spi")

#s = recons3d_sirt(stack, list_proj, 35, 1.0e-2)

#dropImage(s,"sirt.spi")
'''
