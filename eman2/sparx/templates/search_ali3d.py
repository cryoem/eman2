#! /usr/bin/env python

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

from   EMAN2        import *
from   sparx        import *
import sys

if len(sys.argv) != 5:
    print 'search_alivol.py <vol_to_align.hdf> <ref_vol.hdf> <name_vol_aligned.hdf> <radius_mask>'
    sys.exit()

# vars
name_vol = sys.argv[1]
name_ref = sys.argv[2]
name_ali = sys.argv[3]
rad      = int(sys.argv[4])
name_mir = 'tmp_vol_m.hdf'

im = EMData()
im.read_image(name_ref, 0, True)
mask = model_circle(rad, im.get_xsize(), im.get_ysize(), im.get_zsize()) 

# create the mirror of the volume
vol = get_im(name_vol)
vol = mirror(vol)
vol.write_image(name_mir)

# define list of angles ~ 15 angles (phi, theta, psi)
agls      = even_angles(40, 0.0, 179.9, 0.0, 359.9, 'P')
size_agls = len(agls)

cc  = [0]     * size_agls
mir = [False] * size_agls

# for each angles search the alignment of the volumes
for n in xrange(size_agls): # another loop over Psi
    # set the value of phi, theta, and psi in the header of the volume
    im.read_image(name_vol, 0, True)
    im.set_attr_dict({'phi': agls[n][0], 'theta': agls[n][1], 'psi': 0, 's3x': 0, 's3y':0, 's3z': 0, 'scale': 1})
    im.write_image(name_vol, 0, EMUtil.ImageType.IMAGE_HDF, True)

    # search the alignment
    ali_vol_rotate(name_vol, name_ref, 20, rad)

    # read the new value of orientation
    im  = get_im(name_vol)

    # apply the orientation on the volume
    vol_ali = rot_shift3D(im, im.get_attr('phi'), im.get_attr('theta'), im.get_attr('psi'))
    vol_ali.write_image(name_ali)

    # compute the cross correlation between the ref and the new orientation of the vol
    #val_cc = imgstat([name_ali, name_ref], True, '', False, (dic['nx'] / 2) - 2)
    tmp_ali = get_im(name_ali)
    tmp_ref = get_im(name_ref)
    val_cc = ccc(tmp_ali, tmp_ref, mask)

    # to display the result
    print '================================'
    print 'Agls: cross correlation: %3.2f' % val_cc
    print '================================\n'
    
    #---- the same on the mirror struvture ---------------------------
    # set the value of phi, theta, and psi in the header of the volume
    im.read_image(name_mir, 0, True)
    im.set_attr_dict({'phi': agls[n][0], 'theta': agls[n][1], 'psi': 0, 's3x': 0, 's3y':0, 's3z': 0, 'scale': 1})
    im.write_image(name_mir, 0, EMUtil.ImageType.IMAGE_HDF, True)

    # search the alignment
    ali_vol_rotate(name_mir, name_ref, 20, rad)

    # read the new value of orientation
    im  = get_im(name_mir)

    # apply the orientation on the volume
    vol_ali = rot_shift3D(im, im.get_attr('phi'), im.get_attr('theta'), im.get_attr('psi'))
    vol_ali.write_image(name_ali)

    # compute the cross correlation between the ref and the new orientation of the vol
    #n_val_cc = imgstat([name_ali, name_ref], True, '', False, (dic['nx'] / 2) - 2)
    tmp_ali = get_im(name_ali)
    tmp_ref = get_im(name_ref)
    n_val_cc = ccc(tmp_ali, tmp_ref, mask)

    # to display the result
    print '================================'
    print 'Agls mirror: cross correlation: %3.2f' % n_val_cc
    print '================================\n'

    #----- choose the best and store
    if n_val_cc > val_cc:
        cc[n]  = n_val_cc
        mir[n] = True
    else:
        cc[n]  = val_cc
        mir[n] = False


#---- choose the best from the all angles ----
val_max  = -1e10
best_pos = -1
best_mir = False
for n in xrange(size_agls):
    if cc[n] > val_max:
        val_max  = cc[n]
        best_pos = n
        best_mir = mir[n]

#---- assign ----
phi   = agls[best_pos][0]
theta = agls[best_pos][1]
if best_mir: name_vol = name_mir

#---- search the alignment for severals values of psi ----
psi = range(0, 360, 36)
cc  = -1e10

for n in psi:
    # set the value of phi, theta, and psi in the header of the volume
    im.read_image(name_vol, 0, True)
    im.set_attr_dict({'phi': phi, 'theta': theta, 'psi': n, 's3x': 0, 's3y':0, 's3z': 0, 'scale': 1})
    im.write_image(name_vol, 0, EMUtil.ImageType.IMAGE_HDF, True)

    # search the alignment
    ali_vol_rotate(name_vol, name_ref, 20, rad)

    # read the new value of orientation
    im  = get_im(name_vol)

    # apply the orientation on the volume
    vol_ali = rot_shift3D(im, im.get_attr('phi'), im.get_attr('theta'), im.get_attr('psi'))
    vol_ali.write_image(name_ali)

    tmp_ali = get_im(name_ali)
    tmp_ref = get_im(name_ref)
    val_cc = ccc(tmp_ali, tmp_ref, mask)

    # to display the result
    print '================================'
    print 'Psi: cross correlation: %3.2f' % val_cc
    print '================================\n'

    # choose the best
    if val_cc > cc:
        cc       = val_cc
        best_pos = n

#---- assignt the best start orientation on the volume ----
im.read_image(name_vol, 0, True)
im.set_attr_dict({'phi': phi, 'theta': theta, 'psi': best_pos, 's3x': 0, 's3y':0, 's3z': 0, 'scale': 1})
im.write_image(name_vol, 0, EMUtil.ImageType.IMAGE_HDF, True)

# search the alignment
ali_vol_rotate(name_vol, name_ref, 1.0, rad)

# read the new value of orientation
im  = get_im(name_vol)

# apply the orientation on the volume
vol_ali = rot_shift3D(im, im.get_attr('phi'), im.get_attr('theta'), im.get_attr('psi'))
vol_ali.write_image(name_ali)

# print the end infos
print '\n\n'
print 'best cc value:', cc
print '%10.3f\t%10.3f\t%10.3f' % (phi, theta, best_pos)
if best_mir: print 'use mirror structure'


