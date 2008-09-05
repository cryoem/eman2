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

from EMAN2 import *
from sparx import *

from sys import  *
from string import *

'''
  The purpose of this program is to import into SPARX results of SPIDER refinement
  done using a setup that involves defocus groups.  
  data*** are the files with original particles multiplied by the CTF
  trans*** are the transformation parameters
  
  This program alos changes the pixel size (increases) and windows the data to 128^ size (hardwired)
'''

doc_home = "/mnt/fast2/spahn/ttEF-TU/HIGH/"
prj_home = "/mnt/fast1/connell/EFG/DATA-EFG-GTP/"

i_ctfs = open( doc_home + "ctfs_TFCRF083.FM1" )
ignore = i_ctfs.readline()

pixel = 1.22

scale = 1.26/2.6 # used the original pixel size

Cs = 2.0
voltage = 300.0
amp_contrast = 0.1

total_high_proj = 0

proj_out = "/home/christian/EFTU.hdf"

for ii in xrange(1,222) :

    proj_in = "/mnt/fast1/spahn/ttEF-TU/HIGH/DATA-MuCTF/data-MuCTF%3d.FM1" % ii
    proj_in = replace( proj_in, ' ', '0' )

    fn_angl = doc_home + "FINAL/defgrp%3d/angles135.FM1" % ii
    fn_angl = replace( fn_angl, ' ', '0' )

    fn_tran = doc_home + "FINAL/defgrp%3d/trans135.FM1" % ii
    fn_tran = replace( fn_tran, ' ', '0' )

    fn_high = doc_home + "FINAL/defgrp%3d/select.FM1" % ii
    fn_high = replace( fn_high, ' ', '0' )

    nimage  = EMUtil.get_image_count( proj_in )

    f_trans = open( fn_tran, "r" )
    ignore = f_trans.readline()

    f_angle = open( fn_angl, "r" )
    ignore = f_angle.readline()

    defocus = atof( split( i_ctfs.readline() )[2] )

    print "Converting ", proj_in, ", defocus is ", defocus

    for i in xrange(nimage):
        if( i%60 == 0 ) :
            stdout.write( " %6d  " % i )
            stdout.flush()
        trans_line = f_trans.readline()
        angle_line = f_angle.readline()

        data = EMData()
        data.read_image( proj_in, i )


        # use the following code instead of readSpiderDoc because sometimes in trans doc two numbers
        # are connected together like "8.1753E-02-3.1129", readSpiderDoc can not handle this.
        psi   = atof( angle_line[ 7:19] )
        theta = atof( angle_line[19:31] )
        phi   = atof( angle_line[31:43] )

        alpha = atof( trans_line[ 7:19] )
        sx    = atof( trans_line[19:31] )
        sy    = atof( trans_line[31:43] )
        mir   = atof( trans_line[43:55] )
        active = 1

        assb, sxnb, synb, sc = compose_transform2(0.0, sx, sy, 1,  -alpha, 0.,0.,1)
        if( mir > 0.5 ) :
            phi = (540.0 + phi)%360.0
            theta = 180.0 - theta
            psi = (540.0 - psi + assb)%360.0
        else :
            psi = (psi + assb)%360.0

	deci = Util.window(resample(fshift(data,sxnb,synb),scale),128,128,1,0,0,0)
        #the following two line should be used in the case of decimation
        deci.set_attr_dict({'phi':phi, 'theta':theta, 'psi':psi, 's2x':0.0, 's2y':0.0, 'mirror':0.0})
        deci.set_attr_dict({'active':1, 'ctf_applied':1})
        deci.set_attr_dict({'defocus':defocus, 'amp_contrast':amp_contrast, 'voltage':voltage, 'Cs':Cs})
        deci.set_attr_dict( {'Pixel_size':pixel/scale, 'B_factor':0.0} )
        deci.set_attr_dict({'B_factor':0.0})

        deci.write_image( proj_out, total_high_proj )
        total_high_proj += 1

        stdout.write( "." )
        stdout.flush( )

        if( i%60 == 59 ) :
            stdout.write( "\n" )

    stdout.write( "\n" )


print "Totally ", total_high_proj, "  high wrote"

