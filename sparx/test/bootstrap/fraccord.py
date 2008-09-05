#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from mpi import *
from string import replace


       
import sys

prj_stack = sys.argv[1]
vol = getImage( sys.argv[2] )
eigvol = getImage( sys.argv[3] )
nimage = EMUtil.get_image_count( prj_stack )

volft, kb = prep_vol( vol )
eigvolft, eigkb = prep_vol( eigvol )

m28 = model_circle( 28, 75, 75 )

for i in xrange( nimage ) :
        exp_prj = get_im( prj_stack, i )

        nx = exp_prj.get_xsize()
        phi = exp_prj.get_attr( 'phi' )
        theta = exp_prj.get_attr( 'theta' )
        psi = exp_prj.get_attr( 'psi' )
        s2x = exp_prj.get_attr( 's2x' )
        s2y = exp_prj.get_attr( 's2y' )
        defocus = exp_prj.get_attr( 'defocus' )
        wgh = exp_prj.get_attr( 'amp_contrast' )
        Cs = exp_prj.get_attr( 'Cs' )
        voltage = exp_prj.get_attr( 'voltage' )
        pixel = exp_prj.get_attr( 'Pixel_size' )


        shift_params = {"filter_type" : Processor.fourier_filter_types.SHIFT,
			"x_shift" : s2x, "y_shift" : s2y, "z_shift" : 0.0}
        exp_prj =  Processor.EMFourierFilter(exp_prj, shift_params)


        ref_prj = prgs( volft, kb, [phi, theta, psi, 0.0, 0.0] )
        ref_eigprj = prgs( eigvolft, eigkb, [phi, theta, psi, 0.0, 0.0] )
        

        ref_ctfprj = filt_ctf( ref_prj, defocus, Cs, voltage, pixel, wgh )
        ref_ctfeigprj = filt_ctf( ref_eigprj, defocus, Cs, voltage, pixel, wgh )


        diff,a,b = im_diff( ref_ctfprj, exp_prj, m28)

        print 'prjid, defocus, frac cord: %6d %10.3f %10.3e ' % ( i, defocus, diff.cmp("dot", ref_ctfeigprj, {"negative":0, 'mask':m28})  )


