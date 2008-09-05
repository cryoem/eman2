#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from mpi import *
from string import replace, atoi


       
from sys import argv, exit


if len(argv) != 6:
    print "prj_eigvol.py  prj_stack vol_stack eigvol_stack neigvol output_stack"
    exit(-1)


prj_stack = argv[1]
vol = getImage( argv[2] )
feigvol = argv[3]
output  = argv[5]

eigvols = []
neigvol = atoi( argv[4] ) #EMUtil.get_image_count(feigvol)
for i in xrange(neigvol):
    eigvols.append( get_im(feigvol, i) )

volft, kb = prep_vol( vol )

eigvolfts = []
for e in eigvols:
    eigvolfts.append( prep_vol(e) )

m28 = model_circle( 28, 75, 75 )

nimage = EMUtil.get_image_count( prj_stack )

for i in xrange( 103800, nimage ) :
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
        ref_ctfprj = filt_ctf( ref_prj, defocus, Cs, voltage, pixel, wgh )

        diff,a,b = im_diff( ref_ctfprj, exp_prj, m28)

        img = model_blank( len(eigvols) )

        for j in xrange( len(eigvolfts) ) :
            eft = eigvolfts[j]

            ref_eigprj = prgs( eft[0], eft[1], [phi, theta, psi, 0.0, 0.0] )
            ref_ctfeigprj = filt_ctf( ref_eigprj, defocus, Cs, voltage, pixel, wgh )

            d = diff.cmp( "dot", ref_ctfeigprj, {"negative":0, "mask":m28} )

            img.set_value_at( j, 0, 0, d )

        img.write_image( output, i )
        print i, ' done'

