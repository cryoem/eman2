#!/bin/env python
from EMAN2  import *
from sparx  import *
from random import seed,random

from sys import argv

def peak_range( nx, defocus, Cs, voltage, pixel ):
    ctf = ctf_1d( nx, pixel, defocus, voltage, Cs )
    
    for i in xrange( 1, len(ctf)-1 ):
        prev = ctf[i-1]
	curt = ctf[i]
	next = ctf[i+1]

	if curt > prev and curt > next:
	    freq = float(i)/nx
	    print 'peak at ', i, freq
            return [freq-0.03, freq+0.02]

    assert false


stack_in = argv[1]
vol = getImage( argv[2] )
stack_ot = argv[3]

volft, kb = prep_vol( vol )

nimage = EMUtil.get_image_count( stack_in )


scales = []

for i in xrange(nimage):

    exp_prj = get_im( stack_in, i )

    phi = exp_prj.get_attr( 'phi' )
    theta = exp_prj.get_attr( 'theta' )
    psi = exp_prj.get_attr( 'psi' )
    s2x = exp_prj.get_attr( 's2x' )
    s2y = exp_prj.get_attr( 's2y' )

    assert s2x==0.0 and s2y==0.0


    ref_prj = prgs( volft, kb, [phi,theta,psi,0.0,0.0] )

    ref_prj = filt_btwo( ref_prj, 0.01,0.1,0.2)
    

    defocus = exp_prj.get_attr( 'defocus' )
    Cs = exp_prj.get_attr( 'Cs' )
    voltage = exp_prj.get_attr( 'voltage' )
    pixel = exp_prj.get_attr( 'Pixel_size' )

    ref_ctfprj = filt_ctf( ref_prj, defocus, Cs, voltage, pixel )

    nx = ref_prj.get_xsize()

    range = peak_range(  nx, defocus, Cs, voltage, pixel )

    d,a,b = im_diff(filt_tophatb(ref_ctfprj, range[0], range[1], False), filt_tophatb(exp_prj, range[0], range[1], False), 28)

    scales.append( 1/a )


avg_scale = sum(scales) / len(scales)

for i in xrange( len(scales) ):
    print "%10.5f" % ( avg_scale/scales[i] )



