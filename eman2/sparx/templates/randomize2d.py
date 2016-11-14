#!/bin/env python

from EMAN2  import *
from sparx  import *

from sys import exit

# apply random alignment parameters to 2D images
# invert the transformation and store parameters in headers
#  so transform2d will produce properly aligned image series.

from random import random,randint
stack_data = "ang12_2.hdf"
outf = "rang12_2.hdf"
nima = EMUtil.get_image_count(stack_data)
attributes = ['alpha', 'sx', 'sy', 'mirror']
data = EMData()
data.read_image(stack_data, 0, True)
im = data.get_xsize()
kb = kbt(im)

for im in xrange(nima):
	data = EMData()
	data.read_image(stack_data, im)
	sx = (random()-0.5)*10.0
	sy = (random()-0.5)*10.0
	alpha = random()*360.0
	mir = randint(0,1)
	data = rot_shift2D(data,alpha,sx,sy,interpolation_method="gridding")
	#  invert the transformation.
	alphah, sxh, syh, sc = compose_transform2(0.0,-sx, -sy, 1.0 ,-alpha,0.,0.,1.)
	if(mir):
		data = mirror(data)
		alphah, sxh, syh, sc = combine_params2(0.,0.,0.,1.,alphah,sxh, syh,0.)
	set_arb_params(data, [alphah, sxh, syh, mir], attributes)
	data.write_image(outf, im)
