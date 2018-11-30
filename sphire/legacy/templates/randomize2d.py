#!/bin/env python
from __future__ import print_function

from builtins import range
pass#IMPORTIMPORTIMPORT from EMAN2  import *
pass#IMPORTIMPORTIMPORT from sparx  import *

pass#IMPORTIMPORTIMPORT from sys import exit

# apply random alignment parameters to 2D images
# invert the transformation and store parameters in headers
#  so transform2d will produce properly aligned image series.

import EMAN2_cppwrap
import alignment
import fundamentals
import random
import utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import alignment
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import random
pass#IMPORTIMPORTIMPORT import sparx
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import utilities
pass#IMPORTIMPORTIMPORT from random import random,randint
stack_data = "ang12_2.hdf"
outf = "rang12_2.hdf"
nima = EMAN2_cppwrap.EMUtil.get_image_count(stack_data)
attributes = ['alpha', 'sx', 'sy', 'mirror']
data = EMAN2_cppwrap.EMData()
data.read_image(stack_data, 0, True)
im = data.get_xsize()
kb = alignment.kbt(im)

for im in range(nima):
	data = EMAN2_cppwrap.EMData()
	data.read_image(stack_data, im)
	sx = (random()-0.5)*10.0
	sy = (random()-0.5)*10.0
	alpha = random()*360.0
	mir = randint(0,1)
	data = fundamentals.rot_shift2D(data,alpha,sx,sy,interpolation_method="gridding")
	#  invert the transformation.
	alphah, sxh, syh, sc = utilities.compose_transform2(0.0,-sx, -sy, 1.0 ,-alpha,0.,0.,1.)
	if(mir):
		data = fundamentals.mirror(data)
		alphah, sxh, syh, sc = utilities.combine_params2(0.,0.,0.,1.,alphah,sxh, syh,0.)
	utilities.set_arb_params(data, [alphah, sxh, syh, mir], attributes)
	data.write_image(outf, im)
