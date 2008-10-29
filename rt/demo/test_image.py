#!/usr/bin/env python

from EMAN2 import *

e = EMData()
e.set_size(800, 800)

e.process_inplace('testimage.circlesphere', {'radius':200, 'fill':1})
e.write_image('circle.hdf')

e.process_inplace('testimage.squarecube', {'edge_length':200, 'fill':1})
e.write_image('square.hdf')

e.process_inplace('testimage.sinewave.circular', {'wavelength':150})
e.write_image('sinewave_circular.hdf')

e.process_inplace('testimage.sinewave', {'wavelength':150})
e.write_image('sinewave.hdf')

e.process_inplace('testimage.scurve')
e.write_image('scurve.hdf')

e.process_inplace('testimage.gaussian', {'sigma':300})
e.write_image('gaussian_blob.hdf')

e.process_inplace('testimage.puregaussian', {'sigma':300})
e.write_image('gaussian_blob2.hdf')

e.process_inplace('testimage.noise.gauss')
e.write_image('gaussian_noise.hdf')

e.process_inplace('testimage.noise.uniform.rand')
e.write_image('white_noise.hdf')

e.process_inplace('testimage.noise.uniform.rand')
e.write_image('white_noise.hdf')
