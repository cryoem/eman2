#!/bin/env python


from EMAN2 import *

def run():
	imagename = TestUtil.get_debug_image("monomer.mrc")
	
	a=EMData()
	a.read_image(imagename)
	b=a.get_rotated_clip(Transform([24,24,24], Transform.EulerType.EMAN,0,0,0),
						 [32,32,32],1.0)
	b.write_image("imageio_caching_2.mrc")
	b.write_image("imageio_caching_2.mrc")
	
def write_twice():
	e = EMData()
	e.set_size(100,100,1)
	e.write_image("imageio_caching_1.mrc")
	e.write_image("imageio_caching_1.mrc")

write_twice()
run()
