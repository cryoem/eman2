#!/bin/env python

import os
from EMAN2 import *

def run():
	imagename = os.environ['HOME'] + "/images/monomer.mrc"
	
	a=EMData()
	a.read_image(imagename)
	b=a.get_rotated_clip([24,24,24],Rotation(0,0,0,Rotation.Type.EMAN),[32,32,32],1.0)
	b.write_image("z.mrc")
	b.write_image("z.mrc")
	
def write_twice():
	e = EMData()
	e.set_size(100,100,1)
	e.write_image("test1.mrc")
	e.write_image("test1.mrc")

write_twice()
run()
