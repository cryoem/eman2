#!/bin/env python

from EMAN2 import *
import math
import sys
import testlib

def rotate_3d():
	img = TestUtil.get_debug_image("3d.mrc")

	a = EMData()
	a.read_image(img)
	b=a.copy(0)
	b.rotate(0,0,math.pi/4)
	testlib.check_emdata(b, sys.argv[0])

	b=a.copy(0)
	b.rotate(0,0,math.pi/2)
	testlib.check_emdata(b, sys.argv[0])
	
def rotate_2d():
	img = TestUtil.get_debug_image("lattice.mrc")
	
	a = EMData()
	a.read_image(img)
	b=a.copy(0)
	b.rotate(0,0,math.pi/4)
	testlib.check_emdata(b, sys.argv[0])
	
	b=a.copy(0)
	b.rotate(0,0,math.pi/2)
	testlib.check_emdata(b, sys.argv[0])


rotate_2d()
rotate_3d()
