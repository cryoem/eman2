#!/bin/env python

from EMAN2 import *

def test1(img1, img2):
	a=EMData()
	b=EMData()
	a.read_image(img1)
	b.read_image(img2)
	c=a.calc_ccf(b,1)
	c.write_image("test1_ccf2.mrc")


def test2(img1):
	a=EMData()
	a.read_image(img1)
	b=a.copy(0)

	c=a.calc_ccf(b,1)
	c.write_image("test21_ccf2.mrc")
	
	c=a.calc_ccf(b,1,b)
	c.write_image("test22_ccf2.mrc")

img1 = TestUtil.get_debug_image("3d86_1.mrc")
img2 = TestUtil.get_debug_image("3d86_2.mrc")

test1(img1, img2)
test2(img1)

