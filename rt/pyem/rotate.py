#!/bin/env python

from EMAN2 import *
import os
import math

def rotate_3d():
	img = os.environ["HOME"] + "/images/3d.mrc"

	a = EMData()
	a.read_image(img)
	b=a.copy(0)
	b.rotate(0,0,math.pi/4)
	b.write_image("3d1.mrc")
	b=a.copy(0)
	b.rotate(0,0,math.pi/2)
	b.write_image("3d2.mrc")

def rotate_2d():
	img = os.environ["HOME"] + "/images/lattice.mrc"
	
	a = EMData()
	a.read_image(img)
	b=a.copy(0)
	b.rotate(0,0,math.pi/4)
	b.write_image("2d1.mrc")
	b=a.copy(0)
	b.rotate(0,0,math.pi/2)
	b.write_image("2d2.mrc")


rotate_2d()
rotate_3d()
