#! /usr/bin/env python
from sparx import *
from EMAN2 import *
e=test_image()
#check rotation
scale=.7

sx = -2.1
sy = -3.7

#sx=0
#sy=0
for i in xrange(36):
	img   = e.copy()
	alpha = i*10
	out1  = rot_shift2D(img,alpha,sx,sy,"linear",scale)
	out2  = rot_shift2D(img,alpha,sx,sy,"quadratic",scale)
	out3  = rot_shift2D(img,alpha,sx,sy,"gridding",scale)
	print  ccc(out1,out2),ccc(out1,out3),ccc(out3,out2)
	out1.write_image("linear.spi",i)
	out2.write_image("quadratic.spi",i)
	out3.write_image("gridding.spi",i)
