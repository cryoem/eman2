#!/usr/bin/env python
# Uses EMAN2/NumPy to compute a real-space deconvolution matrix for use in non-Cartesian coordinates

from EMAN2 import *
import numpy as np

def mxprint(mx):
	for y in xrange(81):
		for x in xrange(81):
			print "{:3.0f}".format(float(mx[x][y])),
		print ""

def mxprintsm(mx):
	for y in xrange(9):
		for x in xrange(9):
			print "{:7.2f}".format(float(mx[x,y])),
		print ""

mxl=[]
for x in xrange(-4,5):
	for y in xrange(-4,5):
		gau=EMData(9,9,1)
		gau.to_one()
		gau.process_inplace("mask.gaussian",{"outer_radius":1,"dx":x,"dy":y,"dz":0})		# 1/2 width of Gaussian in pixels
		gau.process_inplace("normalize.unitsum")
		a=gau.numpy().reshape(81).copy()
		mxl.append(a)

# mxl is now the convolution matrix for a 9x9 region (with edge issues we ignore)
mx=np.asarray(mxl)
#mxprint(mx)
imx=np.linalg.inv(mx)
#mxprint(imx)


row=imx[40].reshape(9,9)
emd=from_numpy(row)


gau=EMData(9,9,1)
gau.to_one()
gau.process_inplace("mask.gaussian",{"outer_radius":1.9})		# 1/2 width of Gaussian in pixels
gau.process_inplace("normalize.unitsum")
mxprintsm(gau)

print " "
mxprintsm(row)

#display((emd,gau),True)
