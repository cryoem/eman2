#!/bin/env python
from EMAN2 import *
from sparx  import *
from time import time
mode = "F"
a = model_circle(30,300,300)
#info(a)
#aft,kb = prepi(a)
#cimage=aft.rot_scale_conv(0.0, 0.0, 0.0, kb)
#info(cimage)

for i in xrange(100000):
	f = fft(a)
	b = fft(a)


refm = []
img = []

for i in xrange(20):
	refm.append(a.copy())

for i in xrange(2000):
	img.append(a.copy())

print  "START      ",ttime()
start = time()
apmqs(img, refm, 1, 5, 2, 2, 1, 1, mode)
print  "apmq   ",time() - start
