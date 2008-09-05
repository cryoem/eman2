#!/bin/env python

from EMAN2  import *
from sparx  import *

from random import seed,gauss,uniform
from math import sqrt
seed()
j=0
aiq=0.
viq=0.
miq = 0.0
for k in xrange(1000):
	sums = 0.
	sums2 = 0.
	n = 100
	am = -1.0e23
	qs = uniform(1.,1000.)
	for i in xrange(n):
		x = gauss(0.0,qs)
		sums += x
		sums2 += x**2
		if(x>=am): am = x
	a=sqrt((sums2-sums**2/n)/(n-1.0))
	iq = (am-sums/n)/a
	print  iq
	aiq += iq
	viq += iq*iq
	if(iq >= miq):
		miq = iq

print  aiq/1000. ,sqrt((viq-aiq**2/1000.)/(1000.-1.0)), miq
