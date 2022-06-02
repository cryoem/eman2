#!/bin/env python
from builtins import range
from sparx  import *
from EMAN2 import *


a = test_image(0,(1000,1000))
kernel = model_blank(15,15,bckg = 0.21)
print(ttime())
for i in range(100):
	b = rsconvolution(a,kernel)
	#if(i%100 == 0):  print i,ttime()
print(ttime())

a = test_image(0,(2000,2000))
kernel = model_blank(3,3,bckg = 0.33)
print(ttime())
for i in range(100):
	b = rsconvolution(a,kernel)
	#if(i%100 == 0):  print i,ttime()
print(ttime())

a = test_image(0,(128,128))
kernel = model_blank(7,7,bckg = 0.21)
print(ttime())
for i in range(10000):
	b = rsconvolution(a,kernel)
	#if(i%100 == 0):  print i,ttime()
print(ttime())

a = test_image(0,(256,256))
kernel = model_blank(3,3,bckg = 0.33)
print(ttime())
for i in range(10000):
	b = rsconvolution(a,kernel)
	#if(i%100 == 0):  print i,ttime()
print(ttime())
