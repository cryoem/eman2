#!/usr/bin/env python

# clstoclsmx.py - Steve Ludtke, 7/28/08
# This program will convert a set of cls*lst files from EMAN1 into an EMAN2
# style classification matrix (for use with programs like e2classaverage.py

import os
from EMAN2 import *

files=[i for i in os.listdir(".") if i[:3]=="cls" and i[-4:]==".lst"]
files.sort()

mx={}
pm=0
for f in files:
	cls=int(f[3:-4])
	for l in open(f,"r"):
		if l[0]=="#" : continue
		ptcl=int(l.split()[0])
		mx[ptcl]=cls
		if ptcl>pm : pm=ptcl

a=EMData(1,pm+1)
for i in mx:
	a[0,i]=mx[i]

a.write_image("classmx.hdf",0)	# y=ptcl #, value = class #
a.to_one()
a.write_image("classmx.hdf",1)  # weight, one class/ptcl so weights are 1
a.to_zero()
a.write_image("classmx.hdf",2)  # dx,dy,da all zero
a.write_image("classmx.hdf",3)  
a.write_image("classmx.hdf",4) 
