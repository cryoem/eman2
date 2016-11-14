#!/usr/bin/env python

from sys import argv
from EMAN2 import *
from math import *
from numpy import *

img=EMData(argv[1],int(argv[2]))

resultx=[]
resulty=[]
for i in range(2,30):
	csum=0
	n=0
	for ang in arange(360.0/i,360.0,360.0/i):
		imc=img.copy()
		imc.rotate(ang,0,0)
#		display((imc,img))
		csum+=imc.cmp("ccc",img,{"negative":0})
		n+=1
		
	resultx.append(i)
	resulty.append(csum/n)
	print i,csum/n


# This is a plot of peak values vs peak location
plot((resultx,resulty))


