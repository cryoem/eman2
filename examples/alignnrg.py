#!/usr/bin/env python

# This computes the local energy surface for 2-D alignment of a particle vs a class-average


from EMAN2 import *
from sys import argv,exit
from numpy import *

if len(argv)<2 : 
	print "alignnrg <refine#> <it#> <ptcl#> <cmp>"
	exit(1)
argv[1]=int(argv[1])
argv[2]=int(argv[2])
argv[3]=int(argv[3])
cmp=parsemodopt(argv[4])

# get the name of the input particles
db=db_open_dict("bdb:refine_%02d#register"%argv[1],True)
ptcl=db["cmd_dict"]["input"]

# class for the particle & orientation parameters
clmx=EMData.read_images("bdb:refine_%02d#cls_result_%02d"%(argv[1],argv[2]))
projn=int(clmx[0][0,argv[2]])
dx=clmx[2][0,argv[2]]
dy=clmx[3][0,argv[2]]
da=clmx[4][0,argv[2]]
df=int(clmx[5][0,argv[2]])

im1=EMData(ptcl,argv[3])
im2=EMData("bdb:refine_%02d#projections_%02d"%(argv[1],argv[2]),projn)

out=EMData(50,50,40)
x=0
for tx in arange(dx-5,dx+5,0.2):
	y=0
	for ty in arange(dy-5,dy+5,0.2):
		a=0
		for ta in arange(da-10,da+10,0.5):
			xfm=Transform({"type":"2d","tx":tx,"ty":ty,"alpha":ta,"mirror":df})
			im1a=im1.copy()
			im1a.transform(xfm)
			sim=im2.cmp(cmp[0],im1a,cmp[1])
			out[x,y,a]=sim
			a+=1
		y+=1
	x+=1

#display(out)
out.write_image("out_%d.hdf"%argv[3],0)
