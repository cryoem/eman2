#!/bin/env python
# comparepairs.py  08/09/2005  Steven Ludtke

from EMAN2 import *
import sys
from math import *

cmps=[]

n=EMUtil.get_image_count(sys.argv[1])

# read pairs of images and calculate a variety of
# similarity metrics between the pairs
for i in range(0,n,2):
	l=EMData.read_images(sys.argv[1],(i,i+1))
	l[1].process_inplace("eman1.normalize.edgemean")
	l[1].process_inplace("eman1.filter.lowpass.gaussian",{"lowpass":.1})

	cmps.append((
		l[0].cmp("phase",l[1],{}),
		sqrt(l[0].cmp("optvariance",l[1],{"radweight":1})/l[0].get_xsize()),
		sqrt(l[0].cmp("optvariance",l[1],{})),
		l[0].cmp("dot",l[1],{"normalize":1}),
		l[0].cmp("frc",l[1],{}),
		l[0].cmp("variance",l[1],{})))

# calculate the mean for each measure
sm=[0,0,0,0,0,0]
for i in cmps:
	for j,k in enumerate(i):
		sm[j]+=k

sm=[j/len(cmps) for j in sm]

# now print out each similarity measure for each pair
for i,j in enumerate(cmps):
	print "%d.\t"%i,
	for k,l in enumerate(j):
#		print "%1.4f\t"%-(l-sm[k]),
		print "%1.4f\t"%l,
	print ""
