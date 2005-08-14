#!/bin/env python
# sortsimilar.py  08/09/2005  Steven Ludtke
# This will sort projection/image pairs based on their mutual similarity
# more similar pairs will be placed in the output file first

from EMAN2 import *
import sys

cmps=[]

n=EMUtil.get_image_count(sys.argv[1])

# read pairs of images and calculate similarity
# store in a list containing similarity and image number
for i in range(0,n,2):
	l=EMData.read_images(sys.argv[1],(i,i+1))
#	l[1].process("eman1.normalize.edgemean")
#	l[1].process("eman1.filter.lowpass.gaussian",{"lowpass":.08})
#	cmps.append((l[0].cmp("phase",l[1],{})+l[0].cmp("optvariance",l[1],{"radweight":1})/l[0].get_xsize(),i))
	cmps.append((l[0].cmp("optvariance",l[1],{"matchamp":1})/l[0].get_xsize(),i))


#	cmps.append((l[0].cmp("dot",l[1],{"normalize":1}),
#		l[0].cmp("frc",l[1],{}),
#		l[0].cmp("optvariance",l[1],{}),
#		l[0].cmp("phase",l[1],{}),
#		l[0].cmp("quadmindot",l[1],{"normalize":1}),
#		l[0].cmp("variance",l[1],{})))

# sort the list in ascending order. Smaller similarity values are better
cmps.sort()

# write the images to the output file
for i in cmps:
	l=EMData.read_images(sys.argv[1],(i[1],i[1]+1))
	l[0].write_image(sys.argv[2],-1)
	l[1].write_image(sys.argv[2],-1)

