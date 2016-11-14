#!/usr/bin/env python

from EMAN2 import *
from sys import argv,stdout,exit

if len(argv)<4 : 
	print """usage: ringwedgefill <input stack> <output stack> <fill image>
This program will identify and fill in the missing wedge in a stack of subtomograms (with a missing wedge) which have been aligned to an average (with no missing wedge). While in a sense this creates 'Frankenstien particles' it is assumed that having values from an average of many particles is better than having zero in the same areas, at least for certain purposes.

You must provide an average volume from which to draw the missing values and it is critical that this volume be an average of the particles being corrected. The particles should also have been properly normalized prior to averaging to achieve the desired effect.

Specifically, this process is designed to make it possible to run PCA or other classification processes on the particles. Any missing values will all be at the same point near the center of the cloud defined by the entire population, rather than having zero values distorting the cloud. Once particles have been classified, the originals wihout filled wedge can be used for averaging. 
"""
	exit(1)

fill=EMData(argv[3],0)
fillf=fill.do_fft()

n=EMUtil.get_image_count(argv[1])
for i in xrange(n):
	im=EMData(argv[1],i)
	imf=im.do_fft()
	imf.process_inplace("mask.wedgefill",{"fillsource":fillf})
#	imf.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.05})
#	imf.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
	im=imf.do_ift()
	im.process_inplace("normalize")
	im.write_image(argv[2],i)

	if i%10==0:
		print "  %d/%d\r"%(i+1,n),
		sys.stdout.flush()

print "\ndone\n"

