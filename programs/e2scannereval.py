#!/bin/env python

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
from Simplex import Simplex
from bisect import insort

cmp_probe=None
cmp_target=None
tdim=None
pdim=None
c2alt=0
degrad=pi/180.0

def display(img):
	img.write_image("tmploc.mrc")
	os.system("v2 tmploc.mrc")

def main():
	global tdim,pdim
	global cmp_probe,cmp_target
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: %prog [options] input.mrc
	
Designed for use on scanned micrographs. Overlays image with a grid of
2D power spectra."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--box", "-S", type="int", help="size in pixels of the power spectra", default=256)
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input file required")
	
	# read the target and probe
	target=EMData()
	target.read_image(args[0])
	sig=target.get_attr("sigma")
	
	nx=target.get_xsize();
	ny=target.get_ysize();
	
	nbx=nx/int(1.5*options.box)		# number of boxes in x
	sepx=options.box*3/2+(nx%(options.box*3/2))/nbx-1	# separation between boxes
	nby=ny/int(1.5*options.box)
	sepy=options.box*3/2+(ny%int(1.5*options.box))/nby
	
	for x in range(nbx):
		for y in range(nby):
			cl=target.get_clip(Region(x*sepx+(sepx-options.box)/2,y*sepy+(sepy-options.box)/2,options.box,options.box))
			cl.filter("RealToFFT")
			cl*=(float(sig)/float(cl.get_attr("sigma")))
			target.insert_clip(cl,(x*sepx+(sepx-options.box)/2,y*sepy+(sepy-options.box)/2,0))

	target.write_image(args[0][:-3]+"eval.mrc")

if __name__ == "__main__":
    main()
