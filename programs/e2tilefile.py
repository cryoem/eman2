#!/usr/bin/env python

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def main():
	global tdim,pdim
	global cmp_probe,cmp_target
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <tile file>
	
Operates on files containing sets of tiled JPEG images representing larger images. Used for 
interactive web browsing."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--build", type="string", help="Build a new tile file from this image")
	parser.add_option("--tilesize", type="int",default=256, help="Build a new tile file from this image")
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Tile file required")
	try: chains=options.chains
	except: chains=None
		
	if options.build:
		# read the target and probe
		orig=EMData()
		orig.read_image(options.build)
	
		build_tiles(orig,args[0])
		
def build_tiles(img,tilefile,tilesize):
	levels=log(min(img.get_xsize(),img.get_ysize())/tilesize)/log(2.0)
	
	
	
if __name__ == "__main__":
    main()
