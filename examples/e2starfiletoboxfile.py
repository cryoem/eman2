#!/usr/bin/env python

import os
import sys
from EMAN2 import *
from EMAN2star import *

from IPython import embed

def main():

	usage = """Program to convert an ".star" file to .box coordinate format."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--input",help="""The .star file you wish to convert to .box format""",type=str)
	parser.add_argument("--output",help="""The name of the .box file to be written""",type=str)
	parser.add_argument("--boxsize",help="""Specify the boxsize for each particle.""",required=True,type=int)
	(options, args) = parser.parse_args()

	if len(args) != 2: exit(-1)
	
	if not options.output: options.output = args[1]
	if not options.input: options.input = args[0]
	
	if os.path.exists(args[0]): starf = StarFile(args[0])
	else: exit(-1)
	
	if options.output.split('.')[-1] != "box":
		f = options.output.split('.')
		if len(f) > 2: bn = ".".join(f[:-2])
		else: bn = f[-2]
		print("Writing to {base}.box rather than {base}.{ext}".format(base=bn,ext=f[-1]))
		options.output = bn+".box"
	boxf = options.output
	
	bs = int(options.boxsize / 2)
	
	logger = E2init(sys.argv)
	
	hdr = starf.keys()
	
	print(hdr)
	
	mk = [i for i in hdr if "Micrograph" in i]
	xk = [i for i in hdr if "X" in i]
	yk = [i for i in hdr if "Y" in i]
	
	# resolve correct x and y dictionary keys
	
	# case 1: only one X key
	# case 2: multiple xk
	# case 3: no xk
	
	# case 1: only one Y key
	# case 2: multiple yk
	# case 3: no yk
	
	embed()
	exit(-1)
	
	# case 1: one file to one micrograph
	if len(mks) < 1:

		coords =
		coords = [[l[xk],l[yk]] for l in data]
		with open(boxf,'w+') as bf:
			for coord in coords:
				x = int(float(coord[0]))
				y = int(float(coord[1]))
				bf.write("{}\t{}\t{}\t{}\n".format(x-bs/2,y-bs/2,bs,bs))
	
	# case 2: multiple micrographs in file, one header name
	elif len(mks) == 1:
		mgs = list(set(starf[mghdr[0]))
		# make box file for each micrograph
		# add each particle to box file
			# store x,y coords
		# write box files with box size. format is "{}\t{}\t{}\t{}".format(x-bs/2,y-bs/2,bs/2,bs/2)
			
	# case 3: multiple micrograph data in header
	else:
		print("Sorry, I'm not ready to handle this yet")
		exit(-1)
	

	E2end(logger)

	return

if __name__ == "__main__":
	main()
