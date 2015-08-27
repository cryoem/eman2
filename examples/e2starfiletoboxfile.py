#!/usr/bin/env python

import os
import sys
from EMAN2 import *
from EMAN2star import *

from IPython import embed

def main():

	usage = """Program to convert an ".star" file to .box coordinate format."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--input",help="""The .star file you wish to convert to .box format""",type=str,required=True)
	parser.add_argument("--output",help="""The name of the .box file to be written""",type=str)
	parser.add_argument("--boxsize",help="""Specify the boxsize for each particle.""",required=True,type=int)
	# need an argument to append word/phrase or dirname to output
	# need path argument for alternative location
	(options, args) = parser.parse_args()

	if len(args) < 1: exit(-1)
	
	if os.path.exists(args[0]): starf = StarFile(args[0])
	else: exit(-1)
	
	if options.output.split('.')[-1] != "box":
		f = options.output.split('.')
		if len(f) > 2: bn = ".".join(f[:-2])
		else: bn = f[-2]
		print("Writing to {base}.box rather than {base}.{ext}".format(base=bn,ext=f[-1]))
		options.output = bn+".box"
		options.output
	
	bs = int(options.boxsize / 2)
	
	logger = E2init(sys.argv)
	
	hdr = starf.keys()
	
	print(hdr)
	
	mks = [i for i in hdr if "Micrograph" in i]
	xk = [i for i in hdr if "X" in i]
	yk = [i for i in hdr if "Y" in i]
	
	# resolve correct x and y dictionary keys
	
	if len(xk) == 1: # case 1: only one X key
		xk = xk[0]
		if not options.output:
			print("No output file name was specified. Will use the input basename as output.")
			options.output = options.input.split('.')[0] + ".box"
	elif len(xk) > 1: # case 2: multiple xk
		xk = [i for i in xk if "Coordinate" in i][0]
	else: # case 3: no xk
		print("Could not find any keys containing 'X'")
		exit(-1)
	
	if len(yk) == 1: # case 1: only one Y key
		yk = yk[0]
	elif len(yk) > 1: # case 2: multiple yk
		yk = [i for i in yk if "Coordinate" in i][0]
	else: # case 3: no xk
		print("Could not find any keys containing 'Y'")
		exit(-1)
	
	embed()
	
	# case 1: one file to one micrograph
	if len(mks) < 1:
		with open(options.output,'w+') as boxf:
			for x,y in starf[xk],starf[yk]):
				boxf.write("{}\t{}\t{}\t{}\n".format(x-bs/2,y-bs/2,bs,bs))
	# case 2: multiple micrographs in file, one header name
	elif len(mks) == 1:
		mgs = list(set(starf[mks[0]]))
		for mg in mgs:
			boxfile = mg.split('.')[0] + ".box"
			with open(boxfile,'w+') as boxf:
				for x,y in zip(starf[xk],starf[yk]):
					boxf.write("{}\t{}\t{}\t{}\n".format(x-bs/2,y-bs/2,bs,bs))
	# case 3: multiple micrograph data in header
	else:
		print("Sorry, I'm not ready to handle this type of .star file yet")
		exit(-1)
	

	E2end(logger)

	return

if __name__ == "__main__":
	main()
