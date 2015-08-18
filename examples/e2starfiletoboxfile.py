#!/usr/bin/env python

import os
import sys
from EMAN2 import *

def main():

	usage = """Program to convert an "_autopick.star" file to .box coordinate format."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--starfile",help="""The _autopick.star file you wish to convert to .box format""",required=True,type=str)
	parser.add_argument("--boxsize",help="""Specify the boxsize for each particle.""",required=True,type=int)
	parser.add_argument("--noclobber",help="""If specified, this program will not overwrite existing box files.""",action="store_true",default=False)

	(options, args) = parser.parse_args()

	boxfile = options.starfile.replace('_autopick.star','.box')
	bs = int(options.boxsize / 2)

	if options.noclobber and os.path.isfile(boxfile):
		print("The file {} already exists. It will not be overwritten.".format(boxfile))
		sys.exit(1)

	logger = E2init(sys.argv)

	with open(options.starfile) as sf:
		lines = sf.readlines()

	hdr = [l.split()[0] for l in lines if '_' in l][2:]
	xind = [i for i,x in enumerate(hdr) if x == '_rlnCoordinateX'][0]
	yind = [i for i,x in enumerate(hdr) if x == '_rlnCoordinateY'][0]
	data = filter(None,[l.split() for l in lines if '_' not in l if l != '\n'])
	coords = [[l[xind],l[yind]] for l in data]

	with open(boxfile,'w+') as bf:
		for coord in coords:
			x = int(float(coord[0]))
			y = int(float(coord[1]))
			bf.write("{}\t{}\t{}\t{}\n".format(x-bs/2,y-bs/2,bs,bs))

	E2end(logger)

	return

if __name__ == "__main__":
	main()
