#!/usr/bin/env python

import os
import sys
from EMAN2 import *

def main():

	usage = """Program to convert an "_autopick.star" file to .box coordinate format."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--starfile",help="""The _autopick.star file you wish to convert to .box format""",required=True,type=str)
	parser.add_argument("--boxsize",help="""Specify the boxsize for each particle.""",required=True,type=int)
	
	(options, args) = parser.parse_args()
	
	boxfile = options.starfile.replace('_autopick.star','.box')
	bs = int(options.boxsize / 2)

	#if os.path.isfile(boxfile):
	#	print("The file {} already exists. Please move or remove it before running this program.".format(boxfile))
	#	sys.exit(1)
	
	logger = E2init(sys.argv)
	
	with open(options.starfile) as sf:
		lns = filter(None,[l.split()[:2] for l in sf if '_' not in l if l != '\n'])
		lns = [l for l in lns]
	
	with open(boxfile,'w+') as bf:
		for ln in lns:
			bf.write("{}\t{}\t{}\t{}\n".format(ln[0]-bs,ln[1]-bs,bs,bs))
	
	E2end(logger)
	
	return

if __name__ == "__main__":
	main()