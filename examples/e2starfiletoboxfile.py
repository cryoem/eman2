#!/usr/bin/env python

import os
import sys
from EMAN2 import *

def main():

	usage = """Program to convert an "_autopick.star" file to .box coordinate format."""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--starfile",type=str,default=None,help="""The _autopick.star file you wish to convert to .box format""",required=True)
	parser.add_argument("--boxsize",help="""Specify the boxsize for each particle.""", required=True, type=int)
	
	(options, args) = parser.parse_args()
	
	bf = options.starfile.replace('_autopick.star','.box')
	bs = options.boxsize
	
	if os.path.isfile(bf):
		print("The file {} already exists. Please move or remove it before running this program.".format(bf))
		sys.exit(1)
	
	logger = E2init(sys.argv)
	
	with open(options.starfile) as sf:
		lns = filter(None,[l.split()[:2] for l in sf if '_' not in l if l != '\n'])
	
	with open(bf,'w+') as bf:
		for ln in lns:
			x = ln[0] - bs / 2
			y = ln[1] - bs / 2
			bf.write("{}\t{}\t{}\t{}\n".format(x,y,bs,bs))
	
	E2end(logger)
	
	return

if __name__ == "__main__":
	main()