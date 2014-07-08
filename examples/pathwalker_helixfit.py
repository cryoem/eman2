#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Fit helixes using pathwalker results

import EMAN2
from EMAN2 import *
import sys
import numpy.linalg as LA
import random
import math
import numpy as np

def read_fixed(edgefile):
	# Edge file format:
	# 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37
	# 39 40 41 42 43
	# ...
	fixededges = []
	if not edgefile:
		return fixededges
	f = open(edgefile)
	fragments = [i.strip() for i in f.readlines()]
	f.close()
	
	for fragment in fragments:
		fragment = map(int, fragment.split())
		for i in range(len(fragment)-1):
			fixededges.append((fragment[i], fragment[i+1]))
	return fixededges

def read_pdb(filename):
    
    atomnumber=np.array([])
    pdbfile = open(filename, "r")
    lines = pdbfile.readlines()
    pdbfile.close()

    count = 0
    for line in (i for i in lines if i.startswith("ATOM  ")):
		atomnumber=np.append(atomnumber,int(line[22:27]))
		
    return atomnumber

def main():
	
	usage = """ Usage...
	"""
	parser = EMAN2.EMArgumentParser(usage=usage,version=EMAN2.EMANVERSION)
	parser.add_argument("--output", type=str,help="Output pdb file")
	parser.add_argument("--mrcin", type=str,help="mrc file for input",default=None)
	parser.add_argument("--pdbin", type=str,help="pdb file for input",default=None)
	parser.add_argument("--lenthr", type=int,help="length threshold of helixes",default=13)
	parser.add_argument("--denthr", type=float,help="density threshold of helixes",default=4)
	parser.add_argument("--mapwohelix", type=str,help="Write a map without helix density",default=None)
	parser.add_argument("--dirs", type=int,help="Counting from one direction?",default=0)
	parser.add_argument("--edgefile", type=str,help="Existing helixes file",default=None)
	(options, args) = parser.parse_args()
	
	eg=[]
	if options.edgefile<>None:
		edge=read_fixed(options.edgefile)
		eg.append(edge[0][0])
		for i in range(1,len(edge)):
			if edge[i][0]<>edge[i-1][1]:
				eg.append(edge[i-1][1])
				eg.append(edge[i][0])
		eg.append(edge[len(edge)-1][1])
		atomnumber=read_pdb(options.pdbin)
		print eg
		for i in range(len(eg)):
			for j in range(len(atomnumber)):
				if atomnumber[j]==eg[i]:
					eg[i]=j
					break
		print eg
		#exit()
	
	mrc=EMData(options.mrcin)
	atoms=PointArray()
	atoms.read_from_pdb(options.pdbin)
	#atoms.reverse_chain()
	#mrc.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.5})

	hlx=atoms.fit_helix(mrc,options.lenthr,options.denthr,eg,options.dirs)
	
	for i in range(len(hlx)/8):
		print hlx[i*8],hlx[i*8+1]
	
	
	atoms.save_pdb_with_helix(options.output,hlx)
	#atoms.save_to_pdb(options.output)
	if options.mapwohelix<>None:
		atoms.remove_helix_from_map(mrc,hlx)
		mrc.write_image(options.mapwohelix)


if __name__ == '__main__':
	main()
	