#!/usr/bin/env python
# Muyuan Chen 12/2014
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
	parser.add_argument("--mapin", type=str,help="mrc file for input",default=None)
	parser.add_argument("--pdbin", type=str,help="pdb file for input",default=None)
	parser.add_argument("--lenthr", type=int,help="length threshold of helixes",default=13)
	parser.add_argument("--minlen", type=int,help="minimun length helixes",default=10)
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
	
	mrc=EMData(options.mapin)
	atoms=PointArray()
	atoms.read_from_pdb(options.pdbin)
	#atoms.reverse_chain()
	#mrc.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.5})
	allhlx=np.array([0,0,0,0,0,0])
	allhlx=np.vstack((allhlx,allhlx))
	for ii in range(1):
		newhlx=0
		hlx=atoms.fit_helix(mrc,options.lenthr,options.denthr,eg,options.dirs,options.minlen)
		for i in range(len(hlx)/8):
			h=[hlx[i*8+2],hlx[i*8+3],hlx[i*8+4],hlx[i*8+5],hlx[i*8+6],hlx[i*8+7]]
			mindist=min(np.sum((allhlx-h)*(allhlx-h),axis=1))
			print h,mindist
			if mindist>100:
				allhlx=np.vstack((allhlx,h))
				newhlx=1
		if newhlx==0:
			break
				
	allhlx=np.delete(allhlx,[0,1],axis=0)
	print allhlx
	
	hlx=[]
	n=atoms.get_number_points()
	inhlx=-1
	for i in range(1,n-1):
		p1=atoms.get_vector_at(i-1)
		p2=atoms.get_vector_at(i)
		p3=atoms.get_vector_at(i+1)
		l1=p2-p1
		l2=p2-p3
		#print l1.length(),l2.length(),i,
		a=(l1.dot(l2))/(l1.length()*l2.length())
		a=min(a,1)
		a=max(a,-1)
		#print a,
		ang=acos(a)
		#print i,ang
		if inhlx<0:
			if abs(ang-1.585)<.01:
				inhlx=i
				inp=p2
		else:
			if abs(ang-1.585)>.01:
				if i-inhlx>5:
					hlx.extend([inhlx,i,inp[0],inp[1],inp[2],p2[0],p2[1],p2[2]])
				inhlx=-1
		
	#for h in allhlx:
		#print h
		#p1=Vec3f(h[0],h[1],h[2])
		#p2=Vec3f(h[3],h[4],h[5])
		#n1=0
		#n2=0
		#minl1=10000
		#minl2=10000
		#for i in range(atoms.get_number_points()):
			#p=atoms.get_vector_at(i)
			#dp=p-p1
			#l=dp.length()
			#if l<minl1:
				#minl1=l
				#n1=i
		
		#for i in range(max(0,n1-100),min(n1+100,atoms.get_number_points())):
			
			#p=atoms.get_vector_at(i)
			#dp=p-p2
			#l=dp.length()
			#if l<minl2:
				#minl2=l
				#n2=i
			
		#print n1,n2
		#if n2>n1:
			#hlx.extend([n1+1,n2-1,h[0],h[1],h[2],h[3],h[4],h[5]])
		#else:
			#hlx.extend([n2+1,n1-1,h[3],h[4],h[5],h[0],h[1],h[2]])
		#print hlx
			
	#for i in range(len(hlx)/8):
		#print hlx[i*8],hlx[i*8+1]
	
	
	atoms.save_pdb_with_helix(options.output,hlx)
	#atoms.save_to_pdb(options.output)
	
	#hlx=[]
	#for h in allhlx:
		#hlx.extend([0,0,h[0],h[1],h[2],h[3],h[4],h[5]])
	
	print hlx
	if options.mapwohelix<>None:
		atoms.remove_helix_from_map(mrc,hlx)
		mrc.write_image(options.mapwohelix)


if __name__ == '__main__':
	main()
	