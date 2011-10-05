#!/usr/bin/env python

#
# Author: Steven Ludtke, 08/28/2006 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

# ssematch.py    Steven Ludtke, Matt Baker
# This tries to match sequence based secondary structure prediction to
# ssehunter results, based on length and distance constraints in each case

from EMAN2 import *
from math import *
import time
import os
import sys
from pprint import pprint

pl=()

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <ssehunter_file> <sec_struct_file>

	WARNING: This program is experimental and still under development. Contact mbaker@bcm.edu for more information.	
		
	ssehunter file is in ??? format
	sec_struct_file is a string of -,H,E defining per-residue predicted structure"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_argument("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
#	parser.add_argument("--auto","-A",type=str,action="append",help="Autobox using specified method: circle, ref, grid",default=[])
#	parser.add_argument("--threshold","-T",type=float,help="Threshold for keeping particles. 0-4, 0 excludes all, 4 keeps all.",default=2.0)
	parser.add_argument("--maxbad","-M",type=int,help="Maximumum number of unassigned helices",default=2)
	parser.add_argument("--minhelix","-H",type=int,help="Minimum residues in a helix",default=6)
	parser.add_argument("--maxpairerr","-E",type=float,help="Maximum error match between pairs of helices, default=50",default=10.0)
	parser.add_argument("--skelpath","-K",type=str,help="Optional (recommended) output from the e2skelpath.py program")
#	parser.add_argument("--lengthmatchmatrix",type=str,help="Writes an image containing an exhaustive comparison of predicted vs SSE helix lengths as a matrix",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	logid=E2init(sys.argv,options.ppid)
	
	ssematch(args[0],args[1],options)

def ssematch2(ssehfsp,sspredfsp,options):
	sseh=readsseh(ssehfsp)
	sspred=readsspred(sspredfsp,options.minhelix)
	
	for i in sseh[0]: print "%d "%int(i/1.5),
	print
	
	skel=readconnect(options.skelpath,len(sseh[0]))
	try:
		skel=readconnect(options.skelpath,len(sseh[0]))
		sseh=(sseh[0],skel)
		print "Skeletonization results read, %d paths"%len(skel)
#		pprint(skel)
	except: pass
	
# 	if options.lengthmatchmatrix:
# 		lengthmatrix(sspred,sseh,options.lengthmatchmatrix)
	
	print "%d predicted helices    %d helices in density"%(len(sspred),len(sseh[0]))
	for i in sspred: print "%4d "%int(i[0]/1.5),
	print ""

	

def ssematch(ssehfsp,sspredfsp,options):
	"older algorithm that matches pairs of predicted helices to pairs of SSEhunter helices"
	
	sseh=readsseh(ssehfsp)
	sspred=readsspred(sspredfsp,options.minhelix)
	
	for i in sseh[0]: print "%d "%int(i/1.5),
	print
	
	skel=readconnect(options.skelpath,len(sseh[0]))
	try:
		skel=readconnect(options.skelpath,len(sseh[0]))
		sseh=(sseh[0],skel)
		print "Skeletonization results read, %d paths"%len(skel)
#		pprint(skel)
	except: pass
	
#	if options.lengthmatchmatrix:
#		lengthmatrix(sspred,sseh,options.lengthmatchmatrix)
	
	print "%d predicted helices    %d helices in density"%(len(sspred),len(sseh[0]))
	for i in sspred: print "%4d "%int(i[0]/1.5),
	print ""
	
	# get lists of possible pairwise assignments and quality assessment for each
	pairqual={}
	for i in range(len(sspred)-1):
		pairqual[i]=findpairs(i,sspred,sseh,options.maxpairerr)
		print "%4d "%len(pairqual[i]),
	print
		
	# This is where we generate all of the final answers
	all=[]
	recursesoln(pairqual,[],[],[],all,options.maxbad)
	
	out=file("ssematch.out","w")
	for i in all:
		out.write("%f\t%s\t%s\n"%(i[0],str(i[1]),str(i[2])))
	out.close()
	
	print all[0]
	
	print len(all)

def recursesoln(pairqual,tot,soln,ends,all,maxbad):
#	print soln,"\n",ends
	# first round, try all 1st level pairs
	if len(soln)==0 :
		for j,i in enumerate(pairqual[0]): 
			recursesoln(pairqual,[i[0]],[i[1],i[2]],[i[3],i[4]],all,maxbad)
			print "%d/%d"%(j,len(pairqual[0]))
		return
		
	# if we get here, we're done
	if len(soln)==len(pairqual)+1:
		v=sum(tot)/len(tot)
		try:
			if v<min(all)[0]: 
				print v,soln,"\n",ends
				for i in tot: print "%1.2f "%i,
				print

		except: print v,soln,"\n",ends
		all.append((v,soln,ends))
		return
	
	tries=0
	for i in pairqual[len(soln)-1]:
		# three tests. If the previous element isn't undefined, it must match
		# and, the next choice in the series must not already be used
#		if (soln[-1]!=-1 and i[1]!=soln[-1]) or (ends[-1]==i[3])or i[2] in soln : continue		# next one in series already assigned
		if (soln[-1]!=-1 and i[1]!=soln[-1]) or i[2] in soln : continue		# next one in series already assigned
		if ends[-1]==i[3] : 
			mm=170
		else: mm=0
		tries+=1
		try: minq=min(minq,i[0])
		except: minq=i[0]
		recursesoln(pairqual,tot+[i[0]+mm],soln+[i[2]],ends+[i[4]],all,maxbad)
	
	# if we didn't find even one good assignment, we skip this helix (unless we've skipped too many)
	if (tries==0 or minq>32) and soln.count(-1)<maxbad:
		recursesoln(pairqual,tot,soln+[-1],ends+[-1],all,maxbad)
	
def findpairs(p1,sspred,sseh,maxpe):
	"""This will generate a sorted list of possible pair assignments. Assigns the
	predicted helices p1 and p1+1 to two helices from sseh. Returns a sorted list of:
	(error,s1,s2)		lower error is a better match """
	
	ssemin=sseh[1]
	
	poss=[]
	for s1 in range(len(sseh[0])):
		for s2 in range(len(sseh[0])):
			for a in (0,1):
				for b in (0,1):
					if s1==s2 or ssemin[a][b][s1][s2]<0 or ssemin[a][b][s1][s2]>sspred[p1+1][1]+9.0: continue
					# error includes squared length mismatches and a term downweighting long distances between helices
#					err=sqrt((sspred[p1][0]-sseh[0][s1])**2+(sspred[p1+1][0]-sseh[0][s2])**2)
					err=sseh[0][s1]*(tanh(fabs(sspred[p1][0]-sseh[0][s1])-6)+1)/30.0+sseh[0][s2]*(tanh(fabs(sspred[p1+1][0]-sseh[0][s2])-6)+1)/30.0
#					if ssemin[a][b][s1][s2]/sspred[p1+1][1]>.75: err+=(4.0*(ssemin[a][b][s1][s2]/sspred[p1+1][1]-.75))**2
#					err+=2.0*fabs(ssemin[a][b][s1][s2]/sspred[p1+1][1]-1.0)
#					if ssemin[s1][s2][1]/sspred[p1+1][1]>.75: err+=(16.0*(ssemin[s1][s2][1]/sspred[p1+1][1]-.75))**2
					poss.append((err,s1,s2,a,b))
	poss.sort()
	if len(poss)==0: return poss
	for i,v in enumerate(poss):
		if v[0]>maxpe : break
	if i==0: i+=1
	return poss[:i]
	
def lengthmatrix(sspred,sseh,fsp):
	"""Compares each predicted length to each SSE helix length in a matrix
	stored as a 2-D image."""
	out=EMData()
	out.set_size(len(sspred),len(sseh[0]),1)
	
	for p in range(len(sspred)):
		for h in range(len(sseh[0])):
			out.set_value_at(p,h,0,fabs(sspred[p][0]-sseh[0][h]))
	
	out.write_image(fsp)
	
def readsseh(fsp):
	"""reads a ssehunter output file, returns 2 results:
	a list of helix lengths
	a distance matrix (h1,h2):(1|2,1|2,dist)  h1<h2
	
	This should really use skeletonization results to generate the distance
	matrix based on only connected h1,h2 pairs"""
	
	# this makes a list of all lines starting with ALPHA
	lns=[i for i in file(fsp,"r").readlines() if i[:5]=="ALPHA"]
	
	# this makes a vector of 6 floats (x0,y0,z0,x1,y1,z1) for each line
	hlx=[[float(j) for j in i.split()[5:]] for i in lns]
	
	# list of lengths of each helix
	lenlist=[sqrt((i[0]-i[3])**2+(i[1]-i[4])**2+(i[2]-i[5])**2) for i in hlx]
	
	# min length matrix (this needs to be replaced with connectivity info from skeletonization)
	# [h1][h2] -> ((0|1,0|1),len)
	lenmx=[[[] for i in range(len(hlx))] for i in range(len(hlx))]
	for h1 in range(len(hlx)):
		for h2 in range(len(hlx)):
			if h1>=h2 : continue
			k=[(0,0),(0,1),(1,0),(1,1)]
			for i in k:
				l=sqrt((hlx[h1][i[0]*3]-hlx[h2][i[1]*3])**2+(hlx[h1][i[0]*3+1]-hlx[h2][i[1]*3+1])**2+(hlx[h1][i[0]*3+2]-hlx[h2][i[1]*3+2])**2)
				try: 
					if l<lenmx[h1][h2][1] : 
						lenmx[h1][h2]=(i,l)
						lenmx[h2][h1]=(i,l)
				except: 
						lenmx[h1][h2]=(i,l)
						lenmx[h2][h1]=(i,l)
	
	return (lenlist,lenmx)
					
def readconnect(fsp,nel):
	"""reads the results from e2skelpath, a list of connectivity and pathlengths between all ssehunter helices
	returns the same results as the second component of readsseh(), and should replace it when possible
	nel is the number of helices in the SSEhunter file supplied to e2skelpath"""
	
	lns=file(fsp,"r").readlines()
	lns=[i.split() for i in lns]
	lns=[(int(i[1]),int(i[3]),int(i[0]),int(i[2]),float(i[4])) for i in lns]	# reorder as helixA#,helixB#,endA,endB,pathlen
	
	# ret is a 4 dimensional array. The first 2 dimensions are either 0 or 1 and represent the
	# end of the helix (arbitrarily assigned by ssehunter). The second two dimensions
	# represent the helix numbers. ie - helix 8 end 1 connected to helix 4 end 0 would be ret[1][0][8][4] or ret[0][1][4][8]
	ret=[[0,0],[0,0]]

	ret[0][0]=[[-1.0 for i in range(nel)] for i in range(nel)]
	ret[0][1]=[[-1.0 for i in range(nel)] for i in range(nel)]
	ret[1][0]=[[-1.0 for i in range(nel)] for i in range(nel)]
	ret[1][1]=[[-1.0 for i in range(nel)] for i in range(nel)]
	for i in lns:
		ret[i[2]][i[3]][i[0]][i[1]]=i[4]
	
	# Symmetrize the distance matrices
	# it is reasonable that the file results are not symmetric due to the discrete
	# path-tracing algorithm used
	for a in (0,1):
		for b in (0,1):
			for x in range(nel):
				for y in range(x):
					t=max(ret[a][b][x][y],ret[b][a][y][x])		# pick the largest of the 2 symmetric values
					ret[a][b][x][y]=t
					ret[b][a][y][x]=t
	
	return ret

def readsspred(fsp,minhelix):
	"""reads a file containing a sequence of -,H,E characters, 1 per residue.
	strips out extraneous whitespace. Produces an ordered list of helix lengths and the maximum distance
	to the previous helix 
	Probably should have a better measure of probable distance vs # residues."""
	seq=file(fsp,"r").read()
	
	# slow, but doesn't matter much here
	seq="".join([i for i in seq if i in ("H","E","-")])
	
	# construct the list by counting H and not H letter sequences
	lenlist=[]
	inh=0
	h1=0
	for i,l in enumerate(seq):
		if inh:
			if l!="H":
				if i-h0<minhelix:
					inh=0
					continue
				lenlist.append((1.5*(i-h0),3.8*(h0-h1)))
				h1=i-1
				inh=0
		else:
			if l=="H":
				h0=i
				inh=1

	return lenlist

if __name__ == "__main__":
	main()
