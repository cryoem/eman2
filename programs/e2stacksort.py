#!/usr/bin/env python

#
# Author: Steven Ludtke, 01/03/07 (sludtke@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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

# e2stacksort.py  01/03/07  Steven Ludtke
# This program will sort a stack of images based on some similarity criterion


from EMAN2 import *
from math import *
import os
import sys
from time import time
import random

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <input_stack> <output_stack>
	
	This program will sort a stack of images based on some similarity criterion. Note that byptcl, iterative and reverse are mutually exclusive. 

	It can handle file sequences rather than a single stack file as well, in a limited form. If you have files named var3d_000.mrc, var3d_001.mrc, ...
	you can specify input (and output) files as (for example) var3d_03d.mrc. The '3' indicates the number of digits in the value, the '0'
	means there should be leading zeroes in the name, and the 'd' means it is a decimal number.

	Note that there is no low-memory option for this command, and all images in the sequence are read in, so make sure you have enough RAM."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--simcmp",type=str,help="The name of a 'cmp' to be used in comparing the after optional alignment (default=optvariance:keepzero=1:matchfilt=1)", default="optvariance:keepzero=1:matchfilt=1")
	parser.add_argument("--simalign",type=str,help="The name of an 'aligner' to use prior to comparing the images (default=no alignment)", default=None)
	parser.add_argument("--simmask",type=str, help="A file containing a mask to apply prior to comparisons, to focus the sort on one particular region",default=None)
	parser.add_argument("--reverse",action="store_true",default=False,help="Sort in order of least mutual similarity")
	parser.add_argument("--byptcl",action="store_true",default=False,help="Sort in order of number of particles represented in each class-average. No alignment, shrinking, etc. is performed")
	parser.add_argument("--bykurtosis",action="store_true",default=False,help="Sort by image Kurtosis. No alignment, shrinking, etc. is performed")
	parser.add_argument("--byheader",type=str, help="Uses the named header parameter to sort the images",default=None)
	parser.add_argument("--iterative",action="store_true",default=False,help="Iterative approach for achieving a good 'consensus alignment' among the set of particles") 
	parser.add_argument("--useali",action="store_true",default=False,help="Save aligned particles to the output file, note that if used with shrink= this will store the reduced aligned particles")
	parser.add_argument("--center",action="store_true",default=False,help="After alignment, particles are centered via center of mass before comparison")
	parser.add_argument("--nsort",type=int,help="Number of output particles to generate (mainly for reverse mode)",default=0)
	parser.add_argument("--ninput",type=int,help="Number of input particles to read (first n in the file)",default=0)
	parser.add_argument("--shrink",type=int,help="Reduce the particles for comparisons",default=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
#	parser.add_argument("--tilt", "-T", type=float, help="Angular spacing between tilts (fixed)",default=0.0)
#	parser.add_argument("--maxshift","-M", type=int, help="Maximum translational error between images (pixels), default=64",default=64.0)
#	parser.add_argument("--mode",type=str,help="centering mode 'modeshift', 'censym' or 'region,<x>,<y>,<clipsize>,<alisize>",default="censym")
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	
	if options.iterative+options.byptcl+options.reverse>1 :
		parser.error("byptcl, iterative and reverse are mututally exclusive")

	print "Beginning image sort/alignment"
	E2n=E2init(sys.argv,options.ppid)

	if options.simalign : options.simalign=parsemodopt(options.simalign)
	else: options.simalign=[None,None]
	if options.simcmp : options.simcmp=parsemodopt(options.simcmp)
	
	if options.simmask!=None : options.simmask=EMData(options.simmask,0)
	
	
	# read all images
	if "%" not in args[0] :
		a=EMData.read_images(args[0])
	else :
		a=[]
		i=0
		while 1:
			try: im=EMData(args[0]%i,0)
			except: break
			a.append(im)
			i+=1
	
	# technically we read them all then truncate the list. inefficient...
	if options.ninput>0 : a=a[:options.ninput]
			
	if options.nsort<1 : options.nsort=len(a)
	if options.byptcl : 
		b=sortstackptcl(a,options.nsort)
		if options.reverse : b.reverse()
	elif options.bykurtosis:
		b=sortstackkurt(a,options.nsort)
		if options.reverse : b.reverse()
	elif options.byheader!=None:
		b=sortstackheader(a,options.nsort,options.byheader)
		if options.reverse : b.reverse()
	elif options.iterative: b=sortstackiter(a,options.simcmp[0],options.simcmp[1],options.simalign[0],options.simalign[1],options.nsort,options.shrink,options.useali,options.center,options.simmask)
	elif options.reverse: b=sortstackrev(a,options.simcmp[0],options.simcmp[1],options.simalign[0],options.simalign[1],options.nsort,options.shrink,options.useali,options.center,options.simmask)
	else : b=sortstack(a,options.simcmp[0],options.simcmp[1],options.simalign[0],options.simalign[1],options.nsort,options.shrink,options.useali,options.center,options.simmask)
	
	if "%" not in args[1] :
		for i,im in enumerate(b): im.write_image(args[1],i)
	else :
		for i,im in enumerate(b): im.write_image(args[1]%i,0)

	E2end(E2n)

def sortstackiter(stack,cmptype,cmpopts,align,alignopts,nsort,shrink,useali,center,mask):
	"""The goal here is to provide a "consensus orientation" for each particle, despite the impossibility
	of this task for a distribution of projections on a sphere. Since in most cases a "perfect" answer is
	impossible anyway, we avoid computing the full similarity matrix and use this iterative technique
	to get a decent result. """
	stackshrink=[i.copy() for i in stack]
	if (shrink>1) :
		for i in stackshrink: i.process_inplace("math.meanshrink",{"n":shrink})

	if center : 
		for i in stackshrink: i.process_inplace("xform.centerofmass")

	# initialize the connectivity with a random linear chain
	for i,im in enumerate(stack): 
		im.set_attr("align_target",(i+1)%len(stack))
		ima=stackshrink[i].align(align,stackshrink[(i+1)%len(stack)],alignopts)
		if mask!=None : ima.mult(mask)
		c=stackshrink[i].cmp(cmptype,ima,cmpopts)
		im.set_attr("align_qual",c)
		im.set_attr("aligned",0)
		
	# now we iterate to improve the similarity of each particle to the reference its being aligned to
	changes=1
	while (changes) :
		changes=0
		for i,im in enumerate(stack):
			# first we try comparing to the reference of the particle we are currently linked to
			at=stack[im.get_attr("align_target")].get_attr("align_target")
			if at!=i :
				ima=stackshrink[i].align(align,stackshrink[at],alignopts)
#				c=ima.cmp(cmptype,stackshrink[at],cmpopts)
				if mask!=None : ima.mult(mask)
				c=stackshrink[at].cmp(cmptype,ima,cmpopts)
				if c<im.get_attr("align_qual") :
					im.set_attr("align_qual",c)
					im.set_attr("align_target",at)
					changes+=1
				
			# then we also try a random particle from the list
			at=i
			while at==i or at==im.get_attr("align_target") : at=random.randint(0,len(stack)-1)
			ima=stackshrink[i].align(align,stackshrink[at],alignopts)
#			c=ima.cmp(cmptype,stackshrink[at],cmpopts)
			if mask!=None : ima.mult(mask)
			c=stackshrink[at].cmp(cmptype,ima,cmpopts)
			if c<im.get_attr("align_qual"):
				im.set_attr("align_qual",c)
				im.set_attr("align_target",at)
				changes+=1
				
		print changes, "changed"

	# a list for each particle of particles aligning to this particle
	br=[[] for i in stack]
	for i,im in enumerate(stack):
		br[im.get_attr("align_target")].append(i)

	# tack the length in as the first element of the list and sort so the most 'aligned to' particles are first
	bn=[(len(j),i,j) for i,j in enumerate(br)]
	bn.sort(reverse=1)

	# a single particle can stay in its original orientation. We will use the most referenced particle for this
	stack[bn[0][1]].set_attr("aligned",1)
	if center : 
		stack[bn[0][1]].process_inplace("xform.centerofmass")


	ret=[stack[bn[0][1]]]
	# now align the particles in the order we find them in the referenced sublists
	for i in bn:
		for j in i[2]:
			recursealign(stack,j,align,alignopts,ret)

	# This is for debugging. Write out 'clusters' of aligned image matches
#	c=0
#	for i in bn:
#		if len(i[2])>1 :
#			stack[i[1]].write_image("cluster.%d.hdf"%c,-1)
#			for j in i[2]: stack[j].write_image("cluster.%d.hdf"%c,-1)
#			c+=1

	# sort in order of the number of particles aligned to this one
	return ret
	

def recursealign(stack,src,align,alignopts,ret):
	"""This is used to align one particle in stack to another, specified by the "align_target" attribute.
 If the target of 'src' is not already aligned ("aligned" attribute set), recursealign will be called on it first.
 ret is an empty or almost empty list which will be built with aligned particles in order of alignment  """
	if stack[src].get_attr("aligned") : return
	else : stack[src].set_attr("aligned",1)		# we set aligned now to prevent infinite recursion if there is a loop

	trg=stack[src].get_attr("align_target")
	if not stack[trg].get_attr("aligned") : recursealign(stack,trg,align,alignopts,ret)
	
	stack[src]=stack[src].align(align,stack[trg],alignopts)
	stack[src].set_attr("aligned",1)
#	print "aligned", src
	ret.append(stack[src])

def sortstackrev(stack,cmptype,cmpopts,align,alignopts,nsort,shrink,useali,center,mask):
	"""Sorts a list of images in order of LEAST similarity"""
	
	stackshrink=[i.copy() for i in stack]
	if (shrink>1) :
		for i in stackshrink: i.process_inplace("math.meanshrink",{"n":shrink})
	
	ret=[stack[0]]
	rets=[stackshrink[0]]
	if stack[0].get_attr_default("ptcl_repr",0)>0 : check_rep=1    # if the images have ptcl_repr, then we want to ignore those with 0 images represented
	else : check_rep=0
	del stack[0]
	del stackshrink[0]
	print "nsort=",nsort
	while (len(stack)>0 and len(ret)<nsort) :
		best=(0,-1)
		for i in range(len(stackshrink)):
			if check_rep and stackshrink[i].get_attr_default("ptcl_repr",1)<=0 : continue
			c=1.0e38
			cj=-1
			ci=None
			for j,r in enumerate(rets):			# compare to all existing solutions, and use the MOST similar value
				if align : ims=stackshrink[i].align(align,r,alignopts)
				else : ims=stackshrink[i].copy()
				if center : ims.process_inplace("xform.centerofmass")
				if mask!=None : ims.mult(mask)
				cc=r.cmp(cmptype,ims,cmpopts)
				if cc<c : c,cj,ci=cc,j,ims
#			print "\t%d. %1.3g (%d)"%(i,c,cj)
			if c>best[0] or best[1]<0 : best=(c,i,ims)
		if useali : 
			ret.append(best[2])
			rets.append(best[2])
		else : 
			ret.append(stack[best[1]])
			rets.append(stackshrink[best[1]])
		del stack[best[1]]
		del stackshrink[best[1]]
		print "%d.\t%d  (%1.4f)"%(len(ret)-1,best[1],best[0])

	return ret

def sortstack(stack,cmptype,cmpopts,align,alignopts,nsort,shrink,useali,center,mask):
	"""Sorts a list of images based on a standard 'cmp' metric. cmptype is the name
	of a valid cmp type. cmpopts is a dictionary. Returns a new (sorted) stack.
	The original stack is destroyed."""
	
	stackshrink=[i.copy() for i in stack]
	if (shrink>1) :
		for i in stackshrink: i.process_inplace("math.meanshrink",{"n":shrink})
	ret=[stack[0]]
	rets=[stackshrink[0]]
	del stack[0]
	del stackshrink[0]
	while (len(stack)>0 and len(ret)<nsort) :
		best=(1.0e38,-1)
		for i,ims in enumerate(stackshrink):
			if align : ims=ims.align(align,rets[-1],alignopts)
			else: ims=ims.copy()
			if center : ims.process_inplace("xform.centerofmass")
			if mask!=None: ims.mult(mask)
			c=rets[-1].cmp(cmptype,ims,cmpopts)+ims.cmp(cmptype,rets[-1],cmpopts)	# symmetrize results
			if c<best[0] or best[1]<0 : best=(c,i,ims)
		if useali : 
			ret.append(best[2])
			rets.append(best[2])
		else : 
			ret.append(stack[best[1]])
			rets.append(stackshrink[best[1]])
		del stack[best[1]]
		del stackshrink[best[1]]
		print "%d.\t%d  (%1.4f)"%(len(ret)-1,best[1],best[0])

	return ret

def sortstackkurt(stack,nsort):
	"""Sorts a list of images based on the number of particles each image represents"""

	stack.sort(key=lambda B:B.get_attr("kurtosis"),reverse=True)

	return stack[:nsort]

def sortstackheader(stack,nsort,header):
	"""Sorts a list of images based on the number of particles each image represents"""

	stack.sort(key=lambda B:B.get_attr(header))

	return stack[:nsort]


def sortstackptcl(stack,nsort):
	"""Sorts a list of images based on the number of particles each image represents"""

	stack.sort(key=lambda B:B.get_attr("ptcl_repr"),reverse=True)

	return stack[:nsort]
	
if __name__ == "__main__":
    main()
