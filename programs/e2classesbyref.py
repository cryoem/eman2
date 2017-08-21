#!/usr/bin/env python

#
# Author: Steve Ludtke, 7/5/14 
# Copyright (c) 2000- Baylor College of Medicine
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

from math import *
import os
import sys
from EMAN2db import db_check_dict
from EMAN2 import *
import Queue

def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """prog <references> <particles> <classmx> [options]
	
	*** THIS PROGRAM IS NOT YET FUNCTIONAL ***

	This program classifies a set of particles based on a set of references (usually projections). This program makes use of
	bispectral rotational/translational invariants which, aside from computing the invariants, makes the process extremely fast.
	
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--sep", type=int, help="The number of classes a particle can contribute towards (default is 1)", default=1)
	parser.add_argument("--align",type=str,help="specify an aligner to use after classification. Default rotate_translate_tree", default="rotate_translate_tree")
	parser.add_argument("--aligncmp",type=str,help="Similarity metric for the aligner",default="ccc")
	parser.add_argument("--ralign",type=str,help="specify a refine aligner to use after the coarse alignment", default=None)
	parser.add_argument("--raligncmp",type=str,help="Similarity metric for the refine aligner",default="ccc")
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on the local computer")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	if (len(args)<3 ): parser.error("Please specify <references> <particles> <classmx file>")
	
	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.ralign=parsemodopt(options.ralign)
	options.raligncmp=parsemodopt(options.raligncmp)

	
	E2n=E2init(sys.argv, options.ppid)
	
	options.threads+=1		# one extra thread for storing results

	if os.path.exists(args[2]):
		remove_file(args[2])

	nref=EMUtil.get_image_count(args[0])
	nptcl=EMUtil.get_image_count(args[1])

	# get refs and bispectra
	refs=EMData.read_images(args[0])
	refsbsfs=args[0].rsplit(".")[0]+"_bispec.hdf"
	try:
		nrefbs=EMUtil.get_image_count(refsbsfs)
		if nrefbs!=len(refs) : raise Exception
	except:
		print "No good bispecta found for refs. Building"
		com="e2proc2dpar.py {} {} --process filter.highpass.gauss:cutoff_freq=0.01 --process normalize.edgemean --process math.bispectrum.slice:size=32:fp=6".format(args[0],refsbsfs)
		run(com)
	
	refsbs=EMData.read_images(refsbsfs)
	#refsbs=[i.process("filter.highpass.gauss",{"cutoff_freq":0.01}).process("normalize.edgemean").process("math.bispectrum.slice:size=32:fp=6") for i in refs]
	
	# Find particle bispectra
	if "__ctf_flip" in args[1]:
		if "even" in args[1]: bsfs=args[1].split("__ctf_flip")[0]+"__ctf_flip_bispec_even.lst"
		elif "odd" in args[1]: bsfs=args[1].split("__ctf_flip")[0]+"__ctf_flip_bispec_odd.lst"
		else:
			bsfs=args[1].split("__ctf_flip")[0]+"__ctf_flip_bispec.lst"
		nptclbs=EMUtil.get_image_count(bsfs)
		if nptclbs!=nptcl : 
			print nptclbs,nptcl
			raise Exception
		
	# initialize output matrices
	# class, weight, dx,dy,dalpha,flip
	clsmx=[EMData(options.sep,nptcl,1) for i in xrange(6)]
	
	# Actual threads doing the processing
	N=nptcl
	npt=max(min(100,N/(options.threads-2)),1)
	
	jsd=Queue.Queue(0)
	# these start as arguments, but get replaced with actual threads
	thrds=[(jsd,refs,refsbs,args[1],bsfs,options,i,i*npt,min(i*npt+npt,N)) for i in xrange(N/npt+1)]
	
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		if thrtolaunch<len(thrds):
			while (threading.active_count()>=options.threads) : time.sleep(0.1)
			if options.verbose>0 : 
				print "\r Starting thread {}/{}      ".format(thrtolaunch,len(thrds)),
				sys.stdout.flush()
			thrds[thrtolaunch]=threading.Thread(target=clsfn,args=thrds[thrtolaunch])		# replace args
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(0.1)
		
		# return is [N,dict] a dict of image# keyed processed images
		while not jsd.empty():
			rd=jsd.get()
			r=rd[0]
			pt=r[0]
			for i,a in enumerate(r[1:]):
				clsmx[0][i,pt]=a[1]
				clsmx[1][i,pt]=1.0
				clsmx[2][i,pt]=a[2]
				clsmx[3][i,pt]=a[3]
				clsmx[4][i,pt]=a[4]
				clsmx[5][i,pt]=a[5]
			
			if rd[2] :
				thrds[rd[1]].join()
				thrds[rd[1]]=None
			
				if options.verbose>1:
					print "{} done. ".format(rd[1]),

	for i,m in enumerate(clsmx):
		m.write_image(args[2],i)
	

	E2end(E2n)

	print "Classification complete, writing classmx"

def clsfn(jsd,refs,refsbs,ptclfs,ptclbsfs,options,grp,n0,n1):
	from bisect import insort
	
	for i in xrange(n0,n1):
		ptcl=EMData(ptclfs,i)
		ptclbs=EMData(ptclbsfs,i)
		
		# we make a list with the number of total element we want
		best=[(1e30,-1)]*options.sep
		for j,refbs in enumerate(refsbs):
			insort(best,(ptclbs.cmp("ccc",refbs),j))		# insert this comparison in sorted order
			best.pop()								# remove the worst element
		
		ret=[i]
		for b in best:
#			print i,b[0],b[1],options.align,refs[b[1]]["nx"],ptcl["nx"]
			aligned=refs[b[1]].align(options.align[0],ptcl,options.align[1],options.aligncmp[0],options.aligncmp[1])

			if options.ralign!=None: # potentially employ refine alignment
				refine_parms=options.ralign[1]
				refine_parms["xform.align2d"] = aligned.get_attr("xform.align2d")
				refs[b[1]].del_attr("xform.align2d")
				aligned = refs[b[1]].align(options.ralign[0],ptcl,refine_parms,options.raligncmp[0],options.raligncmp[1])
		
			t=aligned["xform.align2d"].inverse()
			prm = t.get_params("2d")
			ret.append((b[0],b[1],prm["tx"],prm["ty"],prm["alpha"],prm["mirror"]))		# bs-sim,cls,tx,ty,alpha,mirror
		
		jsd.put(ret,grp,i==n1-1)	# third value indicates whether this is the final result from this thread
			

def run(command):
	"Mostly here for debugging, allows you to control how commands are executed"

	print "{}: {}".format(time.ctime(time.time()),command)

	ret=launch_childprocess(command)

	# We put the exit here since this is what we'd do in every case anyway. Saves replication of error detection code above.
	if ret !=0 :
		print "Error running: ",command
		sys.exit(1)

	return


if __name__ == "__main__":
    main()
