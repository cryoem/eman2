#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

#
# Author: Steve Ludtke, 1/15/18 
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

from future import standard_library
standard_library.install_aliases()
from builtins import range
from math import *
import os
import sys
from EMAN2 import *
import queue
from numpy import array

def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """prog <particles> <ref_bispec1> [ref_bispec2] ... [options]
	
	** EXPERIMENTAL **
	
	This program computes the similarity between a particle and several sets of projection references
	using bispectral methods, very much like e2classesbyref. The goal is to better understand why perticle
	orientation may be ambiguous over multiple iterations.
	
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ptcln", type=int, help="Particle number in the input file to compute similarities for", default=1)
	parser.add_argument("--output",type=str,help="Output filename (text file)", default="ptclplot.txt")
	#parser.add_argument("--align",type=str,help="specify an aligner to use after classification. Default rotate_translate_tree", default="rotate_translate_tree")
	#parser.add_argument("--aligncmp",type=str,help="Similarity metric for the aligner",default="ccc")
	#parser.add_argument("--ralign",type=str,help="specify a refine aligner to use after the coarse alignment", default=None)
	#parser.add_argument("--raligncmp",type=str,help="Similarity metric for the refine aligner",default="ccc")
	#parser.add_argument("--cmp",type=str,help="Default=auto. The name of a 'cmp' to be used in assessing the aligned images", default="ccc")
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on the local computer")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	if (len(args)<2 ): parser.error("Please specify <particles> <refs1> [refs2] ...")
	
	#options.align=parsemodopt(options.align)
	#options.aligncmp=parsemodopt(options.aligncmp)
	#options.ralign=parsemodopt(options.ralign)
	#options.raligncmp=parsemodopt(options.raligncmp)
	#options.cmp=parsemodopt(options.cmp)

	
	E2n=E2init(sys.argv, options.ppid)
	
	options.threads+=1		# one extra thread for storing results

	nfiles=len(args)-1
	nref=EMUtil.get_image_count(args[1])
	nptcl=EMUtil.get_image_count(args[0])
	if options.ptcln>=nptcl : 
		print("Error: particle number outside range of input file")
		sys.exit(0)

	# get refs and bispectra
	refsbs=[EMData.read_images(args[i]) for i in range(1,len(args))]
	for i in range(1,len(refsbs)):
		if len(refsbs[i])!=len(refsbs[0]) :
			print("ERROR: {} has {} references and {} has {}".format(args[0],len(refsbs[0]),args[i],len(refsbs[i])))
			
	# try to get Euler angles
	pfsp=args[1].replace("_bispec","")
	try:
		eulers=[EMData(pfsp,i,True)["xform.projection"] for i in range(len(refsbs[0]))]
	except:
		print("ERROR: cannot read Eulers from ",pfsp)
		sys.exit(1)
	
	# Find particle bispectra
	if "__ctf_flip" in args[0]:
		if "even" in args[0]: bsfs=args[0].split("__ctf_flip")[0]+"__ctf_flip_bispec_even"
		elif "odd" in args[0]: bsfs=args[0].split("__ctf_flip")[0]+"__ctf_flip_bispec_odd"
		else:
			bsfs=args[0].split("__ctf_flip")[0]+"__ctf_flip_bispec"
	else:
		if "even" in args[1]: bsfs=args[0].split("_even")[0]+"_bispec_even"
		elif "odd" in args[1]: bsfs=args[0].split("_odd")[0]+"_bispec_odd"
		
	bsfs+=args[0][-4:]
	
	# Set up threads
	N=len(refsbs)
	rslt={}
	
	ptcl=EMData(args[0],options.ptcln)
	ptclbs=EMData(bsfs,options.ptcln)
	jsd=queue.Queue(0)
	# these start as arguments, but get replaced with actual threads
	thrds=[(jsd,refsbs,ptcl,ptclbs,options.ptcln,i) for i in range(N)]
	
	
	# standard thread execution loop
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		if thrtolaunch<len(thrds):
			while (threading.active_count()>=options.threads) : time.sleep(0.1)
			if options.verbose>0 : 
				print("\r Starting thread {}/{}      ".format(thrtolaunch,len(thrds)), end=' ')
				sys.stdout.flush()
			thrds[thrtolaunch]=threading.Thread(target=simfn,args=thrds[thrtolaunch])		# replace args
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(0.1)
		
		# return is [N,dict] a dict of image# keyed processed images
		while not jsd.empty():
			rd=jsd.get()
			rslt[rd[0]]=rd[1]
			
			thrds[rd[0]].join()
			thrds[rd[0]]=None
			
			if options.verbose>1:
				print("{} done. ".format(rd[0]), end=' ')
			

	#print(rslt)
	# setup output file
	out=open(options.output,"w")
	out.write("# Particle {} from {}\n# ".format(options.ptcln,args[0]))
	for i in args[1:]: out.write(i+", ")
	out.write("\n")
	
	for i in range(len(refsbs[0])):
		eul=eulers[i].get_rotation()
		out.write("{}\t{}".format(eul["alt"],eul["az"]))
		for j in range(N): out.write("\t{}".format(rslt[j][i]))
		out.write("\n")

	E2end(E2n)

	print("Classification complete, writing classmx")

def simfn(jsd,refsbs_org,ptcl,ptclbs,ptcln,thrn):
	
	ctf=ptcl["ctf"]
	refsbs=[im.process("math.simulatectf",{"voltage":ctf.voltage,"cs":ctf.cs,"defocus":ctf.defocus,"bfactor":ctf.bfactor,"ampcont":ctf.ampcont,"apix":ctf.apix,"phaseflip":0,"bispectrumfp":1}) for im in refsbs_org[thrn]]
	
	ret=[]
	for j,refbs in enumerate(refsbs):
		ret.append(ptclbs.cmp("ccc",refbs))
		
	jsd.put((thrn,ret))
			

def run(command):
	"Mostly here for debugging, allows you to control how commands are executed"

	print("{}: {}".format(time.ctime(time.time()),command))

	ret=launch_childprocess(command)

	# We put the exit here since this is what we'd do in every case anyway. Saves replication of error detection code above.
	if ret !=0 :
		print("Error running: ",command)
		sys.exit(1)

	return


if __name__ == "__main__":
    main()
