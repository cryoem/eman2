#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

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

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from math import *
import os
import sys
from EMAN2db import db_check_dict
from EMAN2 import *
import queue
from numpy import array
import traceback

def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """prog <references> <particles> <classmx> [options]
	
	** EXPERIMENTAL **
	
	This program classifies a set of particles based on a set of references (usually projections). This program makes use of
	rotational/translational invariants which, aside from computing the invariants, makes the process extremely fast.
	
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--sep", type=int, help="The number of classes a particle can contribute towards (default is 1)", default=1)
	parser.add_argument("--align",type=str,help="specify an aligner to use after classification. Default rotate_translate_tree", default="rotate_translate_tree")
	parser.add_argument("--aligncmp",type=str,help="Similarity metric for the aligner",default="ccc")
	parser.add_argument("--ralign",type=str,help="specify a refine aligner to use after the coarse alignment", default=None)
	parser.add_argument("--raligncmp",type=str,help="Similarity metric for the refine aligner",default="ccc")
	parser.add_argument("--cmp",type=str,help="Default=auto. The name of a 'cmp' to be used in assessing the aligned images", default="ccc")
	parser.add_argument("--classmx",type=str,help="Store results in a classmx_xx.hdf style file",default=None)
	parser.add_argument("--classinfo",type=str,help="Store results in a classinfo_xx.json style file",default=None)
	parser.add_argument("--classes",type=str,help="Generate class-averages directly. No bad particle exclusion or iteration. Specify filename.",default=None)
	parser.add_argument("--averager",type=str,help="Averager to use for class-averages",default="ctf.weight")	
        parser.add_argument("--invartype",choices=["auto","bispec","harmonic"],help="Which type of invariants to generate: (bispec,harmonic)",default="auto")

	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on the local computer")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	if (len(args)<2 ): parser.error("Please specify <references> <particles>")
	
	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.ralign=parsemodopt(options.ralign)
	options.raligncmp=parsemodopt(options.raligncmp)
	options.cmp=parsemodopt(options.cmp)

        if options.invartype=="auto" :
                try: options.invartype=str(project(["global.invartype"]))
                except: 
                        print("Warning: no project invariant type spectified, using bispectrum")
                        options.invartype="bispec"

	
	E2n=E2init(sys.argv, options.ppid)
	
	options.threads+=1		# one extra thread for storing results

	nref=EMUtil.get_image_count(args[0])
	nptcl=EMUtil.get_image_count(args[1])

	# get refs and invariants
	refs=EMData.read_images(args[0])
	refsbsfs=args[0].rsplit(".")[0]+"_invar.hdf"
	try:
		nrefbs=EMUtil.get_image_count(refsbsfs)
		if nrefbs!=len(refs) :
			print("Reference invariant file too short :",nrefbs,len(refs))
			raise Exception
	except:
#		traceback.print_exc()
#		print("\nError! No good invariants found for refs. Please rerun CTF generate output and set building.")
#		sys.exit(1)

        	if options.invartype=="bispec" :
			com="e2proc2dpar.py {inp} {out} --process filter.highpass.gauss:cutoff_freq=0.01 --process normalize.edgemean --process mask.soft:outer_radius={maskrad}:width={maskw} --process math.bispectrum.slice:size={bssize}:fp={bsdepth} --threads {threads}".format(
			inp=args[0],out=refsbsfs,maskrad=int(refs[0]["nx"]//2.2),maskw=int(refs[0]["nx"]//15),bssize=bispec_invar_parm[0],bsdepth=bispec_invar_parm[1],threads=options.threads)
		else:
			com="e2proc2dpar.py {inp} {out} --process filter.highpass.gauss:cutoff_freq=0.01 --process normalize.edgemean --process mask.soft:outer_radius={maskrad}:width={maskw} --process math.harmonicpow:fp=1 --threads {threads}".format(
			inp=args[0],out=refsbsfs,maskrad=int(refs[0]["nx"]//2.2),maskw=int(refs[0]["nx"]//15),threads=options.threads)

		run(com)
	
	refsbs=EMData.read_images(refsbsfs)
	#refsbs=[i.process("filter.highpass.gauss",{"cutoff_freq":0.01}).process("normalize.edgemean").process("math.bispectrum.slice:size=32:fp=6") for i in refs]
	
	# Find particle invariants
	if "__ctf_flip" in args[1]:
		if "even" in args[1]: bsfs=args[1].split("__ctf_flip")[0]+"__ctf_flip_invar_even.lst"
		elif "odd" in args[1]: bsfs=args[1].split("__ctf_flip")[0]+"__ctf_flip_invar_odd.lst"
		else:
			bsfs=args[1].split("__ctf_flip")[0]+"__ctf_flip_invar.lst"
		try: nptclbs=EMUtil.get_image_count(bsfs)
		except:
			print("Could not get particle count on ",bsfs)
			sys.exit(1)
		if nptclbs!=nptcl : 
			print(nptclbs,nptcl)
			raise Exception("Particle invariant file has wrong particle count")
	else:
		if "even" in args[1]: bsfs=args[1].split("_even")[0]+"_invar_even.lst"
		elif "odd" in args[1]: bsfs=args[1].split("_odd")[0]+"_invar_odd.lst"
		
	### initialize output files
	
	# class, weight, dx,dy,dalpha,flip
	clsmx=[EMData(options.sep,nptcl,1) for i in range(6)]
	
	# JSON style output, classes keyed by class number
	clsinfo={}
	
	# avgs
	if options.classes!=None: 
		options.averager=parsemodopt(options.averager)
		avgrs=[Averagers.get(options.averager[0],options.averager[1]) for i in range(nref)]
	
	# Set up threads
	N=nptcl
	npt=max(min(100,old_div(N,(options.threads-2))),1)
	
	jsd=queue.Queue(0)
	# these start as arguments, but get replaced with actual threads
	thrds=[(jsd,refs,refsbs,args[1],bsfs,options,i,i*npt,min(i*npt+npt,N)) for i in range(old_div(N,npt)+1)]
	
	# standard thread execution loop
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		if thrtolaunch<len(thrds):
			while (threading.active_count()>=options.threads) : time.sleep(0.1)
			if options.verbose>0 : 
				print("\r Starting thread {}/{}      ".format(thrtolaunch,len(thrds)), end=' ')
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
				clsmx[2][i,pt]=a[3]
				clsmx[3][i,pt]=a[4]
				clsmx[4][i,pt]=a[5]
				clsmx[5][i,pt]=a[6]
				
				if options.classinfo!=None:
					try: clsinfo[a[1]].append((pt,a[0],a[2],a[3],a[4],a[5],a[6]))
					except: clsinfo[a[1]]=[(pt,a[0],a[2],a[3],a[4],a[5],a[6])]
						
				if options.classes!=None:
					avgrs[a[1]].add_image(a[7])
			
			if rd[2] :
				thrds[rd[1]].join()
				thrds[rd[1]]=None
			
				if options.verbose>1:
					print("{} done. ".format(rd[1]), end=' ')

	### Write output files
	if options.classmx!=None:
		if os.path.exists(options.classmx): remove_file(options.classmx)
		for i,m in enumerate(clsmx):
			m.write_image(options.classmx,i)
	
	if options.classinfo!=None:
		if os.path.exists(options.classinfo): remove_file(options.classinfo)
		db=js_open_dict(options.classinfo)
		db["input"]=args[1]
		db["inputbs"]=bsfs
		db["refs"]=args[0]
		db["refsbs"]=refsbsfs
		db["classes"]=clsinfo

	if options.classes!=None:
		if os.path.exists(options.classes): remove_file(options.classes)
		empty=EMData(refs[0]["nx"],refs[0]["ny"],1)
		empty.to_zero()
		empty["ptcl_repr"]=0
		for i,avgr in enumerate(avgrs):
			if i in clsinfo:
				avg=avgr.finish()
#				avg.process_inplace("normalize.circlemean",{"radius":avg["ny"]/2-4})
				avg.process_inplace("normalize.toimage",{"to":refs[i],"fourieramp":1,"ignore_lowsig":0.3})
				avg.process_inplace("mask.soft",{"outer_radius":old_div(avg["ny"],2)-4,"width":3})
#				avg.process_inplace("normalize.toimage",{"to":refs[i],"ignore_lowsig":0.75})
				avg["class_ptcl_idxs"]=[p[0] for p in clsinfo[i]]		# particle indices
				quals=array([p[1] for p in clsinfo[i]])
				avg["class_ptcl_qual"]=quals.mean()
				avg["class_ptcl_qual_sigma"]=quals.std()
#				avg["class_qual"]=avg.cmp("frc",refs[i],{"minres":25,"maxres":10})
#				avg["class_qual"]=avg.cmp("ccc",refs[i])	# since we are doing SNR below now, frc seems unnecessary, particularly since ccc is used in e2classaverage
				avg["class_qual"]=old_div(avg.cmp("frc",refs[i],{"minres":30,"maxres":10}),avg.cmp("frc",refs[i],{"minres":100,"maxres":30}))	# Trying something new 2/7/18. This ratio seems pretty effective at identifying bad class-averages. A bit slow, should consider writing something specifically for this
				
				# We compute a smoothed SSNR curve by comparing to the reference. We keep overwriting ssnr to gradually produce what we're after
				ssnr=avg.calc_fourier_shell_correlation(refs[i])
				third=old_div(len(ssnr),3)
				ssnr=[ssnr[third]]*4+ssnr[third:third*2]+[ssnr[third*2-1]]*4	# we extend the list by replication to make the running average more natural
				ssnr=[old_div(sum(ssnr[j-4:j+5]),9.0) for j in range(4,third+4)]		# smoothing by running average
				ssnr=[old_div(v,(1.0-min(v,.999999))) for v in ssnr]						# convert FSC to pseudo SSNR
				avg["class_ssnr"]=ssnr
				
				avg["class_ptcl_src"]=args[1]
				avg["projection_image"]=args[0]
				avg["projection_image_idx"]=i
				try: avg["xform.projection"]=refs[i]["xform.projection"]
				except: pass
				avg.write_image(options.classes,i)
			else:
				empty.write_image(options.classes,i)
		

	E2end(E2n)

	print("Classification complete, writing classmx")

def clsfn(jsd,refs,refsbs_org,ptclfs,ptclbsfs,options,grp,n0,n1):
	from bisect import insort
	
	retali=(options.classes!=None)
	lastdf=-999.0
	lastdfn=-1
	
	for i in range(n0,n1):
		ptcl=EMData(ptclfs,i)
		ptclbs=EMData(ptclbsfs,i)
		
		ctf=ptcl["ctf"]
		if fabs(ctf.defocus-lastdf)>0.1: 
#			print "New DF {} -> {}   ( {} -> {} )".format(lastdf,ctf.defocus,lastdfn,i)
			# note that with purectf=1 this is replacing the image entirely
			ctfim=ptcl.process("math.simulatectf",{"voltage":ctf.voltage,"cs":ctf.cs,"defocus":ctf.defocus,"bfactor":ctf.bfactor,"ampcont":ctf.ampcont,"apix":ctf.apix,"phaseflip":0,"purectf":1})
			if ptclbs.has_attr("is_harmonic_fp"):
				dfmod=ctfim.process("math.harmonicpow",{"fp":1})
			else:
				dfmod=ctfim.process("math.bispectrum.slice",{"size":bispec_invar_parm[0],"fp":bispec_invar_parm[1]})
			#print(ctfim["nx"],dfmod["nx"],refsbs_org[0]["nx"])
			#print(ctfim["ny"],dfmod["ny"],refsbs_org[0]["ny"])
			refsbs=[im.copy() for im in refsbs_org]
			for im in refsbs: im.mult(dfmod)
#			print "New DF done ",ctf.defocus
			lastdf=ctf.defocus
			lastdfn=i
		
		# we make a list with the number of total elements we want
		best=[(1e30,-1)]*(options.sep*5+5)		# we keep 5+5n possible classifications for each desired final output, then pare this down to the best nsep in the next stage
		for j,refbs in enumerate(refsbs):
			try: insort(best,(ptclbs.cmp("ccc",refbs),j))		# insert this comparison in sorted order
			except:
				ptclbs.write_image("debug.hdf",-1)
				refbs.write_image("debug.hdf",-1)
				sys.exit(1)
			best.pop()								# remove the worst element
		
		ret=[i]
		newbest=[]
		for b in best:
#			print i,b[0],b[1],options.align,refs[b[1]]["nx"],ptcl["nx"]
			aligned=refs[b[1]].align(options.align[0],ptcl,options.align[1],options.aligncmp[0],options.aligncmp[1])

			if options.ralign!=None and options.ralign[0]!=None: # potentially employ refine alignment
				refine_parms=options.ralign[1]
				refine_parms["xform.align2d"] = aligned.get_attr("xform.align2d")
				refs[b[1]].del_attr("xform.align2d")
				aligned = refs[b[1]].align(options.ralign[0],ptcl,refine_parms,options.raligncmp[0],options.raligncmp[1])
		
			c=aligned.cmp(options.cmp[0],ptcl,options.cmp[1])
			t=aligned["xform.align2d"].inverse()
#			t=aligned["xform.align2d"]
			prm = t.get_params("2d")
			if retali: insort(newbest,(c,b[1],b[0],prm["tx"],prm["ty"],prm["alpha"],prm["mirror"],ptcl.process("xform",{"transform":t})))		# cls,bs_sim,tx,ty,alpha,mirror,ptcl
			else: insort(newbest,(c,b[1],b[0],prm["tx"],prm["ty"],prm["alpha"],prm["mirror"]))		# cls,bs_sim,cls,tx,ty,alpha,mirror
		
		for j in range(options.sep):
			ret.append(newbest[j])
		jsd.put((ret,grp,i==n1-1))	# third value indicates whether this is the final result from this thread
			

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
