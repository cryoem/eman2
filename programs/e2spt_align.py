#!/usr/bin/env python
# align all particles to reference and store alignment results
# Author: Steven Ludtke (sludtke@bcm.edu)
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

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
import time
import os
import threading
import queue
from sys import argv,exit
from EMAN2jsondb import JSTask
import numpy as np

def alifn(jsd,fsp,i,a,options):
	t=time.time()
	b=EMData(fsp,i).do_fft()
	b.process_inplace("xform.phaseorigin.tocorner")

	# we align backwards due to symmetry
	if options.verbose>2 : print("Aligning: ",fsp,i)
	c=a.xform_align_nbest("rotate_translate_3d_tree",b,{"verbose":0,"sym":options.sym,"sigmathis":0.1,"sigmato":1.0, "maxres":options.maxres},options.nsoln)
	for cc in c : cc["xform.align3d"]=cc["xform.align3d"].inverse()

	jsd.put((fsp,i,c[0]))
	if options.verbose>1 : print("{}\t{}\t{}\t{}".format(fsp,i,time.time()-t,c[0]["score"]))

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_align.py [options] <subvolume_stack> <reference>
This program is part of the 'new' hierarchy of e2spt_ programs. It performs one iteration of a classical subtomogram refinement, ie -  aligning particles with missing wedge to a reference in 3-D

The reference may be <volume> or <volume>,<n>

If --goldstandard is specified, then even and odd particles will be aligned to different perturbed versions of the reference volume, phase-randomized past the specified resolution."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--goldstandard",type=float,help="If specified, will phase randomize the even and odd references past the specified resolution (in A, not 1/A)",default=0)
	parser.add_argument("--goldcontinue",action="store_true",help="Will use even/odd refs corresponding to specified reference to continue refining without phase randomizing again",default=False)
	#parser.add_argument("--saveali",action="store_true",help="Save a stack file (aliptcls.hdf) containing the aligned subtomograms.",default=False)
	#parser.add_argument("--savealibin",type=int,help="shrink aligned particles before saving",default=1)
	parser.add_argument("--path",type=str,default=None,help="Path to a folder where results should be stored, following standard naming conventions (default = spt_XX)")
	parser.add_argument("--sym",type=str,default="c1",help="Symmetry of the input. Must be aligned in standard orientation to work properly.")
	parser.add_argument("--maxres",type=float,help="Maximum resolution to consider in alignment (in A, not 1/A)",default=0)
	#parser.add_argument("--wtori",type=float,help="Weight for using the prior orientation in the particle header. default is -1, i.e. not used.",default=-1)
	parser.add_argument("--nsoln",type=int,help="number of solutions to keep at low resolution for the aligner",default=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--parallel", type=str,help="Thread/mpi parallelism to use", default=None)
	parser.add_argument("--refine",action="store_true",help="local refinement from xform.align3d in header.",default=False)
	
	parser.add_argument("--refinentry", type=int, help="number of tests for refine mode. default is 8",default=8)
	parser.add_argument("--randphi",action="store_true",help="randomize phi during refine alignment",default=False)
	parser.add_argument("--breaksym",action="store_true",help="symmetry breaking.",default=False)
	parser.add_argument("--breaksymsym",type=str,help="the symmetry to use for breaksym. setting sym to c6 and this to c2 results in a c3 structure. default is the same as sym",default=None)
	parser.add_argument("--rand180",action="store_true",help="randomly add a 180 degree rotation during refine alignment",default=False)
	parser.add_argument("--skipali",action="store_true",help="skip alignment. the program will do nothing. mostly for testing...",default=False)
	parser.add_argument("--maxang",type=float,help="Maximum angular difference for the refine mode. default is 30",default=30)
	parser.add_argument("--maxshift",type=float,help="Maximum shift for the refine mode. default is 16",default=-1)


	(options, args) = parser.parse_args()
	
	#task=SptAlignTask(0,1,2,options)
	#from pickle import dumps,loads,dump,load
	#f=open("task.tmp",'w')
	#dump(task,f)
	#f.close()
	#print(task)
	
	
	#return

	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="spt_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : fls=[0]
		options.path = "spt_{:02d}".format(max(fls)+1)
		try: os.mkdir(options.path)
		except: pass

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : options.iter=1
		else: options.iter=max(fls)+1
		
	if options.parallel==None:
		options.parallel="thread:{}".format(options.threads)
		
	if options.breaksym:
		if options.breaksymsym==None:
			if options.sym=="c1":
				print("cannot break a c1 symmetry. breaksym disabled.")
				options.breaksym=False
			else:
				options.breaksymsym=options.sym
		

	# file may be "name" or "name,#"
	reffile=args[1].split(",")[0]
	try: refn=int(args[1].split(",")[1])
	except: refn=0
	
	NTHREADS=max(options.threads+1,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)
	refnames=[]

	if options.goldcontinue:
		ref=[]
		try:
			refnames=[reffile[:-4]+"_even.hdf", reffile[:-4]+"_odd.hdf"]
			ref.append(EMData(refnames[0],0))
			ref.append(EMData(refnames[1],0))
			
		except:
			print("Error: cannot find one of reference files, eg: ",EMData(reffile[:-4]+"_even.hdf",0))
	else:
		ref=[]
		ref.append(EMData(reffile,refn))
		ref.append(EMData(reffile,refn))

		if options.goldstandard>0 :
			ref[0].process_inplace("filter.lowpass.randomphase",{"cutoff_freq":1.0/options.goldstandard})
			ref[0].process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.goldstandard})
			ref[1].process_inplace("filter.lowpass.randomphase",{"cutoff_freq":1.0/options.goldstandard})
			ref[1].process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.goldstandard})
			refnames=["{}/align_ref_even.hdf".format(options.path), "{}/align_ref_odd.hdf".format(options.path)]
			ref[0].write_image(refnames[0],0)
			ref[1].write_image(refnames[1],0)
			
		else:
			refnames=[reffile, reffile]

	#ref[0]=ref[0].do_fft()
	#ref[0].process_inplace("xform.phaseorigin.tocorner")
	#ref[1]=ref[1].do_fft()
	#ref[1].process_inplace("xform.phaseorigin.tocorner")

	
	#jsd=queue.Queue(0)

	n=-1
	tasks=[]
	if args[0].endswith(".lst") or args[0].endswith(".hdf"):
		#### check if even/odd split exists
		fsps=[args[0][:-4]+"__even.lst",args[0][:-4]+"__odd.lst"]
		
		if os.path.isfile(fsps[0]) and os.path.isfile(fsps[1]):
			print("Using particle list: \n\t {} \n\t {}".format(fsps[0], fsps[1]))
			for eo, f in enumerate(fsps):
				N=EMUtil.get_image_count(f)
				tasks.extend([(f,i,refnames, eo) for i in range(N)])
				
		#### split by even/odd by default
		else:
			N=EMUtil.get_image_count(args[0])
			tasks.extend([(args[0],i,refnames, i%2) for i in range(N)])
			#thrds=[threading.Thread(target=alifn,args=(jsd,args[0],i,ref[i%2],options)) for i in range(N)]
	
	elif args[0].endswith(".json"):
		print("Reading particles from json. This is experimental...")
		js=js_open_dict(args[0])
		keys=sorted(js.keys())
		for k in keys:
			src, ii=eval(k)
			dic=js[k]
			xf=dic["xform.align3d"]
			tasks.append([src, ii, refnames, ii%2, xf])
			


	from EMAN2PAR import EMTaskCustomer
	etc=EMTaskCustomer(options.parallel, module="e2spt_align.SptAlignTask")
	num_cpus = etc.cpu_est()
	options.nowtime=time.time()
	print("{} jobs on {} CPUs".format(len(tasks), num_cpus))
	njob=num_cpus#*4
	
	tids=[]
	for i in range(njob):
		t=tasks[i::njob]
		task=SptAlignTask(t, options)
		tid=etc.send_task(task)
		tids.append(tid)
	
	#print("starting...")
	#for t in tasks:
		#if len(t)>3:
			#task = SptAlignTask(t[0], t[1], t[2], options, t[3])
		#else:
			#task = SptAlignTask(t[0], t[1], t[2], options)
		#tid=etc.send_task(task)
		#tids.append(tid)

	while 1:
		st_vals = etc.check_task(tids)
		if -100 in st_vals:
			print("Error occurs in parallelism. Exit")
			return
		E2progress(logid, np.mean(st_vals)/100.)
		
		if np.min(st_vals) == 100: break
		time.sleep(5)

	#dics=[0]*nptcl
	
	angs={}
	for i in tids:
		rets=etc.get_results(i)[1]
		for ret in rets:
			fsp,n,dic=ret
			if len(dic)==1:
				angs[(fsp,n)]=dic[0]
			else:
				angs[(fsp,n)]=dic
		
	out="{}/particle_parms_{:02d}.json".format(options.path,options.iter)
	if os.path.isfile(out):
		os.remove(out)
	js=js_open_dict(out)
	js.update(angs)
	js.close()

	del etc
	
	
	## here we run the threads and save the results, no actual alignment done here
	#print(len(thrds)," threads")
	#thrtolaunch=0
	#while thrtolaunch<len(thrds) or threading.active_count()>1:
		## If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		## note that it's ok that we wait here forever, since there can't be new results if an existing
		## thread hasn't finished.
		#if thrtolaunch<len(thrds) :
			#while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			#if options.verbose : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
			#thrds[thrtolaunch].start()
			#thrtolaunch+=1
		#else: time.sleep(1)

		#while not jsd.empty():
			#fsp,n,d=jsd.get()
			#angs[(fsp,n)]=d
			#if options.saveali:
				#v=EMData(fsp,n)
				#v.transform(d["xform.align3d"])
				#if options.savealibin>1:
					#v.process_inplace("math.meanshrink",{"n":options.savealibin})
				#v.write_image("{}/aliptcls_{:02d}.hdf".format(options.path, options.iter),n)


	#for t in thrds:
		#t.join()

	E2end(logid)
	
	
def reduce_sym(xf, s):
	sym=parsesym(s)
	xf.invert()
	trans=xf.get_trans()
	xf.set_trans(0,0,0)
	xi=sym.in_which_asym_unit(xf)
	xf.set_trans(trans)
	xf=xf.get_sym(s, -xi)
	xf.invert()
	return xf

class SptAlignTask(JSTask):
	
	
	def __init__(self, data, options):
		
		JSTask.__init__(self,"SptAlign",data,{},"")
		self.options=options
	
	
	def execute(self, callback):
		
		callback(0)
		
		options=self.options
		
		#####[src, ii, refnames, ii%2, xf]
		
		refnames=self.data[0][2]
		refs=[]
		for r in refnames:
			ref=EMData(r,0).do_fft()
			ref.process_inplace("xform.phaseorigin.tocorner")
			refs.append(ref)
			
		if options.breaksym:
			x=Transform()
			nsym=x.get_nsym(options.breaksymsym)
			refasym=[]
			for r in refs:
				rs=[]
				for i in range(nsym):
					rt=r.process("xform",{"transform":x.get_sym(options.breaksymsym, i)})
					rs.append(rt)
				refasym.append(rs)
				
			
		rets=[]
		
		myid=self.data[0][1]
		#print(myid, len(self.data))
		for di,data in enumerate(self.data):
			
			fsp=data[0]
			i=data[1]
			
			if len(data)>4:
				dataxf=data[4]
			else:
				dataxf=None
			
			b=EMData(fsp,i)
			if b["sigma"]==0:
				###skip empty particles.
				c=[{"xform.align3d":Transform(), "score":1}]
				rets.append((fsp,i,c))
				#print("empty particle : {} {}".format(fsp, i))
				continue
				
			b.process_inplace("normalize.edgemean")
			if options.maxres>0:
				b.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1./options.maxres})
			b=b.do_fft()
			b.process_inplace("xform.phaseorigin.tocorner")

			aligndic={"verbose":options.verbose,"sym":options.sym,"sigmathis":0.1,"sigmato":1.0}
			
			if options.refine and (dataxf!=None or b.has_attr("xform.align3d")):
				ntry=options.refinentry
				if dataxf!=None:
					initxf=dataxf
				else:
					initxf=b["xform.align3d"]
				xfs=[initxf]
				
				for ii in range(len(xfs), ntry):
					ixf=initxf.get_params("eman")
					if options.randphi:
						ixf["phi"]=np.random.rand()*360.
					if options.rand180:
						ixf["alt"]=ixf["alt"]+180
					ixf=Transform(ixf)
					
					v=np.random.rand(3)-0.5
					nrm=np.linalg.norm(v)
					if nrm>0:
						v=v/nrm
					else:
						v=(0,0,1)
					xf=Transform({"type":"spin", "n1":v[0], "n2":v[1], "n3":v[2],
							"omega":options.maxang*np.random.randn()/3.0})
					xfs.append(xf*ixf)
				
				## rotate back to the first asym unit
				if options.sym!="c1":
					xfs=[reduce_sym(xf, options.sym) for xf in xfs]
					
				aligndic["initxform"]=xfs
				if options.maxshift<0:
					options.maxshift=16
				aligndic["maxshift"]=options.maxshift
				aligndic["maxang"]=options.maxang
				aligndic["randphi"]=options.randphi
				aligndic["rand180"]=options.rand180
			
			else:
				xfs=[Transform()]
				if options.maxshift>0:
					aligndic["maxshift"]=options.maxshift

			# we align backwards due to symmetry
			if options.verbose>2 : print("Aligning: ",fsp,i)
			
			#print(myid,di,i,time.time()-options.nowtime)
			
	
			ref=refs[data[3]]
			if options.skipali:
				c=[{"xform.align3d":xfs[0].inverse(), "score":-1}]
			else:
				c=ref.xform_align_nbest("rotate_translate_3d_tree",b, aligndic, options.nsoln)
			
			
			
			for cc in c : cc["xform.align3d"]=cc["xform.align3d"].inverse()
			
			if options.breaksym:
				xf=c[0]["xform.align3d"]
				b.process_inplace("xform", {"transform":xf})
				cs=[]
				transxf=[]
				for si in range(nsym):
					ref=refasym[data[3]][si]
					if options.maxshift>0:
						bb=b.align("translational",ref,{"intonly":1, "maxshift":options.maxshift})
						ts=bb["xform.align3d"]
						bb=b.process("xform",{"transform":ts})
					else:
						bb=b
						ts=Transform()
						
					#print(bb["xform.align3d"].get_trans())
					ccc=bb.cmp("fsc.tomo.auto", ref, {"sigmaimgval":3.0, "sigmawithval":0.})
					cs.append(ccc)
					transxf.append(ts)
					
				#print(xf, cs, np.argmin(cs))
				
				
				ci=np.argmin(cs).item()
				txf=transxf[ci]
				x=Transform()
				x=x.get_sym(options.breaksymsym, ci).inverse()

				xf.translate(txf.get_trans())
				xf=x*xf
				#xf=txf.inverse()*xf
				c[0]["xform.align3d"]=xf
				c[0]["score"]=cs[ci]
					
			
			rets.append((fsp,i,c))

			callback(len(rets)*100//len(self.data)-1)

		#print(fsp, i, c[0])
		#callback(100)
		#print(i,c[0]["xform.align3d"])
		#exit()
		return rets
		


if __name__ == "__main__":
	main()

