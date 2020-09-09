#!/usr/bin/env python
# average selected subset of particles

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
import numpy as np
from EMAN2jsondb import JSTask
#def rotfncompete(jsd,avgs,fsp,fspn,a,sym,refs,shrinkrefs,maxtilt,wedgesigma,shrink,maxres,simthr2,verbose):
def rotfncompete(jsd,avgs,fsp,fspn,a,refs, shrinkrefs, options):
	
	"""Averaging thread. 
	avgs are n existing Averagers, 
	fsp,i is the particle being averaged
	a is the Transform to put the particle in the correct orientation
	sym is the symmetry for replication of the particle
	refs are the n alignment references for competitive averaging
	maxtilt can optionally better enforce missing wedge exclusion
	"""
	
	shrink=options.shrinkcompare
	if shrink<2: shrink=0
	b=EMData(fsp,fspn).process("normalize.edgemean")
	if options.maxres>0: b.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,options.maxres)})
	if options.maxtilt<90.0 :
		bf=b.do_fft()
		bf.process_inplace("mask.wedgefill",{"thresh_sigma":0.0,"maxtilt":options.maxtilt})
		b=bf.do_ift()
	b.process_inplace("xform",{"transform":a})
	if shrink : bs=b.process("math.meanshrink",{"n":shrink})
	else: bs=b
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(options.sym)
	best=(1.0e50,None,None)
	for r,ref in enumerate(shrinkrefs):
		for i in range(nsym):
			c=bs.process("xform",{"transform":xf.get_sym(options.sym,i)})
			if options.noali:
				d=c
			else:
				d=c.align("translational",ref)
			score=d.cmp("fsc.tomo.auto",ref,{"sigmaimgval":options.wedgesigma,"sigmawithval":0.5})
			if score<best[0] : best=(score,r,d,i)
	
	if shrink:
		c=b.process("xform",{"transform":xf.get_sym(options.sym,best[3])})
		if options.noali:
			d=c
		else:
			d=c.align("translational",refs[best[1]])
	else :
		d=best[2]
	
	if best[0]<options.simthr2 : 
		avgs[best[1]].add_image(d)
		print("{} -> ref {} sym {}   {}".format(fspn,best[1],best[3],best[0]))
	else: print("** {} -> ref {} sym {}   {}".format(fspn,best[1],best[3],best[0]))
	jsd.put((fspn,best[0],best[1],best[3]))


class SptavgmultTask(JSTask):
	
	def __init__(self, inp, refs, options):
		
		data={"data":inp, "refnames":refs}
		JSTask.__init__(self,"Sptavg",data,{},"")
		self.options=options
	
	
	def execute(self, callback):
		
		callback(0)
		data=self.data["data"]
		refnames=self.data["refnames"]
		options=self.options
		sz=options.boxsz
		
		
		refs=[]
		if options.maskclass:
			maskclass=EMData(options.maskclass)
		
		if options.randnclass<=0:
			for r in refnames:
				if "," in r: ref=EMData(r.split(",")[0],int(r.split(",")[-1]))
				else: ref=EMData(r,0)
				if options.maskclass:
					ref.mult(maskclass)
				ref=ref.do_fft()
				ref.process_inplace("xform.phaseorigin.tocenter")
				#ref.process_inplace("xform.fourierorigin.tocenter")
				refs.append(ref)
			nref=len(refs)
		else:
			nref=options.randnclass
		
		avgrs=[]
		nrmvls=[]
		for i in range(nref):
			normvol=EMData((sz//2+1)*2, sz, sz)
			avgr=Averagers.get("mean.tomo",{"thresh_sigma":options.wedgesigma, "doift":0, "normout":normvol})
			avgrs.append(avgr)
			nrmvls.append(normvol)
		
		stats=[]
		for ii,dt in enumerate(data):
			fsp, i, xf=dt
			#xf.set_trans(np.round(xf.get_trans()).tolist())
			b=EMData(fsp,i)
			b.process_inplace("normalize.edgemean")
			b=b.do_fft()
			if options.maxtilt<90.0 :
				b.process_inplace("mask.wedgefill",{"thresh_sigma":0.0,"maxtilt":options.maxtilt})
			
			b.process_inplace("xform.phaseorigin.tocorner")
			
			x0=Transform()
			nsym=x0.get_nsym(options.sym)
			best=[1.0e50,None,None]
			bxfs=[]
			for k in range(nsym):
				x=x0.get_sym(options.sym, k)
				bxfs.append(b.process("xform",{"transform":x*xf}))
				
			for k, bxf in enumerate(bxfs):
				for r in range(nref):
					if options.randnclass<=0:
						score=bxf.cmp("fsc.tomo.auto", refs[r], {"sigmaimgval":3.0, "sigmawithval":0., "maxres":options.maxres})
					else:
						score=float(-1-.1*np.random.rand())
					if score<best[0] : best=[score,r,k]
			
			nsym1=x0.get_nsym(options.applysym)
			if best[0]<options.simthr2: 
				avg=bxfs[best[2]]
				if nsym1==1:
					avgrs[best[1]].add_image(avg)
				else:
					for k in range(nsym1):
						x=x0.get_sym(options.applysym, k)
						c=avg.process("xform",{"transform":x})
						avgrs[best[1]].add_image(c)
					
			
			
			callback(ii*100//len(data))
			stats.append([i]+best)
			#print(fsp, i, best)
			
		
		#callback(100)
		outputs=[]
		for avg in avgrs:
			outputs.append(avg.finish())
				
				
		return (outputs, nrmvls, stats)
	
def inrange(a,b,c): return a<=b and b<=c

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_average.py <ref1> <ref2> ... [options] 
Note that this program is not part of the original e2spt hierarchy, but is part of an experimental refactoring.

Will read metadata from the specified spt_XX directory, as produced by e2spt_align.py, and average a selected subset of subtomograms in the predetermined orientation.
This version of the program competes each particle against N reference volumes, and only averages it with the best match. Alignment parameters from a previous
e2spt_align run are used to define the coarse orientation, so the references must be similar and in the same orientation. Alignments are translationally adjusted only.
If --sym is specified, each possible symmetric orientation is tested starting with the exisiting alignment parameters, and only the best is kept.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--simthr", default=-0.1,type=float,help="Similarity is smaller for better 'quality' particles. Specify the highest value to include from e2spt_hist.py. Default -0.1")
	parser.add_argument("--simthr2", default=0,type=float,help="Simlarity score for the best matching final alignment. Scaling may be different due to resolution limit. Default 0")
	parser.add_argument("--replace",type=str,default=None,help="Replace the input subtomograms used for alignment with the specified file (used when the aligned particles were masked or filtered)")
	parser.add_argument("--wedgesigma",type=float,help="Threshold for identifying missing data in Fourier space in terms of standard deviation of each Fourier shell. Default 3.0",default=3.0)
	parser.add_argument("--minalt",type=float,help="Minimum alignment altitude to include. Default=0",default=0)
	parser.add_argument("--maxalt",type=float,help="Maximum alignment altitude to include. Deafult=180",default=180)
	parser.add_argument("--maxtilt",type=float,help="Explicitly zeroes data beyond specified tilt angle. Assumes tilt axis exactly on Y and zero tilt in X-Y plane. Default 90 (no limit).",default=90.0)
	parser.add_argument("--maxres",type=float,help="Lowpass filter applied to particles prior to alignment/averaging, resolution in A. Default disabled",default=-1)
	parser.add_argument("--listfile",type=str,help="Specify a filename containing a list of integer particle numbers to include in the average, one per line, first is 0. Additional exclusions may apply.",default=None)
	parser.add_argument("--shrinkcompare",type=int,help="Shrink factor for classification only (for speed)",default=0)
	parser.add_argument("--sym",type=str,help="Symmetry of the input. Must be aligned in standard orientation to work properly. The structure will be expanded from this symmetry to c1.",default="c1")
	parser.add_argument("--applysym",type=str,help="Symmetry to apply to the structure after classification.",default="c1")
	parser.add_argument("--path",type=str,default=None,help="Path to a folder containing current results (default = highest spt_XX)")
	parser.add_argument("--parallel",type=str,default=None,help="parallel mode. Not all functions are implemented yet..")
	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.")
	parser.add_argument("--maskclass",type=str,default=None,help="Mask each reference before classification")
	parser.add_argument("--mask",type=str,default=None,help="Mask applied to final averages")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--noali", action="store_true", default=False ,help="Skip translational alignment.")
	parser.add_argument("--sample",type=int,help="use only N samples.",default=-1)
	parser.add_argument("--randnclass",type=int,help="split into N random classes. ignore refs",default=-1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()


	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="spt_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : 
			print("Error, cannot find any spt_XX folders")
			sys.exit(2)
		options.path = "spt_{:02d}".format(max(fls))
		if options.verbose : print("Working in : ",options.path)
		
	options.path=options.path.strip('/\\')

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print("Cannot find a {}/particle_parms* file".format(options.path))
			sys.exit(2)
		options.iter=max(fls)
		if options.verbose : print("Using iteration ",options.iter)
		angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))
	else:
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print("Cannot find a {}/particle_parms* file".format(options.path))
			sys.exit(2)
		mit=max(fls)
		if options.iter>mit : 
			angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,mit))
			print("WARNING: no particle_parms found for iter {}, using parms from {}".format(options.iter,mit))
		else : angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))

	if options.listfile!=None :
		plist=set([int(i) for i in open(options.listfile,"r")])

	NTHREADS=max(options.threads+1,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)


	# filter the list of particles to include 
	keys=list(angs.keys())
	if options.listfile!=None :
		keys=[i for i in keys if eval(i)[1] in plist]
		if options.verbose : print("{}/{} particles based on list file".format(len(keys),len(list(angs.keys()))))
	
	keys=[k for k in keys if angs[k]["score"]<=options.simthr and inrange(options.minalt,angs[k]["xform.align3d"].get_params("eman")["alt"],options.maxalt)]
	if options.verbose : print("{}/{} particles after filters".format(len(keys),len(list(angs.keys()))))
	
	if options.sample>0:
		rnd=np.arange(len(keys))
		np.random.shuffle(rnd)
		rnd=rnd[:options.sample]
		keys=[keys[i] for i in rnd.tolist()]
		print("using {} samples...".format(len(keys)))

	
	if options.parallel:
		if options.randnclass>0:
			nref=options.randnclass
		else:
			nref=len(args)
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer(options.parallel, module="e2spt_average_multi.SptavgmultTask")
		num_cpus = etc.cpu_est()
		
		print("{} total CPUs available".format(num_cpus))
		data=[] ## even/odd
		for i,k in enumerate(keys):
			src, ii=eval(k)[0],eval(k)[1]
			data.append([src, ii, angs[k]["xform.align3d"]])
		
		
		#### check and save size of particle
		fsp, i, xf=data[0]
		b=EMData(fsp,i, True)
		sz=options.boxsz=b["ny"]
		avgs=[]
		
		tasks=[data[i::num_cpus] for i in range(num_cpus)]
		
		
		print("{} particles in {} jobs".format(len(data), len(tasks) ))

		tids=[]
		for t in tasks:
			task = SptavgmultTask(t, args, options)
			#task.execute(print)
			#print("done")
			#exit()
			tid=etc.send_task(task)
			tids.append(tid)

		while 1:
			st_vals = etc.check_task(tids)
			#print("{:.1f}/{} finished".format(np.mean(st_vals), 100))
			#print(tids)
			if np.min(st_vals) == 100: break
			time.sleep(5)

		#dics=[0]*nptcl
		threeds=[]
		normvols=[]
		stats=[]
		print("collecting from workers...")
		for i in tids:
			threed, norm, stat=etc.get_results(i)[1]
			#print(len(threed), len(norm), threed[0]["sigma"])
			threeds.append(threed)
			normvols.append(norm)
			stats.append(stat)
			
		stats=np.vstack(stats)
		stats=stats[np.argsort(stats[:,0]),:]
		np.savetxt("{}/avg_multi_{:02d}.txt".format(options.path,options.iter), stats)
		lsts=[LSXFile(f"sets/{options.path)}_{options.iter:02d}_{i:02d}" for i in range(nref)]
		for n,score,cls,x in stats:
			lsts[cls].write(-1,n,data[0][0])
		lsts=None
		
		avs=[]
		for r in range(nref):
			output=EMData(sz, sz, sz)
			normvol=EMData((sz//2+1)*2, sz, sz)
			output.to_zero()
			output.do_fft_inplace()
			
			avg=Averagers.get("mean")
			normvol.to_zero()
			
			for k,thd in enumerate(threeds):
				threed=thd[r]
				norm=normvols[k][r]
				threed.process_inplace("math.multamplitude", {"amp":norm})
				avg.add_image(threed)
				normvol.add(norm)
				
			output=avg.finish()
			normvol.process_inplace("math.reciprocal")
			output.process_inplace("math.multamplitude", {"amp":normvol})
			output.process_inplace("xform.phaseorigin.tocenter")
			
			output.do_ift_inplace()
			output.depad()
			avs.append(output)
	else:
		
		n=len(args)
		args=[comma(i) for i in args]
		refs=[EMData(i[0],i[1]) for i in args]
		if options.maxres>0:
			for r in refs: r.process_inplace("filter.lowpass.gauss",{"cutoff_freq":old_div(1.0,options.maxres)})
		if options.maskclass!=None:
			mask=EMData(options.maskclass)
			for r in refs: r.mult(mask)
		
		jsd=queue.Queue(0)

		avgs=[Averagers.get("mean.tomo",{"thresh_sigma":options.wedgesigma}) for i in range(n)]
		
		if options.shrinkcompare>1 :
			shrinkrefs=[r.process("math.meanshrink",{"n":options.shrinkcompare}) for r in refs]
		else:
			shrinkrefs=refs

		# Rotation and insertion are slow, so we do it with threads. 
		# Averager isn't strictly threadsafe, so possibility of slight numerical errors with a lot of threads
		if options.replace != None:
			thrds=[threading.Thread(target=rotfncompete,args=(jsd,avgs,options.replace,eval(k)[1],angs[k]["xform.align3d"],options.sym,refs,shrinkrefs,options.maxtilt,options.wedgesigma,options.shrinkcompare,options.maxres,options.simthr2,options.verbose)) for i,k in enumerate(keys)]

		else:
			#thrds=[threading.Thread(target=rotfncompete,args=(jsd,avgs,eval(k)[0],eval(k)[1],angs[k]["xform.align3d"],options.sym,refs,shrinkrefs,options.maxtilt,options.wedgesigma,options.shrinkcompare,options.maxres,options.simthr2,options.verbose)) for i,k in enumerate(keys)]
			thrds=[threading.Thread(target=rotfncompete,args=(jsd,avgs,eval(k)[0],eval(k)[1],angs[k]["xform.align3d"],refs,shrinkrefs,options)) for i,k in enumerate(keys)]


		print(len(thrds)," threads")
		thrtolaunch=0
		out=open("{}/avg_multi_{:02d}.txt".format(options.path,options.iter),"w")
		while thrtolaunch<len(thrds) or threading.active_count()>1 or not jsd.empty():
			# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
			# note that it's ok that we wait here forever, since there can't be new results if an existing
			# thread hasn't finished.
			if thrtolaunch<len(thrds) :
				while (threading.active_count()==NTHREADS ) : time.sleep(.1)
				if options.verbose : print("Starting thread {}/{}".format(thrtolaunch,len(thrds)))
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else: time.sleep(1)
		
			while not jsd.empty():
				fspn,score,ref,sym=jsd.get()
				out.write("{}\t{}\t{}\t{}\n".format(fspn,score,ref,sym))		# Output columns are img #, best score, # of best ref, # of best sym
				#avg[n%2].add_image(ptcl)
			out.flush()

		for t in thrds:
			t.join()

		avs=[i.finish() for i in avgs]

	if options.mask: mask=EMData(options.mask)
	for i,v in enumerate(avs):
		if options.mask: v.mult(mask)
		v.write_image("{}/threed_{:02d}_{:02d}.hdf".format(options.path,options.iter,i),0)

	print("Done")
	E2end(logid)

def comma(fsp):
	if "," in fsp: return fsp.split(",")[0],int(fsp.split(",")[-1])
	return fsp,0

if __name__ == "__main__":
	main()

