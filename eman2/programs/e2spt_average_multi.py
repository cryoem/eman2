#!/usr/bin/env python
# average selected subset of particles

from EMAN2 import *
import time
import os
import threading
import Queue
from sys import argv,exit

def rotfncompete(jsd,avgs,fsp,fspn,a,sym,refs,shrinkrefs,maxtilt,wedgesigma,shrink,maxres,simthr2,verbose):
	"""Averaging thread. 
	avgs are n existing Averagers, 
	fsp,i is the particle being averaged
	a is the Transform to put the particle in the correct orientation
	sym is the symmetry for replication of the particle
	refs are the n alignment references for competitive averaging
	maxtilt can optionally better enforce missing wedge exclusion
	"""
	if shrink<2: shrink=0
	b=EMData(fsp,fspn)
	if maxres>0: b.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/maxres})
	if maxtilt<90.0 :
		bf=b.do_fft()
		bf.process_inplace("mask.wedgefill",{"thresh_sigma":0.0,"maxtilt":maxtilt})
		b=bf.do_ift()
	b.process_inplace("xform",{"transform":a})
	if shrink : bs=b.process("math.meanshrink",{"n":shrink})
	else: bs=b
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(sym)
	best=(1.0e50,None,None)
	for r,ref in enumerate(shrinkrefs):
		for i in xrange(nsym):
			c=bs.process("xform",{"transform":xf.get_sym(sym,i)})
			d=c.align("translational",ref)
			score=d.cmp("fsc.tomo.auto",ref,{"sigmaimgval":wedgesigma,"sigmawithval":0.5})
			if score<best[0] : best=(score,r,d,i)
	
	if shrink:
		c=b.process("xform",{"transform":xf.get_sym(sym,best[3])})
		d=c.align("translational",refs[best[1]])
	else :
		d=best[2]
	
	if best[0]<simthr2 : 
		avgs[best[1]].add_image(d)
		print "{} -> ref {} sym {}   {}".format(fspn,best[1],best[3],best[0])
	else: print "** {} -> ref {} sym {}   {}".format(fspn,best[1],best[3],best[0])
	jsd.put((fspn,best[0],best[1],best[3]))


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

	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.")
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
	parser.add_argument("--sym",type=str,help="Symmetry of the input. Must be aligned in standard orientation to work properly.",default="c1")
	parser.add_argument("--path",type=str,default=None,help="Path to a folder containing current results (default = highest spt_XX)")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="spt_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : 
			print "Error, cannot find any spt_XX folders"
			sys.exit(2)
		options.path = "spt_{:02d}".format(max(fls))
		if options.verbose : print "Working in : ",options.path

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print "Cannot find a {}/particle_parms* file".format(options.path)
			sys.exit(2)
		options.iter=max(fls)
		if options.verbose : print "Using iteration ",options.iter
		angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))
	else:
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print "Cannot find a {}/particle_parms* file".format(options.path)
			sys.exit(2)
		mit=max(fls)
		if options.iter>mit : 
			angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,mit))
			print "WARNING: no particle_parms found for iter {}, using parms from {}".format(options.iter,mit)
		else : angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))

	n=len(args)
	refs=[EMData(i) for i in args]
	if options.maxres>0:
		for r in refs: r.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/options.maxres})
	
	if options.listfile!=None :
		plist=set([int(i) for i in file(options.listfile,"r")])

	NTHREADS=max(options.threads+1,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)

	jsd=Queue.Queue(0)

	avgs=[Averagers.get("mean.tomo",{"thresh_sigma":options.wedgesigma}) for i in xrange(n)]

	# filter the list of particles to include 
	keys=angs.keys()
	if options.listfile!=None :
		keys=[i for i in keys if eval(i)[1] in plist]
		if options.verbose : print "{}/{} particles based on list file".format(len(keys),len(angs.keys()))
	
	keys=[k for k in keys if angs[k]["score"]<=options.simthr and inrange(options.minalt,angs[k]["xform.align3d"].get_params("eman")["alt"],options.maxalt)]
	if options.verbose : print "{}/{} particles after filters".format(len(keys),len(angs.keys()))
																		 

	if options.shrinkcompare>1 :
		shrinkrefs=[r.process("math.meanshrink",{"n":options.shrinkcompare}) for r in refs]
	else:
		shrinkrefs=refs

	# Rotation and insertion are slow, so we do it with threads. 
	# Averager isn't strictly threadsafe, so possibility of slight numerical errors with a lot of threads
	if options.replace != None:
		thrds=[threading.Thread(target=rotfncompete,args=(jsd,avgs,options.replace,eval(k)[1],angs[k]["xform.align3d"],options.sym,refs,shrinkrefs,options.maxtilt,options.wedgesigma,options.shrinkcompare,options.maxres,options.simthr2,options.verbose)) for i,k in enumerate(keys)]

	else:
		thrds=[threading.Thread(target=rotfncompete,args=(jsd,avgs,eval(k)[0],eval(k)[1],angs[k]["xform.align3d"],options.sym,refs,shrinkrefs,options.maxtilt,options.wedgesigma,options.shrinkcompare,options.maxres,options.simthr2,options.verbose)) for i,k in enumerate(keys)]


	print len(thrds)," threads"
	thrtolaunch=0
	out=file("{}/avg_multi_{:02d}.txt".format(options.path,options.iter),"w")
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			if options.verbose : print "Starting thread {}/{}".format(thrtolaunch,len(thrds))
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
	
	for i,v in enumerate(avs):
		v.write_image("{}/threed_{:02d}_{:02d}.hdf".format(options.path,options.iter,i),0)


	E2end(logid)


if __name__ == "__main__":
	main()

