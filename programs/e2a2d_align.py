#!/usr/bin/env python
# align all particles to reference and store alignment results

from EMAN2 import *
import time
import os
import threading
import Queue
from sys import argv,exit

def ali2dfn(jsd,fsp,il,a,options):
#	t=time.time()
	global frac
	rslt=[]
	for n,i in enumerate(il):
		b=EMData(fsp,i)
		c=b.align(options.align[0],a,options.align[1],options.aligncmp[0],options.aligncmp[1])
		if options.ralign!=None and options.ralign[0]!=None:
			rparms=options.ralign[1]
			rparms["xform.align2d"]=c["xform.align2d"]
			c=b.align(options.ralign[0],a,rparms,options.raligncmp[0],options.raligncmp[1])
			sim=c.cmp(options.cmp[0],a,options.cmp[1])
			rslt.append((i,{"xform.align2d":c["xform.align2d"],"prealign":rparms["xform.align2d"],"score":sim}))
		else: 
			sim=c.cmp(options.cmp[0],a,options.cmp[1])
			rslt.append((i,{"xform.align2d":c["xform.align2d"],"score":sim}))
		if il[0]==0 : frac=n/float(len(il))		# we monitor progress in the first thread

	jsd.put((fsp,rslt))
#	if options.verbose>1 : print "{}\t{}\t{}\t{}".format(fsp,i,time.time()-t,c[0]["score"])

frac=0

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e22da_align.py [options] <sets/file.lst> <reference>[,nref]
	This program will align a stack of 2-D particle images to a single 2-D reference image, recording all of the alignment parameters in a
	JSON file. This is designed as the first step in a longer 2-D processing sequence. 
	
	This could, for example, be used to align a full particle set against a 2-D reference, then select a subset which matches well, as opposed 
	to the normal 3-D refinement process where particles are competed against many references."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	#parser.add_argument("--goldstandard",type=float,help="If specified, will phase randomize the even and odd references past the specified resolution (in A, not 1/A)",default=0)
	#parser.add_argument("--goldcontinue",action="store_true",help="Will use even/odd refs corresponding to specified reference to continue refining without phase randomizing again",default=False)
	parser.add_argument("--saveali",action="store_true",help="Save a stack file (aliptcls_XX.hdf) containing the aligned particles.",default=False)
	parser.add_argument("--path",type=str,default=None,help="Path to a folder where results should be stored, following standard naming conventions (default = 2da_XX)")
	#parser.add_argument("--sym",type=str,default="c1",help="Symmetry of the input. Must be aligned in standard orientation to work properly.")
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average. Default is None.", default="rotate_translate_tree")
	parser.add_argument("--aligncmp",type=str,help="The comparitor used for the --align aligner. Default is dot.",default="ccc")
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	parser.add_argument("--raligncmp",type=str,help="The comparitor used by the second stage aligner.",default="ccc")
	parser.add_argument("--cmp",type=str,help="The comparitor used to generate quality scores for the purpose of particle exclusion in classes, strongly linked to the keep argument.", default="ccc")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.ralign=parsemodopt(options.ralign)
	options.raligncmp=parsemodopt(options.raligncmp)
	options.cmp=parsemodopt(options.cmp)


	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="m2d_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : fls=[0]
		options.path = "m2d_{:02d}".format(max(fls)+1)
		try: os.mkdir(options.path)
		except: pass
		if options.verbose : print "Working in folder: ",options.path

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : options.iter=1
		else: options.iter=max(fls)+1
		if options.verbose: print "Iteration: ",options.iter

	reffile=args[1]
	if "," in reffile:
		reffile,nref=reffile.split(",")
		nref=int(nref)
	else: nref=0
	
	ref=EMData(reffile,nref)
	
	NTHREADS=max(options.threads,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)

	angs={}
	jsd=Queue.Queue(0)

	n=-1
	N=EMUtil.get_image_count(args[0])
	thrds=[threading.Thread(target=ali2dfn,args=(jsd,args[0],xrange(i,N,NTHREADS-1),ref,options)) for i in xrange(NTHREADS-1)]

	# here we run the threads and save the results, no actual alignment done here
	if options.verbose: print len(thrds)," threads"
	thrtolaunch=0
	while thrtolaunch<len(thrds) or threading.active_count()>1:
		# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
		# note that it's ok that we wait here forever, since there can't be new results if an existing
		# thread hasn't finished.
		if thrtolaunch<len(thrds) :
			while (threading.active_count()==NTHREADS ) : time.sleep(.1)
			if options.verbose>1 : print "Starting thread {}/{}".format(thrtolaunch,len(thrds))
			thrds[thrtolaunch].start()
			thrtolaunch+=1
		else: time.sleep(1)
		if options.verbose>1: print "{}% complete".format(100.0*frac)
	
		while not jsd.empty():
			fsp,nds=jsd.get()
			for n,d in nds:
				angs[(fsp,n)]=d
				if options.saveali:
					v=EMData(fsp,n)
					v.transform(d["xform.align2d"])
					v.write_image("{}/aliptcls_{:02d}.hdf".format(options.path,options.iter),n)


	for t in thrds:
		t.join()

	if options.verbose : print "Writing results"

	angsd=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))
	angsd.update(angs)

	if options.verbose : print "Done!"

	E2end(logid)


if __name__ == "__main__":
	main()

