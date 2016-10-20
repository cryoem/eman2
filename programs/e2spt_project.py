#!/usr/bin/env python
# average selected subset of particles

from EMAN2 import *
import time
import os
import threading
import Queue
from sys import argv,exit

def rotfn(jsd,fsp,i,a,verbose):
	b=EMData(fsp,i)
	b.process_inplace("xform",{"transform":a})
	c=b.process("misc.directional_sum",{"axis":"z"})
	c["source_orig_path"]=fsp
	c["source_orig_n"]=i
	jsd.put((i,c))


def inrange(a,b,c): return a<=b and b<=c

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_average.py [options] 
Note that this program is not part of the original e2spt hierarchy, but is part of an experimental refactoring.

Will read metadata from the specified spt_XX directory, as produced by e2spt_align.py, and generate projections of the selected subset of particles after
alignment to the reference, projected along the Z axis.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--simthr", default=-0.1,type=float,help="Similarity is smaller for better 'quality' particles. Specify the highest value to include from e2spt_hist.py. Default -0.1")
	parser.add_argument("--replace",type=str,default=None,help="Replace the input subtomograms used for alignment with the specified file (used when the aligned particles were masked or filtered)")
	parser.add_argument("--minalt",type=float,help="Minimum alignment altitude to include. Default=0",default=0)
	parser.add_argument("--maxalt",type=float,help="Maximum alignment altitude to include. Deafult=180",default=180)
	parser.add_argument("--symalimasked",type=str,default=None,help="This will translationally realign each asymmetric unit to the specified (usually masked) reference ")
	parser.add_argument("--sym",type=str,default=None,help="Symmetry of the input. Must be aligned in standard orientation to work properly.")
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

	NTHREADS=max(options.threads+1,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)

	jsd=Queue.Queue(0)


	# Rotation and insertion are slow, so we do it with threads. 
	# Averager isn't strictly threadsafe, so possibility of slight numerical errors with a lot of threads
	if options.replace != None:
		thrds=[threading.Thread(target=rotfn,args=(jsd,options.replace,eval(k)[1],angs[k]["xform.align3d"],options.verbose)) for i,k in enumerate(angs.keys()) if angs[k]["score"]<=options.simthr and inrange(options.minalt,angs[k]["xform.align3d"].get_params("eman")["alt"],options.maxalt)]

	else:
		thrds=[threading.Thread(target=rotfn,args=(jsd,eval(k)[0],eval(k)[1],angs[k]["xform.align3d"],options.verbose)) for i,k in enumerate(angs.keys()) if angs[k]["score"]<=options.simthr and inrange(options.minalt,angs[k]["xform.align3d"].get_params("eman")["alt"],options.maxalt) ]


	print len(thrds)," threads"
	thrtolaunch=0
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
			n,ptcl=jsd.get()
			ptcl.write_image("{}/projections_{:02d}.hdf".format(options.path,options.iter),-1)
#			avg[n%2].add_image(ptcl)


	for t in thrds:
		t.join()

	E2end(logid)


if __name__ == "__main__":
	main()

