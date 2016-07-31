#!/usr/bin/env python
# average selected subset of particles

from EMAN2 import *
import time
import os
import threading
import Queue
from sys import argv,exit

def rotfn(avg,fsp,i,a,verbose):
	b=EMData(fsp,i)
	b.process_inplace("xform",{"transform":a})
	avg.add_image(b)
	#jsd.put((fsp,i,b))

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_average.py [options] 
Note that this program is not part of the original e2spt hierarchy, but is part of an experimental refactoring.

Will read metadata from the specified spt_XX directory, as produced by e2spt_align.py, and average a selected subset of subtomograms in the predetermined orientation.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--simthr", default=-0.1,type=float,help="Similarity is smaller for better 'quality' particles. Specify the highest value to include from e2spt_hist.py. Default -0.1")
	parser.add_argument("--replace",type=str,default=None,help="Replace the input subtomograms used for alignment with the specified file (used when the aligned particles were masked or filtered)")
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

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print "Cannot find a {}/particle_parms* file".format(options.path)
			sys.exit(2)
		options.iter=max(fls)

	NTHREADS=max(options.threads+1,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)

	angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))
#	jsd=Queue.Queue(0)


	avg=[0,0]
	avg[0]=Averagers.get("mean.tomo") #,{"save_norm":1})
	avg[1]=Averagers.get("mean.tomo")

	# Rotation and insertion are slow, so we do it with threads. 
	# Averager isn't strictly threadsafe, so possibility of slight numerical errors with a lot of threads
	if options.replace != None:
		thrds=[threading.Thread(target=rotfn,args=(avg[i%2],options.replace,eval(k)[1],angs[k]["xform.align3d"],options.verbose)) for i,k in enumerate(angs.keys()) if angs[k]["score"]<=options.simthr]
	else:
		thrds=[threading.Thread(target=rotfn,args=(avg[i%2],eval(k)[0],eval(k)[1],angs[k]["xform.align3d"],options.verbose)) for i,k in enumerate(angs.keys()) if angs[k]["score"]<=options.simthr]


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
	
		#while not jsd.empty():
			#fsp,n,ptcl=jsd.get()
			#avg[n%2].add_image(ptcl)


	for t in thrds:
		t.join()

	ave=avg[0].finish()		#.process("xform.phaseorigin.tocenter").do_ift()
	avo=avg[1].finish()		#.process("xform.phaseorigin.tocenter").do_ift()
	av=ave+avo
	av.mult(0.5)

	evenfile="{}/threed_{:02d}_even.hdf".format(options.path,options.iter)
	oddfile="{}/threed_{:02d}_odd.hdf".format(options.path,options.iter)
	combfile="{}/threed_{:02d}.hdf".format(options.path,options.iter)
	ave.write_image(evenfile,0)
	avo.write_image(oddfile,0)
	av.write_image(combfile,0)

	cmd="e2proc3d.py {evenfile} {path}/fsc_unmasked_{itr:02d}.txt --calcfsc={oddfile}".format(path=options.path,itr=options.iter,evenfile=evenfile,oddfile=oddfile)
	launch_childprocess(cmd)

	launch_childprocess("e2proc3d.py {combfile} {combfile} --process=filter.wiener.byfsc:fscfile={path}/fsc_unmasked_{itr:02d}.txt:snrmult=2".format(path=options.path,itr=options.iter,combfile=combfile))


	E2end(logid)


if __name__ == "__main__":
	main()

