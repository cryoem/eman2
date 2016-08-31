#!/usr/bin/env python
# average selected subset of particles

from EMAN2 import *
import time
import os
import threading
import Queue
from sys import argv,exit

def rotfn(avg,fsp,i,a,maxtilt,verbose):
	b=EMData(fsp,i)
	if maxtilt<90.0 :
		bf=b.do_fft()
		bf.process_inplace("mask.wedgefill",{"thresh_sigma":0.0,"maxtilt":maxtilt})
		b=bf.do_ift()
	b.process_inplace("xform",{"transform":a})
	avg.add_image(b)
	#jsd.put((fsp,i,b))

def rotfnsym(avg,fsp,i,a,sym,masked,maxtilt,verbose):
	b=EMData(fsp,i)
	if maxtilt<90.0 :
		bf=b.do_fft()
		bf.process_inplace("mask.wedgefill",{"thresh_sigma":0.0,"maxtilt":maxtilt})
		b=bf.do_ift()
	b.process_inplace("xform",{"transform":a})
	xf = Transform()
	xf.to_identity()
	nsym=xf.get_nsym(sym)
	for i in xrange(nsym):
		c=b.process("xform",{"transform":xf.get_sym(sym,i)})
		d=c.align("translational",masked)
		avg.add_image(d)
	#jsd.put((fsp,i,b))

def inrange(a,b,c): return a<=b and b<=c

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
	parser.add_argument("--wedgesigma",type=float,help="Threshold for identifying missing data in Fourier space in terms of standard deviation of each Fourier shell. Default 3.0",default=3.0)
	parser.add_argument("--minalt",type=float,help="Minimum alignment altitude to include. Default=0",default=0)
	parser.add_argument("--maxalt",type=float,help="Maximum alignment altitude to include. Deafult=180",default=180)
	parser.add_argument("--maxtilt",type=float,help="Explicitly zeroes data beyond specified tilt angle. Assumes tilt axis exactly on Y and zero tilt in X-Y plane. Default 90 (no limit).",default=90.0)
	parser.add_argument("--listfile",type=str,help="Specify a filename containing a list of integer particle numbers to include in the average, one per line, first is 0. Additional exclusions may apply.",default=None)
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

	if options.listfile!=None :
		plist=set([int(i) for i in file(options.listfile,"r")])

	NTHREADS=max(options.threads+1,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)

#	jsd=Queue.Queue(0)


	avg=[0,0]
	avg[0]=Averagers.get("mean.tomo",{"thresh_sigma":options.wedgesigma}) #,{"save_norm":1})
	avg[1]=Averagers.get("mean.tomo",{"thresh_sigma":options.wedgesigma})

	# filter the list of particles to include 
	keys=angs.keys()
	if options.listfile!=None :
		keys=[i for i in keys if eval(i)[1] in plist]
		if options.verbose : print "{}/{} particles based on list file".format(len(keys),len(angs.keys()))
	
	keys=[k for k in keys if angs[k]["score"]<=options.simthr and inrange(options.minalt,angs[k]["xform.align3d"].get_params("eman")["alt"],options.maxalt)]
	if options.verbose : print "{}/{} particles after filters".format(len(keys),len(angs.keys()))
																		 

	# Rotation and insertion are slow, so we do it with threads. 
	if options.symalimasked!=None:
		if options.replace!=None :
			print "Error: --replace cannot be used with --symalimasked"
			sys.exit(1)
		alimask=EMData(options.symalimasked)
		thrds=[threading.Thread(target=rotfnsym,args=(avg[i%2],eval(k)[0],eval(k)[1],angs[k]["xform.align3d"],options.sym,alimask,options.maxtilt,options.verbose)) for i,k in enumerate(keys)]
	else:
		# Averager isn't strictly threadsafe, so possibility of slight numerical errors with a lot of threads
		if options.replace != None:
			thrds=[threading.Thread(target=rotfn,args=(avg[i%2],options.replace,eval(k)[1],angs[k]["xform.align3d"],options.maxtilt,options.verbose)) for i,k in enumerate(keys)]

		else:
			thrds=[threading.Thread(target=rotfn,args=(avg[i%2],eval(k)[0],eval(k)[1],angs[k]["xform.align3d"],options.maxtilt,options.verbose)) for i,k in enumerate(keys)]


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
	# impose symmetry on even and odd halves if appropriate
	if options.sym!=None and options.sym.lower()!="c1" and options.symalimasked==None:
		ave.process_inplace("xform.applysym",{"averager":"mean.tomo","sym":options.sym})
		avo.process_inplace("xform.applysym",{"averager":"mean.tomo","sym":options.sym})
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

