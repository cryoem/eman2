#!/usr/bin/env python
# average particles together based on existing alignments

from EMAN2 import *
import time
import os
#import threading
#import Queue
from sys import argv,exit

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2a2d_average.py [options] 
	This program takes the results of e2a2d_align and produces averages of subsets of the aligned particles."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	#parser.add_argument("--goldstandard",type=float,help="If specified, will phase randomize the even and odd references past the specified resolution (in A, not 1/A)",default=0)
	#parser.add_argument("--goldcontinue",action="store_true",help="Will use even/odd refs corresponding to specified reference to continue refining without phase randomizing again",default=False)
	parser.add_argument("--saveali",action="store_true",help="Save a stack file (aliptcls_XX.hdf) containing the aligned particles.",default=False)
	parser.add_argument("--path",type=str,default=None,help="Path to a folder with existing e2a2d_align results (default = 2da_XX)")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = highest existing iteration",default=0)
	parser.add_argument("--scrorebands", default=0,type=int,help="If specified will generate averages over N bands of 'score' values, including only particles in each band.")
	parser.add_argument("--scoreprogressive", default=0,type=int,help="If specified will generate progressive averages over N bands of 'score' values, including all particles starting with the best through the progressive bands.")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="m2d_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : 
			print "Error: No m2d_* folders found"
			sys.exit(1)
		options.path = "m2d_{:02d}".format(max(fls))
		if options.verbose : print "Working in folder: ",options.path

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print "Error: No particle_parms* files found in ",options.path
			sys.exit(2)
		else: options.iter=max(fls)
		if options.verbose: print "Iteration: ",options.iter

#	NTHREADS=max(options.threads,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)

	angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))

	sortangs=[(v["score"],v["xform.align2d"],k) for k,v in angs.items()]
	N=len(sortangs)
	
	t0=time.time()
	for i,t in enumerate(sorted(sortangs)):
		if options.verbose==1 and time.time()-t0>1:
			t0=time.time()
			print "{:6d}/{:-6d}".format(i,N)
			
		im=EMData(t[2][0],t[2][1])			# this is from the key of angs{}
		if options.scoreprogressive>0 :
			try: sp.add(im)
			except: sp=im
		
			if i%(N//options.scoreprogressive)==0 :
				sp["ptcl_repr"]=i
				sp.process("normalize.edgemean").write_image("{}/avg_prog_{:02d}.hdf".format(options.path,options.iter),-1)

		if options.scorebands>0 :
			try: sb.add(im)
			except: sb=im
		
			if (i+1)%(N//options.scorebands)==0 :
				sb["ptcl_repr"]=N//options.scorebands
				sb.process("normalize.edgemean").write_image("{}/avg_band_{:02d}.hdf".format(options.path,options.iter),-1)
				sb=None
				

	if options.verbose : print "Done!"

	E2end(logid)


if __name__ == "__main__":
	main()

