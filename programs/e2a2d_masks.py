#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
# average particles together based on existing alignments

from past.utils import old_div
from EMAN2 import *
import time
import os
#import threading
#import Queue
from sys import argv,exit
from e2classaverage import class_average

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2a2d_masks.py <mask1> <mask2> ... [options] 
	This program takes the results of e2a2d_align and computes statistics under one or more masks on the aligned particles."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_argument("--threads", default=4,type=int,help="Number of alignment threads to run in parallel on a single computer", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	#parser.add_argument("--goldstandard",type=float,help="If specified, will phase randomize the even and odd references past the specified resolution (in A, not 1/A)",default=0)
	#parser.add_argument("--goldcontinue",action="store_true",help="Will use even/odd refs corresponding to specified reference to continue refining without phase randomizing again",default=False)
	parser.add_argument("--path",type=str,default=None,help="Path to a folder with existing e2a2d_align results (default = 2da_XX)")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = highest existing iteration",default=0)
	parser.add_argument("--extract",action="store_true",help="If set, will also produce a .txt file for plotting with all of the per-particle statistics",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=1, help="verbose level [0-9], higher number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.path == None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="m2d_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : 
			print("Error: No m2d_* folders found")
			sys.exit(1)
		options.path = "m2d_{:02d}".format(max(fls))
		if options.verbose : print("Working in folder: ",options.path)

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print("Error: No particle_parms* files found in ",options.path)
			sys.exit(2)
		else: options.iter=max(fls)
		if options.verbose: print("Iteration: ",options.iter)

#	NTHREADS=max(options.threads,2)		# we have one thread just writing results

	logid=E2init(sys.argv, options.ppid)

	angs=js_open_dict("{}/particle_parms_{:02d}.json".format(options.path,options.iter))

	sortangs=[(v["score"],v["xform.align2d"],eval(k),k) for k,v in list(angs.items())]		# the eval() here is poor programming practice, but is the easiest way :^(
	N=len(sortangs)
		
	#if options.scorebestset>0:
		#out=LSXFile("sets/{}_{}_{}.lst".format(options.path.split("/")[-1],options.iter,options.scorebestset))
		#for i,t in enumerate(sorted(sortangs)[:options.scorebestset]):
			#imh=EMData(t[2][0],t[2][1],True)
			#try:
				#fsp=imh["data_source"]	# if the data is from a LST file this dereferences it
				#fspn=imh["data_n"]
			#except:
				#fsp=imh["source_path"]
				#fspn=imh["source_n"]
			#out.write(i,fspn,fsp)
			
	masks=[EMData(f,0) for f in args]
	names=[os.path.split(f)[1].rsplit(".",1)[0] for f in args]
	
	# open the extract file and write the description line
	if options.extract:
		outtxt=open(f"{options.path}/maskstats_{options.iter}.txt","w")
		outtxt.write("# score")
		for n in names: outtxt.write(f"; {n}_mean; {n}_sig")
		outtxt.write("\n")
	
	# the main per particle loop
	t0=time.time()
	t1=t0
	for i,t in enumerate(sorted(sortangs,key=lambda x:x[0])):
		if options.verbose==1 and time.time()-t1>1:
			t1=time.time()
			frac=old_div(i,float(N))
			try:
				remain=int(old_div((time.time()-t0),frac)-(time.time()-t0))	# est remaining time in sec
				print("{:6d}/{:-6d}   time remaining: {}:{:02d}     \r".format(i,N,remain//60,remain%60))
			except:
				print("{:6d}/{:-6d}     \r".format(i,N))
				
			sys.stdout.flush()
			
		# read and transform the particle
		im=EMData(t[2][0],t[2][1])			# this is from the key of angs{}
		im.process_inplace("xform",{"transform":t[1]})
		adic=angs[t[3]]
		
		if options.extract : line=[t[0]]
		
		# i is particle #, j is mask number
		for j in len(masks):
			tmp=im*masks[j]
			adic[names[j]+"_mean"]=tmp["mean"]
			adic[names[j]+"_sigma"]=tmp["sigma_nonzero"]	# we want to ignore regions which are missing due to the transform
			if options.extract : 
				line.append(tmp["mean"])
				line.append(tmp["sigma_nonzeo"])
				
		angs.setval(t[3],adic,deferupdate=True)

		if options.extract:
			for l in line: outtxt.write(f"{l}\t")
			outtxt.write(f"# {t[2][1]};{t[2][0]}\n")
		

	angs.update()	# since we deferred above
	
	if options.verbose : print("\nDone!")

	E2end(logid)


if __name__ == "__main__":
	main()

