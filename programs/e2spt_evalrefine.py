#!/usr/bin/env python
#
# Author: Steven Ludtke, 12/21/2020
# Copyright (c) 2010 Baylor College of Medicine
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
from EMAN2 import *
from math import *
import os
import sys
import datetime
from numpy import array
import traceback
import json
from time import time,sleep

try:
	import numpy as np
	import matplotlib
	matplotlib.use("AGG")
#	matplotlib.use("PDF")
	import matplotlib.pyplot as plt
	pltcolors=["k","b","g","r","m","c","darkblue","darkgreen","darkred","darkmagenta","darkcyan","0.5"]
except:
	print("Matplotlib not available, plotting options will not be available")


def main():
	global classmx,nptcl,cmxtx,cmxty,cmxalpha,cmxmirror,eulers,threed,ptclmask,rings,pf,cptcl
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]
	This program performs various assessments of e2spt_refine (or similar) runs, and operates in one of several possible modes:

	--jsonortcmp
		compares the orientation of each particle between two different json files (may be from the same or different runs)
		producing jsonstat.txt

	--timing
		will report how long each refinement took to complete as well as individual tasks within the refinement
		
	--timingbypath
		will report total timing information for each refinement folder, along with useful refinement parameters
	
	--resolution_all
		Computes FSC curves for the final iteration of every refinement folder
		
	--resolution_vsref
		Computes FSC curve for the final iteration of each refine_xx folder vs a provided reference volume. 
		
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--jsonortcmp", default=False, action="store_true", help="Compare the particle orientations from two .json files. Provide the path to 2 json files as arguments to the command.")
	parser.add_argument("--timing", default=False, action="store_true", help="Report on the time required for each step of each refinement run")
	parser.add_argument("--timingbypath", default=False, action="store_true", help="Report on the CPU time required in each refine_xx folder")
	parser.add_argument("--resolution_all", default=False, action="store_true", help="generates resolution plot with the last iteration of all refine_xx directories")
	parser.add_argument("--resolution_vsref", type=str, default=None, help="Computes the FSC between the last iteration of each refine_xx directory and a specified reference map. Map must be aligned, but will be rescaled if necessary.")
	#parser.add_argument("--iter", type=int, default=None, help="If a refine_XX folder is being used, this selects a particular refinement iteration. Otherwise the last complete iteration is used.")
	#parser.add_argument("--mask",type=str,help="Mask to be used to focus --evalptclqual and other options. May be useful for separating heterogeneous data.", default=None)
	parser.add_argument("--sym",type=str,help="Symmetry to be used in searching adjacent unit cells", default="c1")
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful",guitype='intbox', row=9, col=0, rowspan=1, colspan=1, mode='evalptcl[4]')
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()

	xticlocs=[i for i in (.01,.05,.0833,.125,.1667,.2,.25,.3333,.4,.5)]
	xticlbl=["1/100","1/20","1/12","1/8","1/6","1/5","1/4","1/3","1/2.5","1/2"]
	yticlocs=(0.0,.125,.143,.25,.375,.5,.625,.75,.875,1.0)
	yticlbl=("0"," ","0.143","0.25"," ","0.5"," ","0.75"," ","1.0")
	yticlocs2=(0.0,.125,.25,.375,.5,.625,.75,.875,1.0)
	yticlbl2=("0"," ","0.25"," ","0.5"," ","0.75"," ","1.0")


	logid=E2init(sys.argv,options.ppid)

	
	if options.jsonortcmp:
		if len(args)!=2:
			print("Please specify the two json files to compare, eg - e2spt_evalrefine.py --jsonortcmp <json1> <json2>")
			sys.exit(1)
			
		js1=js_open_dict(args[0])
		js2=js_open_dict(args[1])
		
		out=open("jsonstat.txt","w")
		out.write(f"# {args[0]} vs {args[1]}\n")
		out.write(f"# score1, score 2, angle diff, sym angle diff, angle axis X, Y, Z, az1, alt1, phi1, az2, alt2, phi2   # particle ref\n")
		
		for k in js1.keys():
			o1=js1[k]["xform.align3d"]
			try: o2=js2[k]["xform.align3d"]
			except: continue
			d=o1*(o2.inverse())
			dspin=d.get_rotation("spin")	# effectively quaternion
			o1e=o1.get_rotation("eman")
			o2e=o2.get_rotation("eman")
			ke=eval(k)
			dsym=angle_ab_sym(options.sym,o1.inverse(),o2.inverse())
			
			out.write(f'{js1[k]["score"]:1.4f}\t{js2[k]["score"]:1.4f}\t{dspin["omega"]:1.4f}\t{dsym:1.4f}\t{dspin["n1"]:1.4f}\t{dspin["n2"]:1.4f}\t{dspin["n3"]:1.4f}\t{o1e["az"]:1.4f}\t{o1e["alt"]:1.4f}\t{o1e["phi"]:1.4f}\t{o2e["az"]:1.4f}\t{o2e["alt"]:1.4f}\t{o2e["phi"]:1.4f}\t# {ke[0]},{ke[1]}\n')
		
	if options.resolution_all:
		######################
		### Resolution plot
		### broken up into multiple try/except blocks because we need some of the info, even if plotting fails
		plt.title("Gold Standard Resolution")
		plt.xlabel(r"Spatial Frequency (1/$\AA$)")
		plt.ylabel("FSC")

		refines=[i for i in os.listdir(".") if "refine_" in i or "spt_" in i or "subtlt_" in i]
		fscs=[]
		for r in refines:
			try: itr=max([i for i in os.listdir(r) if "fsc_masked" in i and i[-4:]==".txt"])
			except: continue
			fscs.append("{}/{}".format(r,itr))

		fscs.sort(reverse=True)
		maxx=0.01

		# iterate over fsc curves
		for num,f in enumerate(fscs):
			# read the fsc curve
			d=np.loadtxt(f).transpose()

			# plot the curve
			try: plt.plot(d[0],d[1],label=f[:9],color=pltcolors[(num)%12])
			except: pass
			maxx=max(maxx,max(d[0]))

		if max(d[0])<max(xticlocs) :
			xticlocs=[i for i in xticlocs if i<=max(d[0])]
			xticlbl=xticlbl[:len(xticlocs)]
		plt.axhline(0.0,color="k")
		plt.axhline(0.143,color="#306030",linestyle=":")
		plt.axis((0,maxx,-.02,1.02))
		plt.legend(loc="upper right",fontsize="x-small")
		plt.xticks(xticlocs,xticlbl)
		plt.yticks(yticlocs,yticlbl)
		plt.savefig("goldstandard.pdf")
		print("Generated: goldstandard.pdf")
		plt.clf()

		os.system("e2display.py --plot "+" ".join(fscs))

	if options.resolution_vsref!=None:
		plt.title("Map vs Ref FSC")
		plt.xlabel(r"Spatial Frequency (1/$\AA$)")
		plt.ylabel("FSC")

		refines=[i for i in os.listdir(".") if "refine_" in i]
		maps=[]
		for r in refines:
			try: itr=max([i for i in os.listdir(r) if "threed_" in i and "even" not in i and "odd" not in i])
			except: continue
			maps.append("{}/{}".format(r,itr))

		maps.sort()
		
		fscs=[]
		ref=EMData(options.resolution_vsref,0,True)
		for m in maps:
			print(m)
			mi=EMData(m,0,True)
			
			# insure volumes have same sampling and box-size
			if fabs(old_div(ref["apix_x"],mi["apix_x"])-1.0)>.001 or ref["nz"]!=mi["nz"] :
				if options.verbose:
					print("{} and {} do not have the same sampling/box size. Adjusting".format(options.resolution_vsref,m))
				sca=old_div(mi["apix_x"],ref["apix_x"])
				if sca>1 : cmd="e2proc3d.py {} cmp_map.hdf --fouriershrink {} --clip {},{},{} --align translational --alignref {}".format(options.resolution_vsref,sca,mi["nx"],mi["ny"],mi["nz"],m)
				else: cmd="e2proc3d.py {} cmp_map.hdf --clip {},{},{} --scale {}  --align translational --alignref {}".format(options.resolution_vsref,mi["nx"],mi["ny"],mi["nz"],old_div(1.0,sca),m)
				launch_childprocess(cmd)
				if options.verbose>1 : print(cmd)
				refname="cmp_map.hdf"
			else: refname=options.resolution_vsref
			
			# FSC
			outname=m.replace("threed_","fsc_vsref").replace(".hdf",".txt")
			cmd="e2proc3d.py {} {} --calcfsc {}".format(refname,outname,m)
			launch_childprocess(cmd)
			if options.verbose>1 : print(cmd)
			fscs.append(outname)
			
			
		maxx=0.01

		# iterate over fsc curves
		for num,f in enumerate(fscs):
			# read the fsc curve
			d=np.loadtxt(f).transpose()

			# plot the curve
			try: plt.plot(d[0],d[1],label=f[:9],color=pltcolors[(num)%12])
			except: pass
			maxx=max(maxx,max(d[0]))

		if max(d[0])<max(xticlocs) :
			xticlocs=[i for i in xticlocs if i<=max(d[0])]
			xticlbl=xticlbl[:len(xticlocs)]
		plt.axhline(0.0,color="k")
		plt.axhline(0.143,color="#306030",linestyle=":")
		plt.axis((0,maxx,-.02,1.02))
		plt.legend(loc="upper right",fontsize="x-small")
		plt.xticks(xticlocs,xticlbl)
		plt.yticks(yticlocs,yticlbl)
		plt.savefig("vsref.pdf")
		print("Generated: vsref.pdf")
		plt.clf()
		os.system("e2display.py --plot "+" ".join(fscs))

	if options.timingbypath:
		dl=[i for i in os.listdir(".") if "spt_" in i]		# list of all refine_ directories
		dl.sort()

		for d in dl:
			try:
				jsparm=js_open_dict("{}/0_spt_params.json".format(d))
				try: cores=int(jsparm["parallel"].split(":")[1])
				except: cores=int(jsparm["threads"])
				lastiter=max([int(i.split("_")[1][:2]) for i in os.listdir(d) if "threed_" in i and len(i)==13])
				starttime=os.stat(f"{d}/model_input.hdf").st_mtime
				endtime=os.stat(f"{d}/threed_{lastiter:02d}.hdf").st_mtime
				input_ptcl=jsparm["input_ptcls"]
				if input_ptcl[-3:]=="hdf" or input_ptcl[-3:]=="lst" : nptcl=EMUtil.get_image_count(input_ptcl)
				elif input_ptcl[-4:]=="json": nptcl=len(js_open_dict(input_ptcl))
				else: nptcl=-1
				#print lastmap
				box=EMData(f"{d}/model_input.hdf",0,True)["nx"]
				targetres=jsparm["restarget"]
				niter=jsparm["niter"]
				sym=jsparm["sym"]
				tophat=jsparm.get("tophat","None")
#				maxshift=jsparm.get("maxshift","None")
#				maxang=jsparm.get("maxang","None")
				print(f"{d}\t{nptcl}\t{niter} iter\t{cores} cores\t{int((endtime-starttime)//3600)}:{int(((endtime-starttime)%3600)//60):02d}\t{cores*(endtime-starttime)/(3600*lastiter):0.1f} CPU-h/it\ttargres {targetres}\t{tophat}")
							  
			except: 
				if options.verbose: traceback.print_exc()
				print("No timing for ",d)

		print("\nWarning: scaling with number of CPUs can be very nonlinear, particularly with small jobs. The larger the number of particles the larger the number of cores which will produce near-linear speedup.")


	if options.timing:
		#dl=[i for i in os.listdir(".") if "refine_" in i]		# list of all refine_ directories
		#dl.sort()

		hist=[]
		fin=open(".eman2log.txt","r")
		while 1:
			line=fin.readline()
			if len(line)==0 : break

			spl=line.split("\t")
			try: com=spl[4].split()[0].split("/")[-1]
			except : continue

			if com in ("e2refine.py","e2refine_easy.py","e2refinemulti.py","e2project3d.py","e2simmx.py","e2simmx2stage.py","e2classesbyref.py","e2classify.py","e2classaverage.py","e2make3d.py","e2make3dpar.py","e2refine_postprocess.py") : hist.append((com,spl))

		n=0
		while n<len(hist):
			com=hist[n][1][4]
			ttime=timestamp_diff(hist[n][1][0],hist[n][1][1])
			if hist[n][0] in ("e2refine.py","e2refine_easy.py"):
				pl=com.find("--path=")
				parl=com.find("--parallel=")
				print("%s\t%1.2f hours\te2refine %s"%(difftime(ttime),old_div(ttime,3600.0),com[pl+7:].split()[0]), end=' ')
				if parl>0: print(com[parl+11:].split()[0])
				else: print(" ")

			else:
				print("\t%s\t%1.2f hours\t%s"%(difftime(ttime),old_div(ttime,3600.0),hist[n][0]))

			n+=1
			
	E2end(logid)
	sys.exit(0)


if __name__ == "__main__":
	main()
