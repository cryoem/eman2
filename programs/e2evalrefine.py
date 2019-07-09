#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

#
# Author: Steven Ludtke, 4/2/2010
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

# slow, but for printing results should be fine
def safediv(a,b):
	if b==0: return 0
	try: return old_div(a,b)
	except: return 0

#@profile
def pqual(n,ptclincls,jsd,includeproj,verbose):
	"""This computes particle quality for all particles in one class average over both iterations"""
	# The first projection is unmasked, used for scaling
	global classmx,nptcl,cmxtx,cmxty,cmxalpha,cmxmirror,eulers,threed,ptclmask,rings,pf,cptcl
	proj=[t.project("standard",eulers[n]) for t in threed]
	projmask=ptclmask.project("standard",eulers[n])		# projection of the 3-D mask for the reference volume to apply to particles

	alt=eulers[n].get_rotation("eman")["alt"]
	az=eulers[n].get_rotation("eman")["az"]

	result={}
	for it,eo,j in ptclincls:
	#for it in xrange(2):			# note that this is 0,1 not actual iteration
		#for eo in range(2):
			#for j in xrange(nptcl[eo]):
				#if classmx[eo+2*it][0,j]!=n :
##						if options.debug: print "XXX {}\t{}\t{}\t{}".format(i,("even","odd")[eo],j,classmx[eo][0,j])
					#continue		# only proceed if the particle is in this class
				if verbose >= 6: print("{}\t{}\t{}".format(cptcl[eo],("even","odd")[eo],j,it))

				truenum=j*2+eo 	# This is the particle number within the full file

				# the particle itself
				try: ptcl=EMData(cptcl[eo],j)
				except:
					print("Unable to read particle: {} ({})".format(cptcl[eo],j))
					sys.exit(1)
				try: defocus=ptcl["ctf"].defocus
				except: defocus=-1.0

				sums=[0,0]
			
				# Find the transform for this particle (2d) and apply it to the unmasked/masked projections
				ptclxf=Transform({"type":"2d","alpha":cmxalpha[eo+2*it][0,j],"mirror":int(cmxmirror[eo+2*it][0,j]),"tx":cmxtx[eo+2*it][0,j],"ty":cmxty[eo+2*it][0,j]}).inverse()

				projc=proj[it].process("xform",{"transform":ptclxf})	# we transform the projection, not the particle (as in the original classification)

				# This is for visualization with e2display later on
				if includeproj and it==1: projc.write_image(pf,truenum)

				projmaskc=projmask.process("xform",{"transform":ptclxf})
				ptcl.mult(projmaskc)

				# Particle vs projection FSC
				fsc = ptcl.calc_fourier_shell_correlation(projc)

				third = old_div(len(fsc),3)
				fsc=array(fsc[third:third*2])
#					snr=fsc/(1.0-fsc)
				result[(truenum,it)]=[old_div(sum(fsc[rings[k]:rings[k+1]]),(rings[k+1]-rings[k])) for k in range(4)]+[alt,az,n,defocus,ptcl["data_source"],ptcl["data_n"]]		# sum the fsc into 5 range values
#					sums=[sum(snr[rings[k]:rings[k+1]])/(rings[k+1]-rings[k]) for k in xrange(4)]		# sum the fsc into 5 range values

	jsd.put(result)

def main():
	global classmx,nptcl,cmxtx,cmxty,cmxalpha,cmxmirror,eulers,threed,ptclmask,rings,pf,cptcl
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] [refine_xx]
	This program performs various assessments of e2refine_easy (or similar) runs, and operates in one of several possible modes.
	A refine_xx folder name must always be provided, except for resolution_all or resolution_vsref, which operate on all 
	refine_xx folders:

	--timing
		will report how long each refinement took to complete as well as individual tasks within the refinement
		
	--timingbypath
		will report total timing information for each refine_xx folder, along with useful refinement parameters
		
	--resolution
		Computes per-iteration FSCs for a single refine_xx folder
		
	--resolution_all
		Computes FSC curves for the final iteration of every refine_xx folder
		
	--resolution_vsref
		Computes FSC curve for the final iteration of each refine_xx folder vs a provided reference volume. 
		
	--evalptclqual
		This provides the ability to assess the agreement between individual particles used for a refinement and the 
		final refinement itself in several critical ways. Using this option will automatically generate two new sets/ 
		containing "good" and "bad" particles from the population used to run the specified refinement. The assessment
		is made using a set of heuristics. A multicolumn text file is also produced containing detailed results of
		the per-particle comparison (ptclfsc_XX.txt). The first line of this text file contains details about the
		meaning of each column:
		0 - 100-30 it1 (integrated FSC from 100-30 A resolution for second to last iteration)
		1 - 30-18 it1
		2 - 18-10 it1
		3 - 10-4 it1
		4 - 100-30 it2 (integrated FSC from 100-30 A resolution for last iteration)
		5 - 30-18 it2
		6 - 18-10 it2
		7 - 10-4 it2
		8 - it12rmsd (RMSD between it1 and it2 integrated FSCs)
		9-11 - alt1, az1, cls1  (orientation parameters for second to last iteration)
		12-14 - alt2, az2, cls2 ( " for last iteration)
		15 - defocus
		
		e2display.py --plot ptclfsc*txt		(then adjust columns the plot displays, eg- 4 vs 5, 4 vs 6, 0 vs 4, 1 vs 5)

	--evalclassqual
		similar analysis to evalptclqual but for class-averages vs projections
		
	--anisotropy
		Assesses the amount and direction of any magnification anisotropy present in a raw data set by considering
		particles in a range of different orientations. Works best with large particles. Specify a class-average number
		for a class containing many particles (within a particular refine_xx folder). It is a good idea to compare results
		among classes in different highly occupied orientations.

	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="refine_xx",help="Name of a completed refine_xx folder.", default="", guitype='filebox', browser="EMRefine2dTable(withmodal=True,multiselect=False)",  filecheck=False, row=0, col=0,rowspan=1, colspan=2, mode='evalptcl')
	parser.add_argument("--timing", default=False, action="store_true", help="Report on the time required for each step of each refinement run")
	parser.add_argument("--timingbypath", default=False, action="store_true", help="Report on the CPU time required in each refine_xx folder")
	parser.add_argument("--resolution", default=False, action="store_true", help="generates a resolution and convergence plot for a single refinement run.")
	parser.add_argument("--resolution_all", default=False, action="store_true", help="generates resolution plot with the last iteration of all refine_xx directories")
	parser.add_argument("--resolution_vsref", type=str, default=None, help="Computes the FSC between the last iteration of each refine_xx directory and a specified reference map. Map must be aligned, but will be rescaled if necessary.")
	parser.add_argument("--evalptclqual", default=False, action="store_true", help="Evaluates the particle-map agreement using the refine_xx folder name. This may be used to identify bad particles.",guitype='boolbox', row=8, col=1, rowspan=1, colspan=1, mode='evalptcl[True]')
	parser.add_argument("--evalclassqual", default=False, action="store_true", help="Evaluates the class-average-projection agreement using the refine_xx folder name.",guitype='boolbox', row=8, col=2, rowspan=1, colspan=1, mode='evalptcl[False]')
	parser.add_argument("--anisotropy", type=int, default=-1, help="Specify a class-number (more particles better). Will use that class to evaluate magnification anisotropy in the data. ")
	parser.add_argument("--evalclassdetail", default=False, action="store_true", help="If specified with evalclassqual, will generate individual FRC curves for each class average in the even subset")
	parser.add_argument("--includeprojs", default=False, action="store_true", help="If specified with --evalptclqual, projections will be written to disk for easy comparison.",guitype='boolbox', row=8, col=0, rowspan=1, colspan=1, mode='evalptcl[True]')
	parser.add_argument("--iter", type=int, default=None, help="If a refine_XX folder is being used, this selects a particular refinement iteration. Otherwise the last complete iteration is used.")
	parser.add_argument("--mask",type=str,help="Mask to be used to focus --evalptclqual and other options. May be useful for separating heterogeneous data.", default=None)
	parser.add_argument("--sym",type=str,help="Symmetry to be used in searching adjacent unit cells, default from refine_xx parms", default=None)
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful",guitype='intbox', row=9, col=0, rowspan=1, colspan=1, mode='evalptcl[4]')
	#parser.add_argument("--parmcmp",  default=False, action="store_true",help="Compare parameters used in different refinement rounds")
	#parser.add_argument("--parmpair",default=None,type=str,help="Specify iter,iter to compare the parameters used between 2 itertions.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	#options associated with e2refine.py
	#parser.add_argument("--iter", dest = "iter", type = int, default=0, help = "The total number of refinement iterations to perform")
	#parser.add_argument("--check", "-c", dest="check", default=False, action="store_true",help="Checks the contents of the current directory to verify that e2refine.py command will work - checks for the existence of the necessary starting files and checks their dimensions. Performs no work ")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")
	#parser.add_argument("--input", dest="input", default=None,type=str, help="The name of the image containing the particle data")

	(options, args) = parser.parse_args()

	xticlocs=[i for i in (.01,.05,.0833,.125,.1667,.2,.25,.3333,.4,.5)]
	xticlbl=["1/100","1/20","1/12","1/8","1/6","1/5","1/4","1/3","1/2.5","1/2"]
	yticlocs=(0.0,.125,.143,.25,.375,.5,.625,.75,.875,1.0)
	yticlbl=("0"," ","0.143","0.25"," ","0.5"," ","0.75"," ","1.0")
	yticlocs2=(0.0,.125,.25,.375,.5,.625,.75,.875,1.0)
	yticlbl2=("0"," ","0.25"," ","0.5"," ","0.75"," ","1.0")

	if len(args)>0 and args[0][:7]=="refine_":
		jsparm=js_open_dict(args[0]+"/0_refine_parms.json")

		if options.iter==None:
			try:
				options.iter=int(jsparm["last_map"].split("_")[-1][:2])
				if options.sym==None :
					options.sym=str(jsparm["sym"])
			except:
				print("Could not find a completed iteration in ",args[0])
				sys.exit(1)
		
		if options.evalptclqual:
			if iter==1 :
				print("evalptclqual requires at least 2 completed iterations (3 or 4 preferred), and will use the specified --iter and the iteration preceeding it. This is not possible if --iter=1.")
				sys.exit(1)
		
		print("Using --iter=",options.iter)

	logid=E2init(sys.argv,options.ppid)

	if options.anisotropy>=0 :
		print("Anisotropy evaluation mode")

		try:
			pathmx="{}/classmx_{:02d}_even.hdf".format(args[0],options.iter)
			classmx=[EMData(pathmx,0)]
			nptcl=[classmx[0]["ny"]]
			cmxtx=[EMData(pathmx,2)]
			cmxty=[EMData(pathmx,3)]
			cmxalpha=[EMData(pathmx,4)]
			cmxmirror=[EMData(pathmx,5)]

			pathmx="{}/classmx_{:02d}_odd.hdf".format(args[0],options.iter)
			classmx.append(EMData(pathmx,0))
			nptcl.append(classmx[1]["ny"])
			cmxtx.append(EMData(pathmx,2))
			cmxty.append(EMData(pathmx,3))
			cmxalpha.append(EMData(pathmx,4))
			cmxmirror.append(EMData(pathmx,5))
		except:
			traceback.print_exc()
			print("====\nError reading classification matrix. Must be full classification matrix with alignments")
			sys.exit(1)

		if options.verbose: print("{} even and {} odd particles in classmx".format(nptcl[0],nptcl[1]))

		# path to the even/odd particles used for the refinement
		cptcl=jsparm["input"]
		cptcl=[str(i) for i in cptcl]

		# this reads all of the EMData headers from the projections, should be same for even and odd
		pathprj="{}/projections_{:02d}_even.hdf".format(args[0],options.iter)
		nref=EMUtil.get_image_count(pathprj)
		eulers=[EMData(pathprj,i,1)["xform.projection"] for i in range(nref)]

		# The 3D reference volume we are using for subtraction
		threed=EMData("{}/threed_{:02d}.hdf".format(args[0],options.iter),0)

		# The mask applied to the reference volume, used for 2-D masking of particles for better power spectrum matching
		if options.mask: ptclmask=EMData(options.mask,0)
		else: ptclmask=EMData(args[0]+"/mask.hdf",0)
		nx=ptclmask["nx"]
		apix=threed["apix_x"]

		# We expand the mask a bit, since we want to consider problems with "touching" particles
		ptclmask.process_inplace("threshold.binary",{"value":0.2})
		ptclmask.process_inplace("mask.addshells",{"nshells":nx//15})
		ptclmask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})

		ring=(int(old_div(2*nx*apix,100.0)),int(old_div(2*nx*apix,10)))
#		fout=open("ptclsnr.txt".format(i),"w")
		fout=open("aniso_{:02d}.txt".format(options.anisotropy),"w")
		# generate a projection for each particle so we can compare
		for i in [options.anisotropy]:							# this is left as a loop in case we decide to do multiple classes later on
			if options.verbose>1 : print("--- Class %d"%i)

			# The first projection is unmasked, used for scaling
			proj=threed.project("standard",{"transform":eulers[i]})
			projmask=ptclmask.project("standard",eulers[i])		# projection of the 3-D mask for the reference volume to apply to particles

			alt=eulers[i].get_rotation("eman")["alt"]
			az=eulers[i].get_rotation("eman")["az"]
			best=(0,0,1.02)

			for angle in range(0,180,5):
				rt=Transform({"type":"2d","alpha":angle})
				xf=rt*Transform([1.02,0,0,0,0,old_div(1,1.02),0,0,0,0,1,0])*rt.inverse()
				esum=0

				for eo in range(2):
					for j in range(nptcl[eo]):
						if classmx[eo][0,j]!=i :
	#						if options.debug: print "XXX {}\t{}\t{}\t{}".format(i,("even","odd")[eo],j,classmx[eo][0,j])
							continue		# only proceed if the particle is in this class
						if options.verbose: print("{}\t{}\t{}".format(i,("even","odd")[eo],j))

						# the particle itself
						try: ptcl=EMData(cptcl[eo],j)
						except:
							print("Unable to read particle: {} ({})".format(cptcl[eo],j))
							sys.exit(1)
						#try: defocus=ptcl["ctf"].defocus
						#except: defocus=-1.0

						# Find the transform for this particle (2d) and apply it to the unmasked/masked projections
						ptcl.transform(xf)		# anisotropy directly on the particle
						ptclxf=Transform({"type":"2d","alpha":cmxalpha[eo][0,j],"mirror":int(cmxmirror[eo][0,j]),"tx":cmxtx[eo][0,j],"ty":cmxty[eo][0,j]}).inverse()
						projc=proj.process("xform",{"transform":ptclxf})	# we transform the projection, not the particle (as in the original classification)

						projmaskc=projmask.process("xform",{"transform":ptclxf})
						ptcl.mult(projmaskc)


						# Particle vs projection FSC
						fsc = ptcl.calc_fourier_shell_correlation(projc)
						third = old_div(len(fsc),3)
						fsc=array(fsc[third:third*2])
						try: esum+= sum(fsc[ring[0]:ring[1]])
						except:
							print("error")
							print(ring,fsc)
							sys.exit(1)

						best=max(best,(esum,angle,1.02))
	#					snr=fsc/(1.0-fsc)
	#					sums=[sum(fsc[rings[k]:rings[k+1]])/(rings[k+1]-rings[k]) for k in xrange(4)]		# sum the fsc into 5 range values
	#					sums=[sum(snr[rings[k]:rings[k+1]])/(rings[k+1]-rings[k]) for k in xrange(4)]		# sum the fsc into 5 range values
	#					fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t# {};{}\n".format(sums[0],sums[1],sums[2],sums[3],alt,az,i,defocus,j,cptcl[eo]))
				fout.write("{}\t{}\t{}\n".format(angle,1.02,old_div(esum,(nptcl[0]+nptcl[1]))))

			if options.verbose>1 : print("--- Class %d"%i)

			angle=best[1]
			print(best)

			for aniso in range(0,30):
				ai=old_div(aniso,1000.0)+1.0
				rt=Transform({"type":"2d","alpha":angle})
				xf=rt*Transform([ai,0,0,0,0,old_div(1,ai),0,0,0,0,1,0])*rt.inverse()
				esum=0

				for eo in range(2):
					for j in range(nptcl[eo]):
						if classmx[eo][0,j]!=i :
	#						if options.debug: print "XXX {}\t{}\t{}\t{}".format(i,("even","odd")[eo],j,classmx[eo][0,j])
							continue		# only proceed if the particle is in this class
						if options.verbose: print("{}\t{}\t{}".format(i,("even","odd")[eo],j))

						# the particle itself
						try: ptcl=EMData(cptcl[eo],j)
						except:
							print("Unable to read particle: {} ({})".format(cptcl[eo],j))
							sys.exit(1)
						#try: defocus=ptcl["ctf"].defocus
						#except: defocus=-1.0

						# Find the transform for this particle (2d) and apply it to the unmasked/masked projections
						ptcl.transform(xf)		# anisotropy directly on the particle
						ptclxf=Transform({"type":"2d","alpha":cmxalpha[eo][0,j],"mirror":int(cmxmirror[eo][0,j]),"tx":cmxtx[eo][0,j],"ty":cmxty[eo][0,j]}).inverse()
						projc=proj.process("xform",{"transform":ptclxf})	# we transform the projection, not the particle (as in the original classification)

						projmaskc=projmask.process("xform",{"transform":ptclxf})
						ptcl.mult(projmaskc)


						# Particle vs projection FSC
						fsc = ptcl.calc_fourier_shell_correlation(projc)
						third = old_div(len(fsc),3)
						fsc=array(fsc[third:third*2])
						esum+= sum(fsc[ring[0]:ring[1]])

						best=max(best,(esum,angle,ai))

				fout.write("{}\t{}\t{}\n".format(angle,ai,old_div(esum,(nptcl[0]+nptcl[1]))))

			print(best)

	if options.evalptclqual:
#		from multiprocessing import Pool
		import threading,queue
		print("Particle quality evaluation mode")

		jsparm=js_open_dict(args[0]+"/0_refine_parms.json")

		if options.iter==None:
			try:
				options.iter=int(jsparm["last_map"].split("_")[-1][:2])
				options.sym=jsparm["sym"]
			except:
				print("Could not find a completed iteration in ",args[0])
				sys.exit(1)
		
		# This is not great programming process, but greatly simplifies threading, and reduces potential memory usage

		try:
			pathmx=["{}/classmx_{:02d}_even.hdf".format(args[0],options.iter-1),"{}/classmx_{:02d}_odd.hdf".format(args[0],options.iter-1),"{}/classmx_{:02d}_even.hdf".format(args[0],options.iter),"{}/classmx_{:02d}_odd.hdf".format(args[0],options.iter)]
			classmx=[EMData(f,0) for f in pathmx]
			nptcl=[classmx[i]["ny"] for i in range(len(pathmx))]
			cmxtx=[EMData(f,2) for f in pathmx]
			cmxty=[EMData(f,3) for f in pathmx]
			cmxalpha=[EMData(f,4) for f in pathmx]
			cmxmirror=[EMData(f,5) for f in pathmx]

		except:
			traceback.print_exc()
			print("====\nError reading classification matrix. Must be full classification matrix with alignments")
			sys.exit(1)

		if options.verbose: print("{} even and {} odd particles in classmx".format(nptcl[0],nptcl[1]))


		# path to the even/odd particles used for the refinement
		cptcl=jsparm["input"]
		cptcl=[str(i) for i in cptcl]

		# this reads all of the EMData headers from the projections, should be same for even and odd, and assume same across iterations
		pathprj="{}/projections_{:02d}_even.hdf".format(args[0],options.iter)
		nref=EMUtil.get_image_count(pathprj)
		eulers=[EMData(pathprj,i,1)["xform.projection"] for i in range(nref)]

		# The 3D reference volume we are using for subtraction
		threed=[EMData("{}/threed_{:02d}.hdf".format(args[0],options.iter-1),0),EMData("{}/threed_{:02d}.hdf".format(args[0],options.iter),0)]

		# The mask applied to the reference volume, used for 2-D masking of particles for better power spectrum matching, we'll assume the mask doesn't change much
		if options.mask: ptclmask=EMData(options.mask,0)
		else: ptclmask=EMData(args[0]+"/mask.hdf",0)
		nx=ptclmask["nx"]
		apix=threed[0]["apix_x"]

		rings=[int(old_div(2*nx*apix,res)) for res in (100,30,15,8,4)]
		print(("Frequency Bands: {lowest},{low},{mid},{high},{highest}".format(lowest=rings[0],low=rings[1],mid=rings[2],high=rings[3],highest=rings[4])))

		# We expand the mask a bit, since we want to consider problems with "touching" particles
		ptclmask.process_inplace("threshold.binary",{"value":0.2})
		ptclmask.process_inplace("mask.addshells",{"nshells":nx//15})
		ptclmask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})

#		try: os.mkdir("ptclfsc")
#		except: pass

#		fout=open("ptclsnr.txt".format(i),"w")
		try: rfnnum="_".join(args[0].split("_")[1:])
		except: rfnnum="X"
		ptclfsc = "ptclfsc_{}.txt".format(rfnnum)
		# generate a projection for each particle so we can compare

		pf = "ptclfsc_{}_projections.hdf".format(rfnnum)

		if options.sym==None or len(options.sym)==0 : options.sym="c1"
		try:
			sym=parsesym(options.sym).get_syms()
			if options.verbose: print("Using symmetry: ",options.sym)
		except:
			print("Symmetry parsing error! ",options.sym)
			sys.exit(1)
			
		#tfs = []

		tlast=time()
		# Put particles in class lists
		classptcls={}
		for it in range(2):			# note that this is 0,1 not actual iteration
			for eo in range(2):
				for j in range(nptcl[eo]):
					cls=int(classmx[eo+2*it][0,j])
					try: classptcls[cls].append((it,eo,j))
					except: classptcls[cls]=[(it,eo,j)]

		# Create Thread objects
		jsd=queue.Queue(0)
		thrds=[threading.Thread(target=pqual,args=(i,classptcls[i],jsd,options.includeprojs,options.verbose)) for i in list(classptcls.keys())]
		result={}
		thrtolaunch=0

		while thrtolaunch<len(thrds) or threading.active_count()>1:
			# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
			# note that it's ok that we wait here forever, since there can't be new results if an existing
			# thread hasn't finished.
			if thrtolaunch<len(thrds) :
				while (threading.active_count()==options.threads+1 ) : sleep(.01)
				if options.verbose : 
					print(" Starting thread {}/{}      \r".format(thrtolaunch,len(thrds)), end=' ')
					sys.stdout.flush()
				thrds[thrtolaunch].start()
				thrtolaunch+=1
			else:sleep(.25)
		
			while not jsd.empty():
				rd=jsd.get()
				result.update(rd)
		if options.verbose: print("Threads complete             ")

		for t in thrds:
			t.join()
	 
		fout=open(ptclfsc,"w")
		fout.write("# 100-30 it1; 30-15 it1; 15-8 it1; 8-4 it1; 100-30 it2; 30-15 it2; 15-8 it2; 8-4 it2; it12rmsd; (30-8)/(100-30) alt1; az1; cls1; alt2; az2; cls2; defocus; angdiff\n")

		# loop over all particles and print results
		rmsds=[]
		for j in range(nptcl[0]+nptcl[1]):
			try:
				r=result[(j,0)][:8]+result[(j,1)]		# we strip out the filename and number from the first result
			except:
				print("Missing results ptcl:",j, end=' ')
				try:
					print(result[(j,0)], end=' ')
					print(result[(j,1)])
				except: print(" ")
				continue
			jj=old_div(j,2)
			eo=j%2
			# This is the rotation to go from the first orientation to the second, ignoring phi, taking symmetry into account
			adiff=angle_ab_sym(sym,r[4],r[5],r[12],r[13])
			rmsd=sqrt((r[0]-r[8])**2+(r[1]-r[9])**2+(r[2]-r[10])**2+(r[3]-r[11])**2)
			rmsds.append(rmsd)
			if options.includeprojs:
				fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t# {};{};{};{}\n".format(r[0],r[1],r[2],r[3],r[8],r[9],r[10],r[11],rmsd,old_div((r[9]+r[10]),r[8]),r[4],r[5],r[6],r[12],r[13],r[14],r[15],adiff,jj,cptcl[eo],j,pf))
			else:
				fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t# {};{}\n".format(r[0],r[1],r[2],r[3],r[8],r[9],r[10],r[11],rmsd,old_div((r[9]+r[10]),r[8]),r[4],r[5],r[6],r[12],r[13],r[14],r[15],adiff,jj,cptcl[eo]))

		fout.close()

		####  Done writing results to text file, now we generate new sets
		print("Generating new sets")
		an=Analyzers.get("kmeans")
		an.set_params({"ncls":4,"seedmode":1,"minchange":len(rmsds)//100,"verbose":0,"slowseed":0,"mininclass":5})
		quals=[]
		for j in range(nptcl[0]+nptcl[1]):
			d=EMData(3,1,1)
			# We use the first 3 resolution bands, taking the max value from the 2 iterations, and upweight the second by 2x, since higher res info
			# is important, and very low res is impacted by defocus
			d[0]=max(result[(j,0)][0],result[(j,1)][0])
			d[1]=max(result[(j,0)][1],result[(j,1)][1])*2.0
			d[2]=max(result[(j,0)][2],result[(j,1)][2])
			quals.append(d)
		an.insert_images_list(quals)

		centers=an.analyze()
		print("Centers: {}({:1.3f},{:1.3f},{:1.3f}), {}({:1.3f},{:1.3f},{:1.3f}), {}({:1.3f},{:1.3f},{:1.3f}, {}({:1.3f},{:1.3f},{:1.3f})".format(
			centers[0]["ptcl_repr"],centers[0][0],centers[0][1],centers[0][2],
			centers[1]["ptcl_repr"],centers[1][0],centers[1][1],centers[1][2],
			centers[2]["ptcl_repr"],centers[2][0],centers[2][1],centers[2][2],
			centers[3]["ptcl_repr"],centers[3][0],centers[3][1],centers[3][2] ))
		
		badcls=min([(centers[i]["mean"],i) for i in (0,1,2,3)])[1]	# this confusing expression finds the number of the class with the smallest summed vector
		print("Class {} is the bad class".format(badcls))
		
		rmsds=array(rmsds)
		rmsdthresh=rmsds.std()*2.5
		print("Within consistency thr {:0.4f}: {}/{}".format(rmsdthresh,len(rmsds[rmsds<rmsdthresh]),len(rmsds)))

		# write actual classes
		nameg = "sets/pf{}_good_{}".format(rfnnum,cptcl[0].replace("_even","").replace("sets/",""))
		namegb=nameg.split("__")[0]+"__ctf_flip_invar.lst"
		nameb = "sets/pf{}_bad_{}".format(rfnnum,cptcl[0].replace("_even","").replace("sets/",""))
		namebb=nameb.split("__")[0]+"__ctf_flip_invar.lst"
		try: os.unlink(nameg)
		except: pass
		try: os.unlink(nameb)
		except: pass
		outb=LSXFile(nameb)
		outg=LSXFile(nameg)
		ngood=0
		for i,q in enumerate(quals):
			r=result[(i,0)]
			if q["class_id"]==badcls or rmsds[i]>rmsdthresh: 
				outb.write(-1,r[-1],r[-2],"{:6.4f},{:6.4f},{:6.4f},{:6.4f}".format(quals[i][0],quals[i][1],quals[i][2],rmsds[i]))
			else: 
				outg.write(-1,r[-1],r[-2],"{:6.4f},{:6.4f},{:6.4f},{:6.4f}".format(quals[i][0],quals[i][1],quals[i][2],rmsds[i]))
				ngood+=1
		outb.close()
		outg.close()
		
		open(namegb,"w").write(open(nameg,"r").read())		# copy the file
		launch_childprocess("e2proclst.py {} --retype ctf_flip_invar".format(namegb))
		
		open(namebb,"w").write(open(nameb,"r").read())		# copy the file
		launch_childprocess("e2proclst.py {} --retype ctf_flip_invar".format(namebb))
		
		print("{}/{} kept as good".format(ngood,len(quals)))

		print("Evaluation complete.\nParticles best resembling results from {} at low/intermediate resolution have been saved in {} and can be used in further refinements.\nNote that this method will identify the worst particles as bad, regardless of whether they actually are (bad), and that it may be wise to do your own classification on these results instead, as described in the tutorial.".format(args[0],nameg))



	if options.evalclassqual:
		print("Class quality evaluation mode")


		classes=["{}/classes_{:02d}_{}.hdf".format(args[0],options.iter,i) for i in ("even","odd")]
		projections=["{}/projections_{:02d}_{}.hdf".format(args[0],options.iter,i) for i in ("even","odd")]
		n=EMUtil.get_image_count(projections[0])
		timg=EMData(projections[0],0)

		# path to the even/odd particles used for the refinement
		cptcl=jsparm["input"]
		cptcl=[str(i) for i in cptcl]

		# The mask applied to the reference volume, used for 2-D masking of particles for better power spectrum matching
		nx=timg["nx"]
		apix=timg["apix_x"]

		rings=[int(old_div(2*nx*apix,res)) for res in (100,30,12,7,3)]
#		rings=[int(2*nx*apix/res) for res in (100,60,30,15,10,6,4)]
		print(("Frequency Bands: {}".format(rings)))

		classfsc="classfsc_{}_{}.txt".format(args[0][-2:],options.iter)
		fout=open(classfsc,"w")
		fout.write("# 100-30 A; 30-12 A; 12-7 A; 7-3 A; (30-12)/(100-30); Nptcl; az; alt; phi; SSNR/Nptcl 100-30 A; SSNR/Nptcl 30-15 A; SSNR/Nptcl 15-8 A\n")
		for eo in range(2):
			for i in range(n):
				cl=EMData(classes[eo],i)
				pr=EMData(projections[eo],i)
				alt=pr["xform.projection"].get_rotation("eman")["alt"]
				az=pr["xform.projection"].get_rotation("eman")["az"]
				try:
					phi=cl["xform.projection"].get_rotation("eman")["phi"]
				except: phi=0
				
				fsc = cl.calc_fourier_shell_correlation(pr)

				third = old_div(len(fsc),3)
				fsc=array(fsc[third:third*2])
				#sums=[sum(fsc[rings[k]:rings[k+1]])/(rings[k+1]-rings[k]) for k in xrange(4)]		# sum the fsc into 5 range values
				#fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t# {};{};{};{}\n".format(sums[0],sums[1],sums[2],sums[3],cl["ptcl_repr"],alt,az,phi,sums[0]/(1.0001-sums[0])/(cl["ptcl_repr"]+0.01),sums[1]/(1.0001-sums[1])/(cl["ptcl_repr"]+0.01),sums[2]/(1.0001-sums[2])/(cl["ptcl_repr"]+0.01),i,classes[eo],i,projections[eo]))
				sums=[old_div(sum(fsc[rings[k]:rings[k+1]]),(rings[k+1]-rings[k])) for k in range(4)]		# sum the fsc into 5 range values
				fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t# {};{};{};{}\n".format(sums[0],sums[1],sums[2],sums[3],safediv(sums[1],sums[0]),cl["ptcl_repr"],alt,az,phi,safediv(safediv(sums[0],(1.0-sums[0])),cl["ptcl_repr"]),safediv(safediv(sums[1],(1.0-sums[1])),cl["ptcl_repr"]),safediv(safediv(sums[2],(1.0-sums[2])),cl["ptcl_repr"]),i,classes[eo],i,projections[eo]))

				if options.evalclassdetail and eo==0:
					out=open("cfsc{:04d}.txt".format(i),"w")
					fsc=cl.calc_fourier_shell_correlation(pr)
					third=old_div(len(fsc),3)
					ssnr=[fsc[third+1]]*5+fsc[third+1:third*2]+[fsc[third*2-1]]*4		# we extend the list by replication to make the running average more natural
#					print(len(ssnr),third)
					npnt=[fsc[third*2+1]]*5+fsc[third*2+1:third*3]+[fsc[-1]]*4	# number of points in each average
					try:
						ssnr=[old_div(sum([ssnr[k]*npnt[k] for k in range(j-4,j+5)]),sum([npnt[k] for k in range(j-4,j+5)])) for j in range(4,third+4)]			# smoothing by weighted running average
					except:
						ssnr=[0,0]
					ssnr=[old_div(v,(1.0-min(v,.999999))) for v in ssnr]							# convert FSC to pseudo SSNR
					for x,v in enumerate(ssnr): out.write("{}\t{}\n".format(x,v))
					out.close()
				

	if options.resolution:

		if not os.path.isdir(args[0]):
			print("You must provide the name of the refine_XX folder")
			sys.exit(1)

		### Convergenece plot

		plt.title("Convergence plot (not resolution)")
		plt.xlabel(r"Spatial Frequency (1/$\AA$)")
		plt.ylabel("FSC")
		cnvrg=[i for i in os.listdir(args[0]) if "converge_" in i and i[-4:]==".txt"]
		cnvrg.sort(reverse=True)
		nummx=int(cnvrg[0].split("_")[2][:2])
		maxx=0.01
		for c in cnvrg:
			num=int(c.split("_")[2][:2])
			d=np.loadtxt("{}/{}".format(args[0],c)).transpose()
			if c[9:13]=="even" : plt.plot(d[0],d[1],label=c[14:-4],color=pltcolors[(nummx-num)%12])
			else : plt.plot(d[0],d[1],color=pltcolors[(nummx-num)%12])
			maxx=max(maxx,max(d[0]))

		if max(d[0])<max(xticlocs) :
			xticlocs=[i for i in xticlocs if i<=max(d[0])]
			xticlbl=xticlbl[:len(xticlocs)]
		plt.axhline(0.0,color="k")
		plt.axis((0,maxx,-.06,1.02))
		plt.legend(loc="upper right",fontsize="x-small")
		#plt.minorticks_on()
		plt.xticks(xticlocs,xticlbl)
		plt.yticks(yticlocs2,yticlbl2)
		plt.savefig("converge_{}.pdf".format(args[0][-2:]))
		print("Generated : converge_{}.pdf".format(args[0][-2:]))
		plt.clf()

		######################
		### Resolution plot
		### broken up into multiple try/except blocks because we need some of the info, even if plotting fails
		plt.title("Gold Standard Resolution")
		plt.xlabel(r"Spatial Frequency (1/$\AA$)")
		plt.ylabel("FSC")

		fscs=[i for i in os.listdir(args[0]) if "fsc_masked" in i and i[-4:]==".txt"]
		fscs.sort(reverse=True)
		nummx=int(fscs[0].split("_")[2][:2])
		maxx=0.01

		# iterate over fsc curves
		for f in fscs:
			num=int(f.split("_")[2][:2])

			# read the fsc curve
			d=np.loadtxt("{}/{}".format(args[0],f)).transpose()

			# plot the curve
			try: plt.plot(d[0],d[1],label=f[4:],color=pltcolors[(nummx-num)%12])
			except: pass
			maxx=max(maxx,max(d[0]))

			# find the resolution from the first curve (the highest numbered one)
			if f==fscs[0]:
				# find the 0.143 crossing
				for si in range(2,len(d[0])-2):
					if d[1][si-1]>0.143 and d[1][si]<=0.143 :
						frac=old_div((0.143-d[1][si]),(d[1][si-1]-d[1][si]))		# 1.0 if 0.143 at si-1, 0.0 if .143 at si
						lastres=d[0][si]*(1.0-frac)+d[0][si-1]*frac
						try:
							plt.annotate(r"{:1.1f} $\AA$".format(old_div(1.0,lastres)),xy=(lastres,0.143),
								xytext=(old_div((lastres*4+d[0][-1]),5.0),0.2),arrowprops={"width":1,"frac":.1,"headwidth":7,"shrink":.05})
						except: pass
						break
				else : lastres=0

		plt.axhline(0.0,color="k")
		plt.axhline(0.143,color="#306030",linestyle=":")
		plt.axis((0,maxx,-.02,1.02))
		plt.legend(loc="upper right",fontsize="x-small")
		plt.xticks(xticlocs,xticlbl)
		plt.yticks(yticlocs,yticlbl)
		plt.savefig("goldstandard_{}.pdf".format(args[0][-2:]))
		print("Generated: goldstandard_{}.pdf".format(args[0][-2:]))
		plt.clf()

	if options.resolution_all:
		######################
		### Resolution plot
		### broken up into multiple try/except blocks because we need some of the info, even if plotting fails
		plt.title("Gold Standard Resolution")
		plt.xlabel(r"Spatial Frequency (1/$\AA$)")
		plt.ylabel("FSC")

		refines=[i for i in os.listdir(".") if "refine_" in i]
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
		dl=[i for i in os.listdir(".") if "refine_" in i]		# list of all refine_ directories
		dl.sort()

		for d in dl:
			try:
				jsparm=js_open_dict("{}/0_refine_parms.json".format(d))
				try: cores=int(jsparm["parallel"].split(":")[1])
				except: cores=int(jsparm["threads"])
				lastmap=str(jsparm["last_map"])
				lastiter=int(lastmap.split("/")[-1].split("_")[-1][:2])
				firstmap="{}/threed_00_even.hdf".format(d)
				starttime=os.stat(firstmap).st_mtime
				endtime=os.stat(lastmap).st_mtime
				#print lastmap
				box=EMData(lastmap,0,True)["nx"]
				targetres=jsparm["targetres"]
				speed=jsparm["speed"]
				bispec="bispec" if jsparm.getdefault("bispec",False) else " "
				if bispec==" ": bispec="invar" if jsparm.getdefault("invar",False) else " "
				nptcl=EMUtil.get_image_count(str(jsparm["input"][0]))+EMUtil.get_image_count(str(jsparm["input"][1]))

				print("{path}\t{nptcl} ptcls\t{niter} iter\t{cores} cores\t{h:02d}:{m:02d} walltime\t{cpuh:1.1f} CPU-h\t{cpuhpi:1.2f} CPU-h/it\t{bs} box\t{targ:1.1f} targetres\tspd={speed} {bispec}".format(
					path=d,niter=lastiter,cores=cores,h=int((endtime-starttime)//3600),m=int(((endtime-starttime)%3600)//60),
					cpuh=old_div(cores*(endtime-starttime),3600),cpuhpi=old_div(cores*(endtime-starttime),(3600*lastiter)),bs=box,targ=targetres,speed=speed,bispec=bispec,nptcl=nptcl))
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
