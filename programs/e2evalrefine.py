#!/usr/bin/env python

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


from EMAN2 import *
from EMAN2db import db_open_dict
from math import *
import os
import sys
import datetime
from numpy import array
import traceback

try:
	import numpy as np
	import matplotlib
	matplotlib.use("AGG")
#	matplotlib.use("PDF")
	import matplotlib.pyplot as plt
	pltcolors=["k","b","g","r","m","c","darkblue","darkgreen","darkred","darkmagenta","darkcyan","0.5"]
except:
	print "Matplotlib not available, plotting options will not be available"


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] [refine_xx]
	This program is still in its early stages. Eventually will provide a variety of tools for 
	evaluating a single particle reconstruction refinement run."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--evalptclqual", default=False, action="store_true", help="Evaluates the particle-map agreement the refine_xx folder name. This may be used to identify bad particles.")
	parser.add_argument("--anisotropy", type=int, default=-1, help="Specify a class-number (more particles better). Will use that class to evaluate magnification anisotropy in the data. ")
	parser.add_argument("--iter", type=int, default=None, help="If a refine_XX folder is being used, this selects a particular refinement iteration. Otherwise the last complete iteration is used.")
	parser.add_argument("--timing", default=False, action="store_true", help="report on how long each step of the refinement process took during the first iteration of each run")
	parser.add_argument("--resolution", default=False, action="store_true", help="generates a resolution and convergence plot for a single refinement run.")
	parser.add_argument("--resolution_all", default=False, action="store_true", help="generates resolution plot with the last iteration of all refine_xx directories")
	#parser.add_argument("--parmcmp",  default=False, action="store_true",help="Compare parameters used in different refinement rounds")
	#parser.add_argument("--parmpair",default=None,type=str,help="Specify iter,iter to compare the parameters used between 2 itertions.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	#options associated with e2refine.py
	#parser.add_argument("--iter", dest = "iter", type = int, default=0, help = "The total number of refinement iterations to perform")
	#parser.add_argument("--check", "-c", dest="check", default=False, action="store_true",help="Checks the contents of the current directory to verify that e2refine.py command will work - checks for the existence of the necessary starting files and checks their dimensions. Performs no work ")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	#parser.add_argument("--input", dest="input", default=None,type=str, help="The name of the image containing the particle data")

	(options, args) = parser.parse_args()

	xticlocs=[i for i in (.01,.05,.0833,.125,.1667,.2,.25,.3333,.4,.5)]
	xticlbl=["1/100","1/20","1/12","1/8","1/6","1/5","1/4","1/3","1/2.5","1/2"]
	yticlocs=(0.0,.125,.143,.25,.375,.5,.625,.75,.875,1.0)
	yticlbl=("0"," ","0.143","0.25"," ","0.5"," ","0.75"," ","1.0")
	yticlocs2=(0.0,.125,.25,.375,.5,.625,.75,.875,1.0)
	yticlbl2=("0"," ","0.25"," ","0.5"," ","0.75"," ","1.0")

	if args[0][:7]=="refine_":
		jsparm=js_open_dict(args[0]+"/0_refine_parms.json")
		
		if options.iter==None: 
			try: options.iter=int(jsparm["last_map"].split("_")[-1][:2])
			except:
				print "Could not find a completed iteration in ",args[0]
				sys.exit(1)

	if options.anisotropy>=0 :
		print "Anisotropy evaluation mode"

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
			print "====\nError reading classification matrix. Must be full classification matrix with alignments"
			sys.exit(1)

		if options.verbose: print "{} even and {} odd particles in classmx".format(nptcl[0],nptcl[1])

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
		ptclmask=EMData(args[0]+"/mask.hdf",0)
		nx=ptclmask["nx"]
		apix=threed["apix_x"]

		# We expand the mask a bit, since we want to consider problems with "touching" particles
		ptclmask.process_inplace("threshold.binary",{"value":0.2})
		ptclmask.process_inplace("mask.addshells",{"nshells":nx//15})
		ptclmask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})

		ring=(2*nx*apix/100.0,2*nx*apix/10)
#		fout=open("ptclsnr.txt".format(i),"w")
		fout=open("aniso_{:02d}.txt".format(options.anisotropy),"w")
		# generate a projection for each particle so we can compare
		for i in [options.anisotropy]:							# this is left as a loop in case we decide to do multiple classes later on
			if options.verbose>1 : print "--- Class %d"%i

			# The first projection is unmasked, used for scaling
			proj=threed.project("standard",{"transform":eulers[i]})
			projmask=ptclmask.project("standard",eulers[i])		# projection of the 3-D mask for the reference volume to apply to particles

			alt=eulers[i].get_rotation("eman")["alt"]
			az=eulers[i].get_rotation("eman")["az"]
			best=(0,0,1.01)

			for angle in xrange(0,180,5):
				rt=Transform({"type":"2d","alpha":angle})
				xf=rt*Transform([1.01,0,0,0,0,1/1.01,0,0,0,0,1,0])*rt.inverse()
				esum=0
				
				for eo in range(2):
					for j in xrange(nptcl[eo]):
						if classmx[eo][0,j]!=i : 
	#						if options.debug: print "XXX {}\t{}\t{}\t{}".format(i,("even","odd")[eo],j,classmx[eo][0,j])
							continue		# only proceed if the particle is in this class
						if options.verbose: print "{}\t{}\t{}".format(i,("even","odd")[eo],j)

						# the particle itself
						try: ptcl=EMData(cptcl[eo],j)
						except:
							print "Unable to read particle: {} ({})".format(cptcl[eo],j)
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
						third = len(fsc)/3
						fsc=array(fsc[third:third*2])
						esum+= sum(fsc[ring[0]:ring[1]])
						
						best=max(best,(esum,angle,1.01))
	#					snr=fsc/(1.0-fsc)
	#					sums=[sum(fsc[rings[k]:rings[k+1]])/(rings[k+1]-rings[k]) for k in xrange(4)]		# sum the fsc into 5 range values
	#					sums=[sum(snr[rings[k]:rings[k+1]])/(rings[k+1]-rings[k]) for k in xrange(4)]		# sum the fsc into 5 range values
	#					fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t# {};{}\n".format(sums[0],sums[1],sums[2],sums[3],alt,az,i,defocus,j,cptcl[eo]))
				fout.write("{}\t{}\t{}\n".format(angle,1.01,esum))

			if options.verbose>1 : print "--- Class %d"%i

			angle=best[1]
			print best

			for aniso in xrange(0,30):
				ai=aniso/1000.0+1.0
				rt=Transform({"type":"2d","alpha":angle})
				xf=rt*Transform([ai,0,0,0,0,1/ai,0,0,0,0,1,0])*rt.inverse()
				esum=0
				
				for eo in range(2):
					for j in xrange(nptcl[eo]):
						if classmx[eo][0,j]!=i : 
	#						if options.debug: print "XXX {}\t{}\t{}\t{}".format(i,("even","odd")[eo],j,classmx[eo][0,j])
							continue		# only proceed if the particle is in this class
						if options.verbose: print "{}\t{}\t{}".format(i,("even","odd")[eo],j)

						# the particle itself
						try: ptcl=EMData(cptcl[eo],j)
						except:
							print "Unable to read particle: {} ({})".format(cptcl[eo],j)
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
						third = len(fsc)/3
						fsc=array(fsc[third:third*2])
						esum+= sum(fsc[ring[0]:ring[1]])
						
						best=max(best,(esum,angle,ai))

				fout.write("{}\t{}\t{}\n".format(angle,ai,esum))

			print best
		sys.exit(0)



	if options.evalptclqual!=None:
		print "Particle quality evaluation mode"
		
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
			print "====\nError reading classification matrix. Must be full classification matrix with alignments"
			sys.exit(1)

		if options.verbose: print "{} even and {} odd particles in classmx".format(nptcl[0],nptcl[1])

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
		ptclmask=EMData(args[0]+"/mask.hdf",0)
		nx=ptclmask["nx"]
		apix=threed["apix_x"]

		rings=[int(2*nx*apix/res) for res in (100,30,15,8,4)]
		print rings

		# We expand the mask a bit, since we want to consider problems with "touching" particles
		ptclmask.process_inplace("threshold.binary",{"value":0.2})
		ptclmask.process_inplace("mask.addshells",{"nshells":nx//15})
		ptclmask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.25})

#		try: os.mkdir("ptclfsc")
#		except: pass

#		fout=open("ptclsnr.txt".format(i),"w")
		fout=open("ptclfsc_{}.txt".format(args[0][-2:]),"w")
		# generate a projection for each particle so we can compare
		for i in xrange(nref):
			if options.verbose>1 : print "--- Class %d"%i

			# The first projection is unmasked, used for scaling
			proj=threed.project("standard",{"transform":eulers[i]})
			projmask=ptclmask.project("standard",eulers[i])		# projection of the 3-D mask for the reference volume to apply to particles

			alt=eulers[i].get_rotation("eman")["alt"]
			az=eulers[i].get_rotation("eman")["az"]

#			fout=open("ptclfsc/f{:04d}.txt".format(i),"w")
			for eo in range(2):
				for j in xrange(nptcl[eo]):
					if classmx[eo][0,j]!=i : 
#						if options.debug: print "XXX {}\t{}\t{}\t{}".format(i,("even","odd")[eo],j,classmx[eo][0,j])
						continue		# only proceed if the particle is in this class
					if options.verbose: print "{}\t{}\t{}".format(i,("even","odd")[eo],j)

					# the particle itself
					try: ptcl=EMData(cptcl[eo],j)
					except:
						print "Unable to read particle: {} ({})".format(cptcl[eo],j)
						sys.exit(1)
					try: defocus=ptcl["ctf"].defocus
					except: defocus=-1.0

					# Find the transform for this particle (2d) and apply it to the unmasked/masked projections
					ptclxf=Transform({"type":"2d","alpha":cmxalpha[eo][0,j],"mirror":int(cmxmirror[eo][0,j]),"tx":cmxtx[eo][0,j],"ty":cmxty[eo][0,j]}).inverse()
					projc=proj.process("xform",{"transform":ptclxf})	# we transform the projection, not the particle (as in the original classification)

					projmaskc=projmask.process("xform",{"transform":ptclxf})
					ptcl.mult(projmaskc)

					#ptcl.write_image("tst.hdf",0)
					#projc[0].write_image("tst.hdf",1)

					# Particle vs projection FSC
					fsc = ptcl.calc_fourier_shell_correlation(projc)
					third = len(fsc)/3
					fsc=array(fsc[third:third*2])
#					snr=fsc/(1.0-fsc)
					sums=[sum(fsc[rings[k]:rings[k+1]])/(rings[k+1]-rings[k]) for k in xrange(4)]		# sum the fsc into 5 range values
#					sums=[sum(snr[rings[k]:rings[k+1]])/(rings[k+1]-rings[k]) for k in xrange(4)]		# sum the fsc into 5 range values
					fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t# {};{}\n".format(sums[0],sums[1],sums[2],sums[3],alt,az,i,defocus,j,cptcl[eo]))
					#xaxis = fsc[0:third]
					#fsc = fsc[third:2*third]
##					saxis = [x/apix for x in xaxis]
##					Util.save_data(saxis[1],saxis[1]-saxis[0],fsc[1:-1],args[1])
					#Util.save_data(xaxis[1],xaxis[1]-xaxis[0],fsc[1:-1],"ptclfsc/f{:04d}_{:1d}_{:06d}.txt".format(i,eo,j))

		print "Results in ptclfsc_{}.txt".format(args[0][-2:])
		sys.exit(0)

	if options.resolution:

		if not os.path.isdir(args[0]):
			print "You must provide the name of the refine_XX folder"
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
		print "Generated : converge_{}.pdf".format(args[0][-2:])
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
				for si in xrange(2,len(d[0])-2):
					if d[1][si-1]>0.143 and d[1][si]<=0.143 :
						frac=(0.143-d[1][si])/(d[1][si-1]-d[1][si])		# 1.0 if 0.143 at si-1, 0.0 if .143 at si
						lastres=d[0][si]*(1.0-frac)+d[0][si-1]*frac
						try:
							plt.annotate(r"{:1.1f} $\AA$".format(1.0/lastres),xy=(lastres,0.143),
								xytext=((lastres*4+d[0][-1])/5.0,0.2),arrowprops={"width":1,"frac":.1,"headwidth":7,"shrink":.05})
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
		print "Generated: goldstandard_{}.pdf".format(args[0][-2:])
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
			itr=max([i for i in os.listdir(r) if "fsc_masked" in i and i[-4:]==".txt"])
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
		print "Generated: goldstandard.pdf"
		plt.clf()


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

			if com in ("e2refine.py","e2refine_easy.py","e2refinemulti.py","e2project3d.py","e2simmx.py","e2simmx2stage.py","e2classify.py","e2classaverage.py","e2make3d.py","e2make3dpar.py","e2refine_postprocess.py") : hist.append((com,spl))
			
		n=0
		while n<len(hist):
			com=hist[n][1][4]
			ttime=timestamp_diff(hist[n][1][0],hist[n][1][1])
			if hist[n][0] in ("e2refine.py","e2refine_easy.py"): 
				pl=com.find("--path=")
				parl=com.find("--parallel=")
				print "%s\t%1.2f hours\te2refine %s"%(difftime(ttime),ttime/3600.0,com[pl+7:].split()[0]),
				if parl>0: print com[parl+11:].split()[0]
				else: print " "

			else:
				print "\t%s\t%1.2f hours\t%s"%(difftime(ttime),ttime/3600.0,hist[n][0])
			
			n+=1
			

		
if __name__ == "__main__":
    main()
