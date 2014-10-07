#!/usr/bin/env python

#
# Author: Steven Ludtke, 02/15/2010 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

import os
import sys
import random
import time
import string
import math
from os import system
from os import unlink
from sys import argv
from EMAN2 import *
from EMAN2db import db_open_dict

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <raw particle file> <class mx> <projections>

This program can use either a 3-D binary mask, or a pair of liganded/unliganded volumes to classify particle data into 2
(or 4) groups. e2refine.py must be run on the data first. The method using a pair of volumes is much more accurate in
separating particles.

This is part of a multi-step process for separating ligand bound from ligand free particles, but relies on sufficient
ligand/no-ligand contrast in individual images:
 1) refine the entire data set with e2refine.py, such that most classes have at least 10-20 particles in them
 2) construct a volume the same size as your reconstruction containing a binary mask, with 1 in the region where the ligand would be
2a) -or- prepare 2 volume files representing the particle with or without associated ligand
 3) run this program

<raw particle file> should be the same file used in the refinement
<class mx> is one of the classification matrix files from the refinement
<projections> contains the projections used for class mx

Typical usage:
e2classifyligand.py sets/myset_even.lst refine_01/classmx_04_even.hdf refine_01/projections_04_even.hdf --ref1 ref3d1.hdf --ref2 ref3d2.hdf --cmp=ccc --plotout=cmp.txt --pairmask --splitparticles -v 1
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ref1",type=str,default=None,help="Rather than using a mask, ref1/ref2 permit using a pair of volumes for classification.")
	parser.add_argument("--ref2",type=str,default=None,help="Rather than using a mask, ref1/ref2 permit using a pair of volumes for classification.")
	parser.add_argument("--pairmask",action="store_true",default=False,help="Will use the ref1/ref2 pair to generate a mask which is applied after subtracting ref1 from the particle")
	parser.add_argument("--alistacks",type=float,help="If sum of cmp results is less than the spefified value, will save the aligned particle to a per-class stack",default=-1.0e10)
	parser.add_argument("--cmp",type=str,help="The name of a 'cmp' to be used when pairmask is not specified", default="ccc")
	parser.add_argument("--process",type=str,default=None,help="A processor to apply to the particle data before classifying")
	parser.add_argument("--plotout",type=str,default="plot_ligand.txt",help="Name of a text file for the classification plot.")
	parser.add_argument("--badgroup",action="store_true",default=False,help="Split the data into 4 groups rather than 2. The extra two groups contain particles more likely to be bad.")
	parser.add_argument("--badqualsig",type=float,help="When identifying 'bad' particles, particles with similarities >mean+sigma*badqualsig will be considered bad. Default 0.5", default=.5)
	parser.add_argument("--badsepsig",type=float,help="When identifying 'bad' particles, if s1/s2 are the similarities to reference 1/2, then those where |s1-s2| < sigma*badsepsig will be excluded. Default 0.25 ", default=0.25)
	parser.add_argument("--postfix",type=str,default="",help="This string will be appended to each set name to help differentiate the results from multiple runs")
	parser.add_argument("--splitparticles",action="store_true",default=False,help="Specify this to write new files containing the classified particles")
	parser.add_argument("--tstcls",type=int,help="Will generate tst.hdf containing test images for a specified class-number",default=-1)
	parser.add_argument("--debug",action="store_true",default=False,help="Enable debugging mode with verbose output and image display. Not suitable for real runs.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	#parser.add_argument("--ncls","-N",type=int,help="Number of classes to generate",default=-1)
	#parser.add_argument("--average","-A",action="store_true",help="Average the particles within each class",default=False)

	(options, args) = parser.parse_args()
#	if len(args)<4 : parser.error("Please specify <raw particle file> <class mx> <projections> <mask> <output file> ")


	try:
		classmx=EMData(args[1],0)
		nptcl=classmx["ny"]
		cmxtx=EMData(args[1],2)
		cmxty=EMData(args[1],3)
		cmxalpha=EMData(args[1],4)
		cmxmirror=EMData(args[1],5)
	except:
		print "Error reading classification matrix. Must be full classification matrix with alignments"
		sys.exit(1)

	# this reads all of the EMData headers from the projections
	nref=EMUtil.get_image_count(args[2])
	eulers=[EMData(args[2],i,1)["xform.projection"] for i in range(nref)]

	if not options.ref1 or not options.ref2:
		print "Ref1 anf Ref2 are required arguments"
		sys.exit(1)

	try:
		ref1=EMData(options.ref1)
		ref2=EMData(options.ref2)
		ref2.process_inplace("normalize.toimage",{"to":ref1})
		softmask=EMData(ref1["nx"],ref2["ny"],1)
		softmask.to_one()
		softmask.process_inplace("mask.soft",{"outer_radius":ref1["nx"]/2-4,"width":3})
	except:
		print "Error reading ref1/ref2"
		sys.exit(1)

	logid=E2init(sys.argv, options.ppid)

	# now we loop over each class, and assess the masked region for each particle in terms of
	# sigma of the image. Note that we don't have a list of which particle is in each class,
	# but rather a list of which class each particle is in, so we do this a bit inefficiently for now
	out=file(options.plotout,"w")
	statall={}	# keyed by particle number, contains (statm,statr,statr2) for each particle, with Null if the right options weren't specified
	for i in range(nref):
		if options.tstcls>=0 and i!=options.tstcls : continue
		if options.verbose>1 : print "--- Class %d"%i

		proj=ref1.project("standard",{"transform":eulers[i]})
		proj.process_inplace("normalize.circlemean")
		proj.mult(softmask)
		proj2=ref2.project("standard",{"transform":eulers[i]})
		proj2.process_inplace("normalize.circlemean")
		proj2.mult(softmask)

		if options.pairmask :
			projd=proj2-proj	# should already be normalized	TODO: see if renormalizing in 2d helps ?
			thr=projd["mean"]+projd["sigma"]*2.0		# find a threshold for "important" differences
			projm=projd.process("threshold.binary",{"value":thr})		# a binary mask, hopefully with just the "important parts" of the difference
			projm.process_inplace("mask.addshells",{"nshells":int(projm["nx"]/50)})
			projd.mult(projm)

			# If we normalize to one projection or the other it will bias the values towards that projection, so we average them
			projref=proj+proj2
			projref.mult(0.5)

			# projection copies with the mask we found above
			projc=proj.copy()
			projc.mult(projm)
			projc2=proj2.copy()
			projc2.mult(projm)



		# Computes the union of non-zero points in both projections
		projmask=proj.process("threshold.notzero")
		projmask.add(proj2.process("threshold.notzero"))
		projmask.process_inplace("threshold.notzero")

		#p1=[]
		#p2=[]
		#p3=[]

		# note that if classmx has multiple classifications for a particle, we use only the first
		simcmp=parsemodopt(options.cmp)
		statn=[]
		statr=[]
		statr2=[]
		statm=[]
		nalis=0
		for j in range(nptcl):
			if classmx[0,j]!=i : continue		# only proceed if the particle is in this class

			try: ptcl=EMData(args[0],j)
			except: 
				print "Cannot read particle {} from {}. This should not happen. Using correct input file?".format(j,args[0])
				continue
			if options.process!=None :
				popt=parsemodopt(options.process)
				ptcl.process_inplace(popt[0],popt[1])
			ptclxf=Transform({"type":"2d","alpha":cmxalpha[0,j],"mirror":int(cmxmirror[0,j]),"tx":cmxtx[0,j],"ty":cmxty[0,j]})

			statn.append(j)

			# In this method we look at the pair of projections and identify an appropriate 2-D mask image to use for classification purposes
			if options.pairmask :
				# TODO : maybe expand the mask a bit ?
				#projm.process_inplace("xform",{"transform":ptclxf})

				# While the original classification moved the reference, not the particle, in this case
				# doing that proves not to work very well due to edge artifacts
				ptcl2=ptcl.process("normalize.edgemean")
				ptcl2.process_inplace("xform",{"transform":ptclxf})
				ptcl2.mult(softmask)
#				ptcl2.process_inplace("filter.lowpass.gauss",{"cutoff_freq":0.05})
				ptcl2.process_inplace("filter.matchto",{"to":projref})
##				ptcl2.process_inplace("normalize.toimage",{"to":proj,"high_threshold":proj["mean"]+proj["sigma"]*4.0,"low_threshold":proj["mean"]+proj["sigma"]})		# worked with the old normalize.toimage
				#ptcl2.process_inplace("normalize.toimage",{"to":proj,"high_threshold":proj["mean"]+proj["sigma"]*4.0,"low_threshold":proj["mean"]-proj["sigma"]/2.0})	# test for the new normalize.toimage
				#ptcl3=ptcl2.copy()
				#ptcl3.sub(proj)
				ptcl3=ptcl2.process("math.sub.optimal",{"ref":proj})		# use the new Fourier subtract instead of the rigamarole above


				if options.tstcls==i :
					proj["xform.align2d"]=ptclxf			# just so it shows up in the display options
					proj.write_image("tst.hdf",-1)
					proj2.write_image("tst.hdf",-1)
					ptcl.write_image("tst.hdf",-1)
					ptcl2["xform.align2d"]=ptclxf
					ptcl2.write_image("tst.hdf",-1)
					ptcl3.write_image("tst.hdf",-1)
					projm.write_image("tst.hdf",-1)
					try: sump.add(ptcl2)
					except: sump=ptcl2.copy()

				cmp1=ptcl2.cmp(simcmp[0],proj, simcmp[1])
				cmp2=ptcl2.cmp(simcmp[0],proj2,simcmp[1])
				ptcl3.mult(projm)
				cmp3=ptcl3.cmp("ccc",projd)
				cmp4=ptcl3.cmp(simcmp[0],projc, simcmp[1])
				cmp5=ptcl3.cmp(simcmp[0],projc2,simcmp[1])

				if options.tstcls==i :
					ptcl2=ptcl.process("xform",{"transform":ptclxf.inverse()})
					ptcl2.mult(projmask)
					if cmp5<cmp4:
						try: avgim1.add(ptcl2)
						except: avgim1=ptcl2
					else:
						try: avgim2.add(ptcl2)
						except: avgim2=ptcl2

				if options.verbose>1 : print j,cmp1+cmp2,cmp2-cmp1,cmp3,cmp5-cmp4
				statr.append(cmp3)
				statr2.append((cmp1+cmp2,cmp3,cmp4-cmp5))
#				statr2.append((cmp1+cmp2,cmp2-cmp1,cmp5-cmp4)

			else:
				ptcl2=ptcl.process("filter.matchto",{"to":proj+proj2})
				projc=proj.process("xform",{"transform":ptclxf})		# we transform the projection, not the particle (as in the original classification)
				projc2=proj2.process("xform",{"transform":ptclxf})
				projc.process_inplace("normalize.toimage",{"to":ptcl2})
				projc2.process_inplace("normalize.toimage",{"to":ptcl2})
				cmp1=ptcl2.cmp(simcmp[0],projc, simcmp[1])
				cmp2=ptcl2.cmp(simcmp[0],projc2,simcmp[1])
				result=cmp1-cmp2
				if options.debug: display((ptcl2,projc,projc2,proj))

				statr.append(result)
				statr2.append((cmp1+cmp2,cmp1-cmp2,0))

				if cmp1+cmp2<options.alistacks :
					ptcl2=ptcl.process("xform",{"transform":ptclxf.inverse()})
					ptcl2.mult(projmask)
					ptcl2.write_image("aligned_{}.hdf".format(i),nalis+2)

					if nalis==0 :
						projc=proj.process("normalize.toimage",{"to":ptcl2})
						projc2=proj2.process("normalize.toimage",{"to":ptcl2})
						projc.write_image("aligned_{}.hdf".format(i),0)
						projc2.write_image("aligned_{}.hdf".format(i),1)

					nalis+=1

				if options.tstcls==i :
					ptcl2=ptcl.process("xform",{"transform":ptclxf.inverse()})
					ptcl2.mult(projmask)
					if cmp1>cmp2:
						try: avgim1.add(ptcl2)
						except: avgim1=ptcl2
					else:
						try: avgim2.add(ptcl2)
						except: avgim2=ptcl2


		if options.tstcls==i :
			try:
				avgim1.process_inplace("normalize.edgemean")
				avgim1.write_image("tst2.hdf",0)
			except: pass
			try:
				avgim2.process_inplace("normalize.edgemean")
				avgim2.write_image("tst2.hdf",1)
			except: pass

			try: sump.write_image("tst3.hdf",0)
			except: pass


		if len(statr)==0 and len(statm)==0 :
			if options.verbose>1 : print "No particles"
			continue


		for i,j in enumerate(statn):
			if options.verbose>1 : print j,statr[i]
			out.write("%f\t%f\t%f\t%d\n"%(statr2[i][0],statr2[i][1],statr2[i][2],j))
			statall[j]=(None,statr[i],statr2[i])

		if options.verbose==1 :
			print "  %d/%d        \r"%(len(statall),nptcl),
			sys.stdout.flush()


	# Now we do the actual classification using the 'statall' array containing all of the results

	if options.splitparticles :
		counts=[0,0,0,0,0,0,0,0]		# image counts, all class 1, all class 2,good class 1, good class 2, bad class 1, bad class2, low mask, high mask
		if options.ref1 and options.ref2 :
			if options.badgroup:
				# These aren't very efficient, but shouldn't matter much
				avgq=  sum([statall[i][2][0] for i in statall])/float(len(statall))
				avgqsq=sqrt(sum([statall[i][2][0]**2 for i in statall])/float(len(statall))-avgq**2)
				avgcsq=sqrt(sum([statall[i][1]**2 for i in statall])/float(len(statall)))
				for i in statall:
					# we consider a particle 'bad' if it's mean quality with both refs is greater than the mean+sigma/2
					# or if the quality difference between images is within 1/4 sigma of zero
					if statall[i][2][0]>avgq+avgqsq*options.badqualsig or fabs(statall[i][1])<avgcsq*options.badsepsig :
						if statall[i][1]<0 :
							write_particle(args[0],"_ref1_bad"+options.postfix,i)		# if the particle was more similar to ref1
							counts[4]+=1
						else :
							write_particle(args[0],"_ref2_bad"+options.postfix,i)						# if the particle was more similar to ref2
							counts[5]+=1
					else :
						if statall[i][1]<0 :
							write_particle(args[0],"_ref1_good"+options.postfix,i)		# if the particle was more similar to ref1
							counts[2]+=1
						else :
							write_particle(args[0],"_ref2_good"+options.postfix,i)						# if the particle was more similar to ref2
							counts[3]+=1


			for i in statall:
				if statall[i][1]<0 :
					write_particle(args[0],"_ref1"+options.postfix,i)		# if the particle was more similar to ref1
					counts[0]+=1
				else :
					write_particle(args[0],"_ref2"+options.postfix,i)		# if the particle was more similar to ref2
					counts[1]+=1
		#if options.maskfile:
			#for i in statall:
				#if statall[i][0]==None : continue
				#if statall[i][0]<0:
					#write_particle(args[0],"_low"+options.postfix,i)		# low mask density
					#counts[6]+=1
				#else :
					#write_particle(args[0],"_high"+options.postfix,i)		# high mask density
					#counts[7]+=1

		print counts

	E2end(logid)

# silly hack
glob_inls=None
glob_outls={}

def write_particle(source,postfix,n):
	"""Writes the output file correctly regardless of initial file type. For databases it makes virtual stacks."""
	if source[-4:].lower()==".lst" :
		global glob_inls,glob_outls
		
		if glob_inls==None:
			glob_inls=LSXFile(source)
			
		if not glob_outls.has_key(postfix):
			glob_outls[postfix]=LSXFile(source[:-4]+postfix+".lst")
		
		ent=glob_inls.read(n)
		glob_outls[postfix].write(-1,ent[0],ent[1],ent[2])
	else:
		im=EMData(source,n)
		im.write_image(source[:-4]+postfix+source[-4:],-1)



if __name__ == "__main__":
	main()
