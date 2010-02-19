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
from optparse import OptionParser

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <raw particle file> <class mx> <projections> 
	
This will use a volume containing a binary mask and use it to classify a set of particles already refined using e2refine.py.

This is part of a multi-step process for separating ligand bound from ligand free particles, but relies on sufficient
ligand/no-ligand contrast in individual images:
1) refine the entire data set with e2refine.py, such that most classes have at least 10-20 particles in them
2) construct a volume the same size as your reconstruction containing a binary mask, with 1 in the region where the ligand would be
3) run this program 

<raw particle file> should be the same file used in the refinement
<class mx> is one of the classification matrix files from the refinement
<projections> contains the projections used for class mx
<mask> is the 3-D mask file identifying the region of interest
<output file> is a prefix for the classified data"""

	parser = OptionParser(usage=usage,version=EMANVERSION)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_option("--maskfile",type="string",default=None,help="File containing a 3-D binary mask to use for data separation")
	parser.add_option("--ref1",type="string",default=None,help="Rather than using a mask, ref1/ref2 permit using a pair of volumes for classification.")
	parser.add_option("--ref2",type="string",default=None,help="Rather than using a mask, ref1/ref2 permit using a pair of volumes for classification.")
	parser.add_option("--cmp",type="string",help="The name of a 'cmp' to be used in conjunction with ref1/2", default="dot:normalize=1")
	parser.add_option("--process",type="string",default=None,help="A processor to apply to the particle data before classifying")
	parser.add_option("--badgroup",action="store_true",default=False,help="Split the data into 4 groups rather than 2. The extra two groups contain particles more likely to be bad.")
	parser.add_option("--badqualsig",type="float",help="When identifying 'bad' particles, particles with similarities >mean+sigma*badqualsig will be considered bad. Default 0.5", default=.5)
	parser.add_option("--badsepsig",type="float",help="When identifying 'bad' particles, if s1/s2 are the similarities to reference 1/2, then those where |s1-s2| < sigma*badsepsig will be excluded. Default 0.25 ", default=0.25)
	
	#parser.add_option("--ncls","-N",type="int",help="Number of classes to generate",default=-1)
	#parser.add_option("--average","-A",action="store_true",help="Average the particles within each class",default=False)

	(options, args) = parser.parse_args()
	if len(args)<4 : parser.error("Please specify <raw particle file> <class mx> <projections> <mask> <output file> ")
	

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
	
	if options.maskfile:
		try:
			mask=EMData(options.maskfile)
		except:
			print "Error reading 3-D mask file"
			sys.exit(1)
	if options.ref1 and options.ref2 :
		try:
			ref1=EMData(options.ref1)
			ref2=EMData(options.ref2)
		except:
			print "Error reading ref1/ref2"
			sys.exit(1)
	elif not options.maskfile:
		print "You must specify either a maskfile or ref1/ref2"
		sys.exit(1)
		
	logid=E2init(sys.argv)
	
	# now we loop over each class, and assess the masked region for each particle in terms of
	# sigma of the image. Note that we don't have a list of which particle is in each class,
	# but rather a list of which class each particle is in, so we do this a bit inefficiently for now
	out=file("plot.ligand.txt","w")
	statall={}	# keyed by particle number, contains (statm,statr,statr2) for each particle, with Null if the right options weren't specified
	for i in range(nref):
		if options.maskfile : 
			projm=mask.project("standard",{"transform":eulers[i]})
			
		if options.ref1 and options.ref2 :
#			print i,eulers[i]
			proj=ref1.project("standard",{"transform":eulers[i]})
			proj.process_inplace("normalize")
			proj2=ref2.project("standard",{"transform":eulers[i]})
			proj2.process_inplace("normalize")
		

		#p1=[]
		#p2=[]
		#p3=[]
			
		# note that if classmx has multiple classifications for a particle, we use only the first
		simcmp=parsemodopt(options.cmp)
		statn=[]
		statr=[]
		statr2=[]
		statm=[]
		for j in range(nptcl):
			if classmx[0,j]!=i : continue		# only proceed if the particle is in this class
			
			ptcl=EMData(args[0],j)
			if options.process!=None :
				popt=parsemodopt(options.process)
				ptcl.process_inplace(popt[0],popt[1])
			ptclxf=Transform({"type":"2d","alpha":cmxalpha[0,j],"mirror":int(cmxmirror[0,j]),"tx":cmxtx[0,j],"ty":cmxty[0,j]}).inverse()
			
			statn.append(j)
			if options.ref1 and options.ref2 :
				projc=proj.process("math.transform",{"transform":ptclxf})		# we transform the mask projection, not the particle (as in the original classification)
				projc2=proj2.process("math.transform",{"transform":ptclxf})
				cmp1=ptcl.cmp(simcmp[0],projc, simcmp[1])
				cmp2=ptcl.cmp(simcmp[0],projc2,simcmp[1])
				result=cmp1-cmp2
				
				statr.append(result)
				statr2.append((cmp1+cmp2,cmp1-cmp2))
#				statr2.append((cmp1+cmp2,cmp1-cmp2,ptcl["ctf"].defocus))
				
				#print j,result			# number of standard deviations above or below surrounding mean
				#out.write("%d\t%f\n"%(j,result))
				
				# generate a list of aligned particles in each group
				#ali=EMData(args[0],j)
				#ali.process_inplace("math.transform",{"transform":ptclxf.inverse()})
				#if result > .01 : p3.append(ali)
				#elif result <-.01 : p1.append(ali)
				#else : p2.append(ali)
				
			if options.maskfile:
				projmc=projm.process("math.transform",{"transform":ptclxf})		# we transform the mask projection, not the particle (as in the original classification)
				ptcl2=ptcl.copy()
				projmc2=projmc.copy()  # 1-mask
				projmc2.mult(-1.0)
				projmc2.add(1.0)
				
				sig=ptcl["sigma"]
				ptcl.mult(projmc)
				ptcl2.mult(projmc2)
				mean1=ptcl["mean"]/projmc["mean"]
				statm.append(mean1)
#				mean2=ptcl2["mean"]/projc2["mean"]
#				result=(mean1-mean2)/sig
#				print j,result		# number of standard deviations above or below surrounding mean
#				out.write("%d\t%f\n"%(j,result))
				
				# generate a list of aligned particles in each group
				#ali=EMData(args[0],j)
				#ali.process_inplace("math.transform",{"transform":ptclxf.inverse()})
				#if result > 1.5 : p3.append(ali)
				#elif result <.5 : p1.append(ali)
				#else : p2.append(ali)
				
		if options.maskfile :
			# classify particles in this specific orientation
			if len(statm)<20 : continue		# we don't try to classify unless we have at least 5 particles in this orientation
			avg=sum(statm)/len(statm)
			statm=[i-avg for i in statm]
		
			if options.ref1 and options.ref2 :
				for i,j in enumerate(statn):
					if options.verbose>1 : print j,statm[i],statr[i]
					out.write("%f\t%f\n"%(statr[i],statm[i]))
					statall[j]=(statm[i],statr[i],statr2[i])
				
			else :
				for i,j in enumerate(statn):
					if options.verbose>1 : print j,statm[i]
					out.write("%d\t%f\n"%(j,statm[i]))
					statall[j]=(statm[i],None,None)
		
		elif options.ref1 and options.ref2:
			for i,j in enumerate(statn):
				if options.verbose>1 : print j,statr[i]
#				out.write("%d\t%f\n"%(j,statr[i]))
				out.write("%f\t%f\n"%(statr2[i][0],statr2[i][1]))
				statall[j]=(None,statr[i],statr2[i])
			
		if options.verbose==1 :
			print "  %d/%d        \r"%(len(statall),nptcl),
			sys.stdout.flush()
		
		#print len(p1),len(p2),len(p3)
		
		#if len(p1)==0 : p1=[test_image()]
		#if len(p2)==0 : p2=[test_image()]
		#if len(p3)==0 : p3=[test_image()]
		#display((sum(p1)/len(p1),sum(p2)/len(p2),sum(p3)/len(p3),proj))
			
			# for debugging
			#oproj=EMData(args[2],i).process("math.transform",{"transform":ptclxf})
			#ptcl.process_inplace("normalize.toimage",{"to":oproj})
			#display((ptcl,projc,oproj))
			#time.sleep(2)

			#print 

	# Now we do the actual classification using the 'statall' array containing all of the results
	
	# here we compute a sigma from a running average
	#out=file("plotsig.txt","w")
	#out2=file("plotmean.txt","w")
	#s=[(statall[i][2][0],statall[i][1]) for i in statall]		# contains c1+c2 vs c1-c2
	#s.sort()
	#for i in range(50,len(s)-50):
		#sq=0
		#sm=0
		#x=0
		#for j in range(i-50,i+50):
			#sq+=s[j][1]**2
			#sm+=s[j][1]
			#x+=s[j][0]
		#x/=100.0
		#sm/=100.0
##		sq=sq/100.0-sm**2
		#sq=sqrt(sq/100.0)			# deviation from zero, not from mean
		#out .write("%f\t%f\n"%(x,sq))
		#out2.write("%f\t%f\n"%(x,sm))
		
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
						write_particle(args[0],"_ref1_bad",i)		# if the particle was more similar to ref1
						counts[4]+=1
					else : 
						write_particle(args[0],"_ref2_bad",i)						# if the particle was more similar to ref2
						counts[5]+=1
				else :
					if statall[i][1]<0 : 
						write_particle(args[0],"_ref1_good",i)		# if the particle was more similar to ref1
						counts[2]+=1
					else : 
						write_particle(args[0],"_ref2_good",i)						# if the particle was more similar to ref2
						counts[3]+=1


		for i in statall:
			if statall[i][1]<0 : 
				write_particle(args[0],"_ref1",i)		# if the particle was more similar to ref1
				counts[0]+=1
			else : 
				write_particle(args[0],"_ref2",i)						# if the particle was more similar to ref2
				counts[1]+=1
	if options.maskfile:
		for i in statall:
			if statall[i][0]==None : continue
			if statall[i][0]<0
				write_particle(args[0],"_low",i)		# low mask density
				counts[6]+=1
			else : 
				write_particle(args[0],"_high",i)		# high mask density
				counts[7]+=1

	print counts
	
	E2end(logid)

def write_particle(source,postfix,n):
	"""Writes the output file correctly regardless of initial file type. For databases it makes virtual stacks."""
	if source[:4].lower()=="bdb:" :
		indb=db_open_dict(source,ro=True)
		outdb=db_open_dict(source+postfix)
		im=indb.get_header(n)
		im["data_path"]=indb.get_data_path(n)		# this makes the file virtual
		outdb[len(outdb)]=im
	else:
		im=EMData(source,n)
		im.write_image(source[:-4]+postfix+source[-4:],-1)
		
	

if __name__ == "__main__":
	main()
