#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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

###	e2classifykmeans.py	Steven Ludtke	3/4/2006
### Program for classifying raw 2d or 3d data by kmeans

import os
import sys
import random
import time
import string
import math
import glob
from os import system
from os import unlink
from sys import argv
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <input stack>
	
Performs k-means classification on a stack of aligned input images or invariants generated
from the input images. eg - if images are unaligned one option is to use the --fp option in 
e2proc2d.py to generate invariants which can then be classified. 

This program should not be confused with the k-means processor which provides a mechanism
for grouping pixels in a single map into regions. This program classifies similar images
together."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ncls","-N",type=int,help="Number of classes to generate",default=-1)
	parser.add_argument("--average","-A",action="store_true",help="Average the particles within each class",default=False)
	parser.add_argument("--sigma",action="store_true",help="with --average, this will also produce standard deviation images for each average",default=False)
	parser.add_argument("--onein",action="store_true",help="Read 1-d input images from a single 2-D image (oneout in e2basis.py)",default=False)
	parser.add_argument("--oneinali",action="store_true",help="Read 1-d input images from a single 2-D image where the first 4 elements on each row are da,dx,dy,flip",default=False)
	parser.add_argument("--normavg",action="store_true",help="Normalize averages",default=False)
	parser.add_argument("--clsmx",type=str,default=None,help="Standard EMAN2 output suitable for use with e2classaverage, etc.")
	parser.add_argument("--clsfiles","-C",action="store_true",help="Write EMAN 1 style cls files with members of each class",default=False)
#	parser.add_argument("--listout","-L",action="store_true",help="Output the results to 'class.list",default=False)
	parser.add_argument("--mininclass",type=int,help="Try to eliminate classes with fewer than specified members. Default=2",default=2)
	parser.add_argument("--original","-O",type=str,help="If the input stack was derived from another stack, you can provide the name of the original stack here",default=None)
	parser.add_argument("--exclude", type=str,default=None,help="The named file should contain a set of integers, each representing an image from the input file to exclude.")
	parser.add_argument("--minchange", type=int,default=-1,help="Minimum number of particles that change group before deicding to terminate. Default = len(data)/(#cls*25)")
	parser.add_argument("--fastseed", action="store_true", default=False,help="Will seed the k-means loop quickly, but may produce lest consistent results.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")
	
	print "Classify by k-means"
	logid=E2init(sys.argv, options.ppid)
	if options.onein :
		d=EMData(args[0],0)
		xs=d.get_xsize()
		data=[]
		for i in range(d.get_ysize()):
			data.append(d.get_clip(Region(0,i,xs,1)))
	elif options.oneinali :
		d=EMData(args[0],0)
		xs=d.get_xsize()-4
		data=[]
		for i in range(d.get_ysize()):
			data.append(d.get_clip(Region(4,i,xs,1)))
			data[-1].set_attr("ref_dx",d.get_value_at(0,i))
			data[-1].set_attr("ref_dy",d.get_value_at(1,i))
			data[-1].set_attr("ref_da",d.get_value_at(2,i))
			data[-1].set_attr("ref_flip",d.get_value_at(3,i))
	else: data=EMData.read_images(args[0])
	nimg=len(data)						# we need this for the classification matrix when exclude is used
	filen=range(len(data))				# when exclude is used, this will map to actual file image numbers

	if options.exclude: 
		try:
			excl=file(options.exclude,"r").readlines()
			excl=[int(i) for i in excl]
			excl.sort(reverse=True)
			for i in excl : 
				del data[i]
				del filen[i]
		except: print "Warning: exclude file failed"		# it's ok if this fails

	print len(data)," images to classify."

	if options.minchange<=0 : options.minchange=len(data)/(options.ncls*25)+1
	if options.fastseed : slowseed=0
	else : slowseed=1
	an=Analyzers.get("kmeans")
	an.set_params({"ncls":options.ncls,"minchange":options.minchange,"verbose":1,"slowseed":slowseed,"calcsigmamean":options.sigma,"mininclass":options.mininclass})
	
	an.insert_images_list(data)
	centers=an.analyze()
	
	nrep=[i.get_attr("ptcl_repr") for i in centers[:options.ncls]]
	maxcls=max(nrep)
	for n,i in enumerate(nrep):
		print "%d) %s (%d)"%(n,"#"*int(i*72/maxcls),i)
		
	classes=[[] for i in range(options.ncls)]
	for n,i in enumerate(data):
		classes[i.get_attr("class_id")].append(n)
		
	# This is the old python version of the algorithm, functional but slow
	# left here in case someone needs something they can tweak
	
	## start with Ncls random images
	#centers=[]		# the average images for each class
	#for i in range(Ncls): centers.append(data[random.randint(0,len(data)-1)])
	
	#iter=40
	#npcold=[0]*Ncls
	
	#while (iter>0) :
		#iter-=1
	
		#classes=[]					# list of particle #'s in each class
		#for i in range(Ncls): classes.append([])
	
		#for i in range(len(data)):			# put each particle in a class
			#best=(1.0e30,-1)
			#for j in range(len(centers)):	# check for best matching center
				#c=1.0-centers[j].cmp("dot",data[i],{"normalize":1,"negative":0})
				#if (c<best[0]) : best=(c,j)
			#classes[best[1]].append((i,best[0]))
	
		## make new averages
		#print "\nIteration ",40-iter
		#todel=-1
		#for j in range(len(centers)):
			#print "%3d. %4d\t(%d)"%(j,len(classes[j]),npcold[j])
			#if (len(classes[j])==0 ) :
				#centers[j]=data[random.randint(0,len(data)-1)]		# reseed empty classes with random image
			#elif options.nosingle and len(classes[j])==1:
				#centers[j]=data[random.randint(0,len(data)-1)]		# reseed empty classes with random image
				#todel=classes[j][0][0]		# delete the particle that was in its own class later
				#iter+=1
			#else :
				#centers[j]=data[classes[j][0][0]].copy()
				#for i in range(1,len(classes[j])):
					#centers[j]+=data[classes[j][i][0]]
				#if options.normavg : centers[j].process_inplace("normalize")
				#else: centers[j]/=len(classes[j])-1
				
		#if todel!=-1 : del data[todel]
				
		#npc=map(lambda x:len(x),classes)		# produces a list with the number of particles in each class
		#if (npc==npcold) : break
		#npcold=npc
	
	if (options.average) :
		if (centers[0].get_zsize()>1) :
			for i in range(options.ncls):
				centers[i].write_image("avg.%04d.mrc"%i,0)
				if options.sigma: centers[i+options.ncls].write_image("avgsig.%04d.mrc"%i,0)
		else:
			# write the class-averages to avg.hed
			for i in range(options.ncls):
				centers[i].write_image("avg.hdf",-1)
				if options.sigma: centers[i+options.ncls].write_image("avgsig.hdf",-1)
			
			# if original images specified, also write those averages to avg.orig.hed
			if options.original :
				for j in range(options.ncls):
					avg=EMData(options.original,filen[classes[j][0]])
					if options.sigma: 
						sig=EMData(avg.get_xsize(),avg.get_ysize(),avg.get_zsize())
						sig.to_zero()
						sig.addsquare(EMData(options.original,filen[classes[j][0]]))

					for i in range(1,len(classes[j])):
						avg+=EMData(options.original,filen[classes[j][i]])
						if options.sigma : sig.addsquare(EMData(options.original,filen[classes[j][i]]))

					avg/=len(classes[j])
					avg.write_image("avg.orig.hdf",-1)

					# sigma = sqrt(sumsq / N - avg ^2)/sqrt(N)
					if options.sigma :
						sig.mult(1.0/len(classes[j]))
						sig.subsquare(avg)
						sig.process("math.sqrt")
						sig.mult(1.0/math.sqrt(len(classes[j])))
						sig.write_image("avgsig.orig.hdf",-1)
					
		
	if (options.clsfiles) :
		map(os.remove, glob.glob('cls????.lst'))
		stackname=argv[1]
		if options.original : stackname=options.original
		for j in range(options.ncls):
			out=open("cls%04d.lst"%j,"w")
			out.write("#LST\n")
			for i in range(len(classes[j])):
				out.write("%d\t%s\n"%(filen[classes[j][i]],stackname))
				im=EMData(stackname,filen[classes[j][i]])
				im.write_image("cls%04d.hdf"%j,-1)
			out.close()
	
	# Write an EMAN2 standard classification matrix. Particles run along y
	# each class a particle is in takes a slot in x. There are then a set of
	# 6 images containing class #, a weight, and dx,dy,dangle,flip
	if (options.clsmx) :
		clsnum=EMData(1,nimg,1)
		weight=EMData(1,nimg,1)
		clsnum.to_zero
		clsnum+= -1			# class numbers are initialized to -1 in case we're using exclude
		
		dx=EMData(1,len(data),1)
		dy=EMData(1,len(data),1)
		dang=EMData(1,len(data),1)
		flip=EMData(1,len(data),1)
	
		weight.to_one()
		dx.to_zero()
		for n,i in enumerate(data):
			clsnum[filen[n]]=float(i.get_attr("class_id"))
			try :
				dx[filen[n]]=float(i.get_attr("ref_dx"))
				dy[filen[n]]=float(i.get_attr("ref_dy"))
				dang[filen[n]]=float(i.get_attr("ref_da"))
				flip[filen[n]]=float(i.get_attr("ref_flip"))
			except:
				dx[filen[n]]=0
				dy[filen[n]]=0
				dang[filen[n]]=0
				flip[filen[n]]=0
		
		remove_image(options.clsmx)
		clsnum.write_image(options.clsmx,0)
		weight.write_image(options.clsmx,1)
		dx.write_image(options.clsmx,2)
		dy.write_image(options.clsmx,3)
		dang.write_image(options.clsmx,4)
		flip.write_image(options.clsmx,5)
	
	E2end(logid)

if __name__ == "__main__":
	main()
