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
from os import system
from os import unlink
from sys import argv
from EMAN2 import *
from optparse import OptionParser

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input stack>
	
Performs k-means classification on a stack of aligned input images. If images are unaligned
one option is to use the --fp option in e2proc2d.py to generate invariants which can then
be classified. """

	parser = OptionParser(usage=usage,version=EMANVERSION)
	parser.add_option("--ncls","-N",type="int",help="Number of classes to generate",default=-1)
	parser.add_option("--average","-A",action="store_true",help="Average the particles within each class",default=False)
	parser.add_option("--clsfiles","-C",action="store_true",help="Write cls files with members of each class",default=False)
	parser.add_option("--listout","-L",action="store_true",help="Output the results to 'class.list",default=False)
	parser.add_option("--nosingle","-X",action="store_true",help="Try to eliminate classes with only 1 member",default=False)
	parser.add_option("--original","-O",type="string",help="If the input stack was derived from another stack, you can provide the name of the original stack here",default=None)

	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	data=EMData.read_images(args[0])

	print len(data)," images to classify."

	Ncls=options.ncls

	# start with Ncls random images
	centers=[]		# the average images for each class
	for i in range(Ncls): centers.append(data[random.randint(0,len(data)-1)])
	
	iter=40
	npcold=[0]*Ncls
	
	while (iter>0) :
		iter-=1
	
		classes=[]					# list of particle #'s in each class
		for i in range(Ncls): classes.append([])
	
		for i in range(len(data)):			# put each particle in a class
			best=(1.0e30,-1)
			for j in range(len(centers)):	# check for best matching center
				c=1.0-centers[j].cmp("dot",data[i],{"normalize":1,"negative":0})
				if (c<best[0]) : best=(c,j)
			classes[best[1]].append((i,best[0]))
	
		# make new averages
		print "\nIteration ",40-iter
		todel=-1
		for j in range(len(centers)):
			print "%3d. %4d\t(%d)"%(j,len(classes[j]),npcold[j])
			if (len(classes[j])==0 ) :
				centers[j]=data[random.randint(0,len(data)-1)]		# reseed empty classes with random image
			elif options.nosingle and len(classes[j])==1:
				centers[j]=data[random.randint(0,len(data)-1)]		# reseed empty classes with random image
				todel=classes[j][0][0]		# delete the particle that was in its own class later
				iter+=1
			else :
				centers[j]=data[classes[j][0][0]].copy()
				for i in range(1,len(classes[j])):
					centers[j]+=data[classes[j][i][0]]
				centers[j].process_inplace("normalize")
				
		if todel!=-1 : del data[todel]
				
		npc=map(lambda x:len(x),classes)		# produces a list with the number of particles in each class
		if (npc==npcold) : break
		npcold=npc
	
	if (options.average != None) :
		if (centers[0].get_zsize()>1) :
			for i in range(len(centers)):
				centers[i].write_image("avg.%04d.mrc"%i,0)
		else:
			for i in range(len(centers)):
				centers[i].write_image("avg.hed",-1)
		
	if (options.clsfiles!=None) :
		os.system("rm -f cls????.lst")
		stackname=argv[0]
		if options.original : stackname=options.original
		for j in range(len(classes)):
			out=open("cls%04d.lst"%j,"w")
			out.write("#LST\n")
			for i in range(len(classes[j])):
				out.write("%d\t%s\n"%(classes[j][i][0],stackname))
			out.close()

if __name__ == "__main__":
	main()

