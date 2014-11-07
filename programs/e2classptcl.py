#!/usr/bin/env python

#
# Author: Steven Ludtke, 09/08/2009 (sludtke@bcm.edu)
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

###	e2classptcl.py	Steven Ludtke	9/8/2009
### Extract particles from class-averages

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
	usage = """prog <classes file> <output> [options]
	
Extracts particles associated with specific class-averages and combines them into a new stack/vstack
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--vstack",action="store_true",help="Will output to a bdb virtual stack instead of copying the image data. Input images must have been BDB for this to work.",default=False)
	parser.add_argument("--classes",type=str,help="Comma separated list of class-numbers to extract particles for",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Please provide a file containing class-averages and an output file")
	try :
		options.classes=[int(i) for i in options.classes.split(",")]
	except:
		print "Please specify --classes with a list of comma-separated class numbers"
		sys.exit(1)
	
	logid=E2init(sys.argv,options.ppid)

	if options.vstack : extract_classav_vstack(args[0],options.classes,args[1],verbose=1)
	else : extract_classav(args[0],options.classes,args[1],verbose=1)

	E2end(logid)

def extract_classav_vstack(inpath,classes,outpath,verbose=0):
	"""Extracts particles from class-averages into a virtual stack. inpath is the path to a file with class-averages,
	outpath is the output BDB, and classes is a list/tuple of class average numbers. Returns the number of images copied"""
	
	vstack=db_open_dict(outpath)		# open output db

	outn=0
	for avn in classes:
		av=EMData(inpath,avn)
		
		try:
			imgsrc=av["class_ptcl_src"]
			imgns=av["class_ptcl_idxs"]
		except:
			raise Exception,"Particle doesn't have source image info (%s,d)"%(inpath,avn)

		# If we are writing to a virtual stack
		try:
			src=db_open_dict(imgsrc)		# this is the source database file
		except:
			raise Exception,"Cannot open source images as BDB (%s)"%imgsrc
		
		for n in imgns:
			# here we make the new entry in the vstack
			d=src.get(n,nodata=1).get_attr_dict()		# get the metadata for the image
			d["data_path"]=src.get_data_path(n)			# this is how we avoid copying the image data
			vstack[outn]=d
			outn+=1

		if verbose>0 : print "Class %d: %d particles"%(avn,len(imgns))
		
	
	if verbose>0 : print "%d total particles written to %s"(outn,outpath)
	
	return outn

def extract_classav(inpath,classes,outpath,verbose=0):
	"""Extracts particles from class-averages into an arbitrary image stack. inpath is the path to a file with class-averages,
	outpath is the output file, and classes is a list/tuple of class average numbers"""

	outn=0
	for avn in classes:
		av=EMData(inpath,avn)
		
		try:
			imgsrc=av["class_ptcl_src"]
			imgns=av["class_ptcl_idxs"]
		except:
			raise Exception,"Particle doesn't have source image info (%s,d)"%(inpath,avn)

		# If we are writing to a virtual stack
		
		for n in imgns:
			im=EMData(imgsrc,n)
			im.write_image(outpath,outn)
			outn+=1
		
		if verbose>0 : print "Class %d: %d particles"%(avn,len(imgns))
		
	
	if verbose>0 : print "%d total particles written to %s"(outn,outpath)



if __name__ == "__main__":
	main()
