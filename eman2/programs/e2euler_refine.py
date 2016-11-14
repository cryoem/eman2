#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/14/2016 (sludtke@bcm.edu),
# Copyright (c) 2000-2008 Baylor College of Medicine
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
from copy import deepcopy
from math import ceil
import os
import sys
import math
import random
import traceback
import numpy as np
from Simplex import Simplex

def get_usage():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options]
	Takes as input/output a set of class-averages or particles with Euler angles in their header, and a 3-D map.
	Optimizes orientation based on agreement with the map. Initial orientations should be close to correct.

	"""
	return usage

def print_usage():

	usage = get_usage()
	print "usage " + usage;
	print "Please run '" + progname + " -h' for detailed options"

def main():
	parser = EMArgumentParser(usage=get_usage())

	parser.add_argument("--input", default=None, help="The input projections. Project should usually have the xform.projection header attribute, which is used for slice insertion")
	parser.add_argument("--ref_volume",default=None, help="The reference volume for orientation determination.")
	parser.add_argument("--input_model", default=None, help="If the class-averages have the model_id parameter (produced by e2refinemulti.py), this will use only class-averages with the specified model_id for the reconstruction.")
	parser.add_argument("--no_wt", action="store_true", dest="no_wt", default=False, help="This argument turns automatic weighting off causing all images to be weighted by 1. If this argument is not specified images inserted into the reconstructed volume are weighted by the number of particles that contributed to them (i.e. as in class averages), which is extracted from the image header (as the ptcl_repr attribute).")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer. This is the only parallelism supported by e2make3dpar", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")

	# Database Metadata storage
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.input==None or options.ref_volume==None :
		print "Error: --input and --ref_volume are both required parameters"
		sys.exit(1)

	if options.input_model!=None : options.input_model=int(options.input_model)

	print "e2euler_refine.py"
	logger=E2init(sys.argv,options.ppid)

	# get basic image parameters, find the first good image
	tmp=EMData()
	tmp.read_image(options.input,0,True)
	if not options.no_wt :
		try:
			n=1
			while tmp["ptcl_repr"]==0 :
				tmp.read_image(options.input,n,True)
				n+=1
		except: pass

	nx=tmp["nx"]
	ny=tmp["ny"]
	apix=tmp["apix_x"]


	if options.verbose>0: print "Image dimensions %d x %d"%(nx,ny)

	data=initialize_data(options.input,options.input_model,options.no_wt)

	ptclcount=sum([i["weight"] for i in data])

	# Read the volume for comparison
	refvol=EMData(options.ref_volume,0)

	#########################################################
	# The actual orientation refinement

	threads=[threading.Thread(target=reorient,args=(data[i::options.threads],refvol,options.verbose-1)) for i in xrange(options.threads)]

	for i,t in enumerate(threads):
		if options.verbose>1: print "started thread ",i
		t.start()

	for t in threads: t.join()

	if options.verbose>0 : print "Finished"

	E2end(logger)

	print "Exiting"

def initialize_data(inputfile,inputmodel,no_weights):
	"""returned list will contain dictionaries containing metadata about each image, and optionally the image itself
	xform - a Transform object
	weight - a floating point weight
	filename - image source filename, not required if data is set
	filenum - image number within filename
	data - an EMData object, if not present or None, filename/num will be used
	nx,ny - size of a single image in a stack file"""

	n_input=EMUtil.get_image_count(inputfile)
	nx,ny,nslice= gimme_image_dimensions3D(inputfile)
	if n_input==1 and nslice>1 and tltfile==None : raise Exception,"Require tlt file to work with volumetric stacks"
	print n_input," input images"

	data=[]

	tmp=EMData()
	for i in xrange(n_input):
		tmp.read_image(inputfile,i,True)
		#else : tmp=get_processed_image(inputfile,i,-1,preprocess,pad)

		# these rely only on the header
		try: elem={"xform":tmp["xform.projection"]}
		except : continue
			#raise Exception,"Image %d doesn't have orientation information in its header"%i

		# skip any particles targeted at a different model
		if inputmodel != None and tmp["model_id"]!=inputmodel : continue

		if no_weights: elem["weight"]=1.0
		else :
			try: elem["weight"]=float(tmp["ptcl_repr"])
			except: elem["weight"]=1.0
			# This is bad if you have actual empty classes...
			#if elem["weight"]<=0 :
				#print "Warning, weight %1.2f on particle %d. Setting to 1.0"%(elem["weight"],i)
				#elem["weight"]=1.0

		try: elem["quality"]=float(tmp["class_qual"])
		except:
			try: elem["quality"]=1.0/(elem["weight"]+.00001)
			except: elem["quality"]=1.0
		elem["filename"]=inputfile
		elem["filenum"]=i
		#if not lowmem:
			#elem["data"]=tmp
			#tmp=EMData()

		data.append(elem)

	print "Using %d images"%len(data)
	return data

def qual(eul,dat):
	xf=Transform({"type":"eman","az":eul[0],"alt":eul[1],"phi":eul[2]})
	prj=dat[0].project("standard",xf)
	ret=dat[1].cmp("frc",prj,{"minres":50.0,"maxres":8.0})
#	print "{:1.4f}\t{:1.2f}\t{:1.2f}\t{:1.2f}".format(ret,eul[0],eul[1],eul[2])
	
	return ret

def reorient(data,refvol,verbose=0):

	output=None		# deletes the results from the previous iteration if any

	if verbose>0:print "Checking Slices (%d)"%len(data)

	ptcl=0
	for i,elem in enumerate(data):
		# get the image to insert
		img=EMData(elem["filename"],elem["filenum"])
		if img["sigma"]==0 : continue

		rd=elem["xform"].get_rotation("eman")
		if verbose>0 : print " %d/%d\r"%(i,len(data)),
		sys.stdout.flush()

		simp=Simplex(qual,[rd["az"],rd["alt"],rd["phi"]],[5,5,5],data=[refvol,img])
		final=simp.minimize(maxiters=500,epsilon=.0001,monitor=0)[0]

		newort=Transform({"type":"eman","az":final[0],"alt":final[1],"phi":final[2]})
		
		img=EMData(elem["filename"],elem["filenum"])
		img["xform.projection.old"]=img["xform.projection"]
		img["xform.projection"]=newort
		
		if verbose>1 : print "\n{:1.2f} {:1.2f} {:1.2f}\t{:1.2f} {:1.2f} {:1.2f} ".format(rd["az"],rd["alt"],rd["phi"],final[0],final[1],final[2])
		img.write_image(elem["filename"],elem["filenum"])

	return

if __name__=="__main__":
	main()
