#!/usr/bin/env python

#
# Author: Steven Ludtke, 1/21/2008 (sludtke@bcm.edu)
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
# e2basis.py  01/21/2008  Steven Ludtke

from EMAN2 import *
from math import *
import time
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <command> <file A> ...
	
Performs various options with basis sets, such as the orthogonal
basis produced by e2msa.py. 

varimax <basis input> <basis output>
	Performs a varimax rotation on an input basis set
	
project <basis input> <image input> <projection output>
	Projects a set of images into the input basis subspace. The default
	is to normalize the individual basis vectors, but not the final resulting
	projection. The projections are stored as a 1-D image stack.
	
projectrot <basis input> <image input> <simmx input> <projection output>
	Same as project, except it will rotate/translate the particles based on the
	best match found in simmx before projection."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--normproj",action="store_true",help="Normalize the projections resulting from 'project', such that the length of each vector is 1",default=False)
	parser.add_argument("--normcomponent",action="store_true",help="Normalize the dot product for each component of the output vector. If the basis spans the input vector, then the projected vector length will be 1, otherwise it will be less than 1.",default=False)
	parser.add_argument("--normalize",type=str,help="Normalize the input images using the named processor. Specify None to disable.",default="normalize.unitlen")
	parser.add_argument("--maskfile","-M",type=str,help="File containing a mask to apply to the particles before normalization", default=None)
	parser.add_argument("--mean1",action="store_true",help="Indicates that the first image in the basis set is actually the mean image, which should be subtracted prior to projection. Output from e2msa requires this flag.")
	parser.add_argument("--recalcmean",action="store_true",help="This will recompute the mean from the input set and subtract before projection. Useful if a different normalization is used than in the original basis file.")
	parser.add_argument("--oneout",action="store_true",help="Output is a single 2-D image rather than a set of 1-D images",default=False)
	parser.add_argument("--nbasis","-n",type=int,help="Will use the first n basis images from the input, excluding the mean if present",default=-1)

	parser.add_argument("--basislist","-z",type=str,help="List of basis vectors to use, comma delimited.",default=None)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	#parser.add_argument("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	#parser.add_argument("--boxsize","-B",type=int,help="Box size in pixels",default=-1)
	#parser.add_argument("--dbin","-D",type=str,help="Filename to read an existing box database from",default=None)
	
	(options, args) = parser.parse_args()
#	if len(args)>0 : parser.error("e2basis.py takes no arguments, only options")

	logid=E2init(sys.argv,options.ppid)
	
	if options.normalize.lower()=="none" : options.normalize=None
	
	# second parameter is always the input basis set
	if options.basislist:
		basislist=options.basislist.split(",")
		basislist=[int(i) for i in basislist]
		basis=EMData.read_images(args[1],basislist)
	else:
		if options.nbasis>1 : basis=EMData.read_images(args[1],range(options.nbasis+1))
		else :basis=EMData.read_images(args[1])
		
	if options.mean1 : 
		mean=basis[0]
		del basis[0]
	else : 
		if options.nbasis>1 : del basis[-1]
		mean=None

	if options.verbose>1 : print "Using %d basis vectors"%len(basis)
	
	if options.maskfile : maskfile=EMData("options.maskfile",0)
	else : maskfile=None
	
	# Project an image stack into a basis subspace
	if args[0]=="project" :
		# normalize the basis vectors to unit length
		for b in basis: b.process_inplace("normalize.unitlen")
		
		# outer loop over images to be projected
		n=EMUtil.get_image_count(args[2])
		
		if options.oneout: 
			proj=EMData(len(basis),n,1)
			
			if options.recalcmean:
				mean=EMData(args[2],0)
				for i in range(1,n):
					im=EMData(args[2],i)
					if options.normalize!=None:
						try: im.process_inplace(options.normalize)
						except: print "Warning: Normalization failed"
					mean+=im
				mean/=float(n)
			
			for i in range(n):
				im=EMData(args[2],i)
				if maskfile : im*=maskfile
				if options.normalize!=None:
					try: im.process_inplace(options.normalize)
					except: print "Warning: Normalization failed"
				if mean : im-=mean
			
				# inner loop over the basis images
				l=0
				for j,b in enumerate(basis):
					proj.set_value_at(j,i,0,im.cmp("dot",b,{"normalize":options.normcomponent,"negative":0}))
					l+=proj[j,i]**2
				
				if options.normproj :
					l=sqrt(l)
					for j in range(len(basis)): proj[j,i]/=l
						
						
			
			proj["isvector"]=1
			proj.write_image(args[3],0)
		else:
			if options.recalcmean:
				mean=EMData(args[2],0)
				for i in range(1,n):
					im=EMData(args[2],i)
					mean+=im
				mean/=float(n)

			for i in range(n):
				im=EMData(args[2],i)
				if maskfile!=None : im*=maskfile
				if options.normalize!=None:
					try: im.process_inplace(options.normalize)
					except: print "Warning: Normalization failed"
				if mean!=None : im-=mean

				proj=EMData(len(basis),1,1)
			
				# inner loop over the basis images
				l=0
				for j,b in enumerate(basis):
					proj.set_value_at(j,0,0,im.cmp("dot",b,{"normalize":options.normcomponent,"negative":0}))
					l+=proj[j,0]**2

				if options.normproj :
					l=sqrt(l)
					for j in range(len(basis)): proj[j,0]/=l
					
				proj["isvector"]=1
				proj.write_image(args[3],i)

	
	# Project rotated images into a basis subspace
	elif args[0]=="projectrot" :
		if options.verbose>1 : print "Entering projectrot routine"
		
		# Just read the whole similarity matrix in, since it generally shouldn't be THAT big
		simmx=EMData(args[3],0)
		simdx=EMData(args[3],1)
		simdy=EMData(args[3],2)
		simda=EMData(args[3],3)
		simflip=EMData(args[3],4)
		
		# outer loop over images to be projected
		n=EMUtil.get_image_count(args[2])
		if options.oneout:
			proj=EMData(len(basis)+4,n,1)
			if options.recalcmean:
				mean=EMData(args[2],0)
				for i in range(1,n):
					im=EMData(args[2],i)
					mean+=im
				mean/=float(n)

			for i in range(n):
				if options.verbose >1 : 
					print "  %5d\r"%i,
					sys.stdout.flush()
				elif options.verbose!=0 and i%100==0:
					print "  %5d\r"%i,
					sys.stdout.flush()
				im=EMData(args[2],i)
				
				# find the best orienteation from the similarity matrix, and apply the transformation
				best=(1.0e23,0,0,0,0)
				
				for j in range(simmx.get_xsize()): 
					if simmx.get(j,i)<best[0] : best=(simmx.get(j,i),simdx.get(j,i),simdy.get(j,i),simda.get(j,i),simflip.get(j,i))
				
#				im.rotate_translate(best[3],0,0,best[1],best[2],0)
				im.transform(Transform({"type":"2d","alpha":best[3],"tx":best[1],"ty":best[2],"mirror":int(best[4])}))
#				print best[3],best[1],best[2],best[4]
#				im.write_image("alib.hdf",i)

				if maskfile!=None : im*=maskfile
				if options.normalize!=None:
					try: im.process_inplace(options.normalize)
					except: print "Warning: Normalization failed"
				if mean!=None : im-=mean
				
				# inner loop over the basis images to generate the components of the projection vector
				l=0
				for j,b in enumerate(basis):
					proj.set_value_at(j+4,i,0,im.cmp("dot",b,{"normalize":options.normcomponent,"negative":0}))
					l+=proj[j+4,i]**2
#					im["dot"]=proj[j+4,i]
#					im.write_image("dbug.hdf",-1)
#					b["dot"]=proj[j+4,i]
#					b.write_image("dbug.hdf",-1)

				if options.normproj :
					l=sqrt(l)
					for j in range(len(basis)): proj[j+4,i]/=l

				proj.set_value_at(0,i,best[1])
				proj.set_value_at(1,i,best[2])
				proj.set_value_at(2,i,best[3])
				proj.set_value_at(3,i,best[4])
				
				proj["isvector"]=1
			proj.write_image(args[4],0)
		else:
			if options.recalcmean:
				mean=EMData(args[2],0)
				for i in range(1,n):
					im=EMData(args[2],i)
					mean+=im
				mean/=float(n)

			for i in range(n):
				if options.verbose >1 : 
					print "  %5d\r"%i,
					sys.stdout.flush()
				elif options.verbose and i%100==0:
					print "  %5d\r"%i,
					sys.stdout.flush()
				im=EMData(args[2],i)
				
				# find the best orienteation from the similarity matrix, and apply the transformation
				best=(1.0e23,0,0,0,0)
				
				for j in range(simmx.get_xsize()): 
					if simmx.get(j,i)<best[0] : best=(simmx.get(j,i),simdx.get(j,i),simdy.get(j,i),simda.get(j,i),simflip.get(j,i))
				
#				im.rotate_translate(best[1],0,0,best[2],best[3],0)
				im.transform(Transform({"type":"2d","alpha":best[3],"tx":best[1],"ty":best[2],"mirror":int(best[4])}))

				if maskfile!=None : im*=maskfile
				if options.normalize!=None:
					try: im.process_inplace(options.normalize)
					except: print "Warning: Normazation failed"
				if mean!=None : im-=mean
				
				proj=EMData(len(basis),1,1)
			
				# inner loop over the basis images to generate the components of the projection vector
				l=0
				for j,b in enumerate(basis):
					proj.set_value_at(j,0,0,im.cmp("dot",b,{"normalize":options.normcomponent,"negative":0}))
					l+=proj[j,0]

				if options.normproj :
					l=sqrt(l)
					for j in range(len(basis)): proj[j,0]/=l
				
				proj.set_attr("ref_da",best[1])
				proj.set_attr("ref_dx",best[2])
				proj.set_attr("ref_dy",best[3])
				proj.set_attr("ref_flip",best[4])
				proj["isvector"]=1
				proj.write_image(args[4],i)
		if options.verbose>1 : print "Projectrot complete"
	
	# Apply the varimax rotation to a set of basis vectors
	elif args[0]=="varimax" :
		mask=basis[0].copy()
		mask.to_one()
		pca=Analyzers.get("varimax",{"mask":mask})
		
		for im in basis:
			pca.insert_image(im)
		
		results=pca.analyze()
		for im in results: im.write_image(args[2],-1)
	else: print "Valid commands are project, varimax and projectrot"
	
	E2end(logid)

if __name__== "__main__":
	main()
	
