#!/usr/bin/env python

#
# Author: Wen Jiang, 04/10/2003 (jiang12@purdue.edu)
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
#Beginning MSA
# e2msa.py  01/20/2008  Steven Ludtke
# Rewritten version which just does PCA, no classification
# uses Chao Yang's new PCA implementation in Analyzer 

from EMAN2 import *
from math import *
import time
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <input stack> <output basis>
At present, this program performs only Principle Components Analysis (PCA),
but it may be expanded to include other Multivariate Statistical Analysis (MSA)
algorithms.

---
Performs multivariate statistical analysis on a stack of images. Writes
a set of Eigenimages which can be uses as a basis set for reducing
the dimensionality of a data set (noise reduction). Typically this
basis set is then used to reproject the data (e2basis.py) and
classify the data based on the projected vectors. If the
output file supports arbitrary metadata (like HDF), Eigenvalues
are stored in the 'eigval' parameter in each image.

Note: The mean value is subtracted from each image prior to Eigenimage
calculation. The mean image is stored as the first image in the output
file, though it is not part of the orthonormal basis when
handled this way."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--nbasis","-n",type=int,help="Number of basis images to generate.",default=20)
	parser.add_argument("--maskfile","-M",type=str,help="File containing a mask defining the pixels to include in the Eigenimages")
	parser.add_argument("--mask",type=int,help="Mask radius",default=-2)
	parser.add_argument("--varimax",action="store_true",help="Perform a 'rotation' of the basis set to produce a varimax basis",default=False)
#	parser.add_argument("--lowmem","-L",action="store_true",help="Try to use less memory, with a possible speed penalty",default=False)
	parser.add_argument("--simmx",type=str,help="Will use transformations from simmx on each particle prior to analysis")
	parser.add_argument("--normalize",action="store_true",help="Perform a careful normalization of input images before MSA. Otherwise normalization is not modified until after mean subtraction.",default=False)
	parser.add_argument("--gsl",action="store_true",help="Use gsl SVD algorithm",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	#parser.add_argument("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	#parser.add_argument("--boxsize","-B",type=int,help="Box size in pixels",default=-1)
	#parser.add_argument("--dbin","-D",type=str,help="Filename to read an existing box database from",default=None)
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output filenames required")

	logid=E2init(sys.argv,options.ppid)
	
	if options.verbose>0 : print "Beginning MSA"
	if options.gsl : mode="svd_gsl"
	else : mode="pca_large"
	#elif options.lowmem : mode="pca_large"
	#else : mode="pca"

	if options.mask>-2 :
		mask=EMData(args[0],0)
		ny=mask.get_ysize()
		mask.to_one()
		if options.mask<=1 : mask.process_inplace("mask.sharp",{"outer_radius":ny/2-1})
		else : mask.process_inplace("mask.sharp",{"outer_radius":options.mask})
	else :
		try:
			# try to read in the mask file
			mask=EMData(options.maskfile,0)
		except:
			# default mask is to use all pixels
			mask=EMData(args[0],0)
			mask.to_one()
	
	if options.simmx : out=msa_simmx(args[0],options.simmx,mask,options.nbasis,options.varimax,mode,options.normalize)
	else : out=msa(args[0],mask,options.nbasis,options.varimax,mode,options.normalize)
	
	if options.verbose>0 : print "MSA complete"
	for j,i in enumerate(out):
		try : 
			if options.verbose>0 : print "Eigenvalue: ",i.get_attr("eigval")
		except: pass
		i.write_image(args[1],j)
		
	E2end(logid)

def msa_simmx(images,simmxpath,mask,nbasis,varimax,mode,normalize=True):
	"""Perform principle component analysis (in this context similar to Multivariate Statistical Analysis (MSA) or
Singular Value Decomposition (SVD). 'images' is a filename containing a stack of images to analyze which is coordinated
with simmx. 'simmx' contains the result of an all-vs-all alignment which will be used to transform the orientation
of each image before MSA. 'mask' is an EMImage with a binary mask defining the region to analyze (must be the same size as the input
images. input images will be masked and normalized in-place. The mean value is subtracted from each image prior to
calling the PCA routine. The first returned image is the mean value, and not an Eigenvector. It will have an
'eigval' of 0.  If 'varimax' is set, the final basis set will be 'rotated' to produce a varimax basis. Mode must be one of
pca,pca_large or svd_gsl"""

	simmx=[EMData(simmxpath,i) for i in range(5)]


	n=EMUtil.get_image_count(images)
	if mode=="svd_gsl" : pca=Analyzers.get(mode,{"mask":mask,"nvec":nbasis,"nimg":n})
	else : pca=Analyzers.get(mode,{"mask":mask,"nvec":nbasis})

	mean=EMData(images,0)
	for i in range(1,n):
		im=EMData(images,i)
		xf=get_xform(i,simmx)
		im.transform(xf)
#		im*=mask
		if normalize: im.process_inplace("normalize.unitlen")
#		im.write_image("ali.hdf",i-1)
		mean+=im
	mean.mult(1.0/float(n))
#	mean.mult(mask)
	
	for i in range(n):
		im=EMData(images,i)
		xf=get_xform(i,simmx)
		im.transform(xf)		
		if normalize: im.process_inplace("normalize.toimage",{"to":mean})
		im-=mean
		im*=mask
		im.process_inplace("normalize.unitlen")
		pca.insert_image(im)
			
	mean*=mask
	results=pca.analyze()
	for im in results: im.mult(mask)
	
	if varimax:
		pca=Analyzers.get("varimax",{"mask":mask})
		
		for im in results:
			pca.insert_image(im)
		
		results=pca.analyze()
		for im in results: im.mult(mask)

	for im in results:
		if im["mean"]<0 : im.mult(-1.0)

	mean["eigval"]=0
	mean.process_inplace("normalize.unitlen")
	results.insert(0,mean)
	return results 

def get_xform(n,simmx):
	"""Will produce a Transform representing the best alignment from a similarity matrix for particle n
	simmx is a list with the 5 images from the simmx file"""
	
	# find the best orienteation from the similarity matrix, and apply the transformation
	best=(1.0e23,0,0,0,0)
	
	for j in range(simmx[0].get_xsize()): 
		if simmx[0].get(j,n)<best[0] : best=(simmx[0].get(j,n),simmx[1].get(j,n),simmx[2].get(j,n),simmx[3].get(j,n),simmx[4].get(j,n))
	
	ret=Transform({"type":"2d","alpha":best[3],"tx":best[1],"ty":best[2],"mirror":int(best[4])})
	
	return ret.inverse()

def msa(images,mask,nbasis,varimax,mode,normalize=True):
	"""Perform principal component analysis (in this context similar to Multivariate Statistical Analysis (MSA) or
Singular Value Decomposition (SVD). 'images' is either a list of EMData or a filename containing a stack of images
to analyze. 'mask' is an EMImage with a binary mask defining the region to analyze (must be the same size as the input
images. input images will be masked and normalized in-place. The mean value is subtracted from each image prior to
calling the PCA routine. The first returned image is the mean value, and not an Eigenvector. It will have an
'eigval' of 0.  If 'varimax' is set, the final basis set will be 'rotated' to produce a varimax basis. Mode must be one of
pca,pca_large or svd_gsl"""
	
	
	if isinstance(images,str) :
		n=EMUtil.get_image_count(images)
		if mode=="svd_gsl" : pca=Analyzers.get(mode,{"mask":mask,"nvec":nbasis,"nimg":n})
		else : pca=Analyzers.get(mode,{"mask":mask,"nvec":nbasis})

		mean=EMData(images,0)
		mean.to_zero()
		for i in range(n):
			im=EMData(images,i)
#			im*=mask
			if normalize : im.process_inplace("normalize.unitlen")
			mean+=im
		mean.mult(1.0/float(n))
#		mean.mult(mask)
		
		for i in range(n):
			im=EMData(images,i)
			if normalize: im.process_inplace("normalize.toimage",{"to":mean})
			im-=mean
			im*=mask
			im.process_inplace("normalize.unitlen")
			pca.insert_image(im)
	else:
		n = len(images)
		if mode=="svd_gsl" : pca=Analyzers.get(mode,{"mask":mask,"nvec":nbasis,"nimg":n})
		else : pca=Analyzers.get(mode,{"mask":mask,"nvec":nbasis})
		
		mean=images[0]
		mean.to_zero()
		for im in images:
#			im*=mask
			if normalize : im.process_inplace("normalize.unitlen")
			mean+=im
		mean/=float(n)		
		
		for im in images:
			if normalize: im.process_inplace("normalize.toimage",{"to":mean})
			im-=mean
			im*=mask
			im.process_inplace("normalize.unitlen")
			pca.insert_image(im)
			
	results=pca.analyze()
	for im in results: im.mult(mask)
	
	if varimax:
		pca=Analyzers.get("varimax",{"mask":mask})
		
		for im in results:
			pca.insert_image(im)
		
		results=pca.analyze()
		for im in results: im.mult(mask)

	for im in results:
		if im["mean"]<0 : im.mult(-1.0)

	mean["eigval"]=0
	mean*=mask
	mean.process_inplace("normalize.unitlen")
	results.insert(0,mean)
	return results 

if __name__== "__main__":
	main()
	
