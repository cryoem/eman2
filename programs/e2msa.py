#!/usr/bin/env python
from __future__ import print_function

#
# Author: Steve Ludtke, 06/26/18 (sludtke@bcm.edu) rewrote using sklearn
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
import numpy as np
import sklearn.decomposition as skdc
import time
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <input stack> <output basis>
This too provides a variety of dimensionality reduction methods. This new version
uses scikit.learn, which provides a greater variety of algorithms, but must load 
all data into memory. If working with a large file, you may want to consider using
--step to operate on a limited subset of the data.

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

	parser.add_argument("--mode",type=str,help="Mode should be one of: pca",default="pca")
	parser.add_argument("--nbasis","-n",type=int,help="Number of basis images to generate.",default=20)
	parser.add_argument("--maskfile","-M",type=str,help="File containing a mask defining the pixels to include in the Eigenimages")
	parser.add_argument("--mask",type=int,help="Mask radius, negative values imply ny/2+1+mask, --mask=0 disables, --maskfile overrides",default=0)
	parser.add_argument("--simmx",type=str,help="Will use transformations from simmx on each particle prior to analysis")
	parser.add_argument("--normalize",action="store_true",help="Perform a careful normalization of input images before MSA. Otherwise normalization is not modified until after mean subtraction.",default=False)
	parser.add_argument("--step",type=str,default="0,1",help="Specify <init>,<step>. Processes only a subset of the input data. For example, 0,2 would process only the even numbered particles")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	#parser.add_argument("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	#parser.add_argument("--boxsize","-B",type=int,help="Box size in pixels",default=-1)
	#parser.add_argument("--dbin","-D",type=str,help="Filename to read an existing box database from",default=None)

	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output filenames required")

	logid=E2init(sys.argv,options.ppid)

	if options.verbose>0 : print("Beginning MSA")

	try : step = int(options.step.split(",")[0]),int(options.step.split(",")[1])	# convert strings to tuple
	except:
		print("Invalid --step specification")
		sys.exit(1)


	# setup mask image
	if options.maskfile:
		mask=EMData(options.maskfile,0)
		if mask["mean_nonzero"]!=1.0 :
			print("ERROR: maskfile must be a binary mask (1/0 only)")
			sys.exit(1)
	else :
		# default is no masking
		mask=EMData(args[0],0)
		mask.to_one()
		# negative values handled by mask.sharp
		if options.mask!=0:
			mask.process_inplace("mask.sharp",{"outer_radius":options.mask})
	
	# Memory usage warning >2G raw data
	n=(EMUtil.get_image_count(args[0])-step[0])//step[1]
	nval=mask["square_sum"]
#	print(args[0],n,nval)
	if options.verbose or n*nval>500000000: print("Estimated memory usage (mb): ",n*nval*4/2**20)
	
	# Read all image data into numpy array
	if options.simmx : data=simmx_get(args[0],options.simmx,mask,step)
	else : data=normal_get(args[0],mask,step)
	
	if options.normalize:
		for i in range(len(data)): data[i]/=np.linalg.norm(data[i])

	# first output image is the mean of the input vectors, which has been subtracted from each vector
	mean=np.mean(data,0)
	for i in range(len(data)): data[i]-=mean
	from_numpy(mean).process("misc.mask.pack",{"mask":mask,"unpack":1}).write_image(args[1],0)
	
	# This is where the actual action takes place!
	if options.mode=="pca":
		msa=skdc.PCA(n_components=options.nbasis)
#		print(data.shape)
		msa.fit(data)
		
#	print(msa.components_.shape)
#	c=from_numpy(msa.components_.copy()).write_image("z.hdf",0)

	if options.verbose>0 : print("MSA complete")

	# write other basis vectors
	for i,v in enumerate(msa.components_):
		im=from_numpy(v.copy()).process("misc.mask.pack",{"mask":mask,"unpack":1})
		im["eigval"]=float(msa.singular_values_[i])
		im["explvarfrac"]=float(msa.explained_variance_ratio_[i])
		im.write_image(args[1],i+1)
		if options.verbose : print("Explained variance: ",im["explvarfrac"],"\tSingular Value: ",im["eigval"])

	E2end(logid)

def simmx_get(images,simmxpath,mask,step):
	"""returns an array of transformed masked images as arrays for PCA"""

	simmx=[EMData(simmxpath,i) for i in range(5)]


	n=(EMUtil.get_image_count(images)-step[0])/step[1]

	ret=EMData(int(mask["square_sum"]),n)
	for i in range(n):
		im=EMData(images,i*step[1]+step[0])
		xf=get_xform(i,simmx)
		im.transform(xf)
		imm=im.process("misc.mask.pack",{"mask":mask})
		ret.insert_clip(imm,(0,i,0))
		
	return to_numpy(ret).copy()

def normal_get(images,mask,step):
	"""returns an array of transformed masked images as arrays for PCA"""

	n=(EMUtil.get_image_count(images)-step[0])/step[1]

	ret=EMData(int(mask["square_sum"]),n)
	for i in range(n):
		im=EMData(images,i*step[1]+step[0])
		imm=im.process("misc.mask.pack",{"mask":mask})
		ret.insert_clip(imm,(0,i,0))
	
	return to_numpy(ret).copy()
	
def get_xform(n,simmx):
	"""Will produce a Transform representing the best alignment from a similarity matrix for particle n
	simmx is a list with the 5 images from the simmx file"""

	# find the best orienteation from the similarity matrix, and apply the transformation
	best=(1.0e23,0,0,0,0)

	for j in range(simmx[0].get_xsize()):
		if simmx[0].get(j,n)<best[0] : best=(simmx[0].get(j,n),simmx[1].get(j,n),simmx[2].get(j,n),simmx[3].get(j,n),simmx[4].get(j,n))

	ret=Transform({"type":"2d","alpha":best[3],"tx":best[1],"ty":best[2],"mirror":int(best[4])})

	return ret.inverse()


if __name__== "__main__":
	main()

