#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
# align all particles to reference and store alignment results
# Author: Steven Ludtke (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
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

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range

from EMAN2 import *
import numpy as np
import matplotlib.pyplot as plt
import sklearn.decomposition as skdc
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.signal import find_peaks_cwt
import time
import os
from sys import argv,exit

def mkpath(base):
	num = 0
	while os.path.isdir("{}_{:02d}".format(base,num)):
		num += 1
	os.mkdir("{}_{:02d}".format(base,num))

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_pcasplit.py [options] <spt_XX> <reference>"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--iter",type=int,required=True,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--path",type=str,required=True,default=None,help="Path to a folder where results should be stored, following standard naming conventions (default = spt_XX)")
	parser.add_argument("--sym",type=str,help="Apply this symmetry.",default=None)
	parser.add_argument("--res",type=float,help="Filter particles to this resolution (in Angstroms) before classification",default=30.0)
	#parser.add_argument("--outpath",type=str,help="",default=None)
	parser.add_argument("--mask",type=str,help="Apply this symmetry.",default=None)
	parser.add_argument("--nowedgefill",action='store_true',help="Do not fill the missing wedge before classification.",default=False)
	parser.add_argument("--keepthresh",type=float,help="Center PCA outliers beyond this value before performing K-means clustering. Default is 0.2.",default=0.2)
	parser.add_argument("--nbasis",type=int,required=True,help="Number of PCA basis vectors. Default is 4.",default=4)
	parser.add_argument("--nclass",type=int,required=True,help="Number of classes. Default is 2.",default=2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()

	if options.path == None:
		print("You must specify the path to an existing spt_XX refinement directory.")
		sys.exit(1)

	if options.iter<=0 :
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print("No subtomogram alignment was complered for this refinement. Please try another")
			sys.exit(1)
		options.iter=max(fls)

	options.outpath = mkpath("sptcls")

	logid=E2init(sys.argv, options.ppid)

	threed = EMData("{}/threed_{:.02d}.hdf".format(options.path,options.iter))
	
	if options.mask: msk = EMData(options.msk)
	else: msk = EMData("{}/mask_tight.hdf".format(options.path))

	refparms=js_open_dict("{}_{:02d}/0_spt_params.json".format(options.path,options.iter))
	inptcls = refparms["input_ptcls"]

	parmsfile = "{}/particle_parms_{:02d}.json".format(options.path,options.iter)
	js=js_open_dict(parmsfile)

	nptcl=len(js.keys())
	pname=eval(js.keys()[0])[0]
	xfs=[]
	coord=[]
	srcs=[]

	refft = threed.do_fft()
	imgs00 = []
	for i in range(nptcl):

		if options.verbose:
			sys.stdout.write("\r{}/{} particles".format(i+1,nptcl))
			sys.stdout.flush()

		k="('{}', {})".format(pname, i)
		xf=js[k]['xform.align3d']
		e=EMData(pname, i)
		e.transform(xf)
		e["score"]=js[k]['score']
		
		if options.nowedgefill: enew = e
		else:
			eft=e.do_fft()
			eft.process_inplace("mask.wedgefill",{"fillsource":refft, "thresh_sigma":0.0})
			enew=eft.do_ift()
		
		enew.mult(msk)
		enew.process_inplace("math.meanshrink",{"n":2})
		enew.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1/options.res})

		if options.sym != 'c1':
			enew.process_inplace("xform.applysym",{"averager":"mean.tomo", "sym":options.sym})
		enew["xform.align3d"]=xf

		enew.write_image("{}/aliptcls.hdf".format(options.outpath), i)
		imgs00.append(enew)

	js.close()

	#imgs00=EMData.read_images("{}/aliptcls.hdf".format(options.outpath))
	nptcl = len(imgs00)

	imgsnp=np.array([m.numpy().copy() for m in imgs00])

	pca=PCA(options.nbasis)
	pout=pca.fit_transform(imgsnp.reshape((nptcl, -1)))

	pout[abs(pout)>options.keepthresh]=0

	basisfile = "{}/pca_basis.hdf".format(options.outpath)
	for i,c in enumerate(pca.components_):
		egmap=c.reshape(imgsnp[0].shape)
		m=from_numpy(egmap.copy())
		m.write_image(basisfile,i)

	kmeans = KMeans(n_clusters=options.ncls).fit(pout)
	lb=kmeans.labels_

	if options.verbose: print("Class: Particle count")
	for i in range(kmeans.n_clusters):
		if options.verbose: print("{}: {}".format(lb,np.sum(lb==i)))

	# subtomogram average particles from each class
	avgr=[Averagers.get("mean.tomo") for i in range(kmeans.n_clusters)]
	for i in range(nptcl):
		avgr[lb[i]].add_image(imgs00[i])

	classfile = "{}/classes3d.hdf".format(options.outpath)
	for i in range(kmeans.n_clusters):
		e=avgr[i].finish()
		e.process_inplace("normalize")
		e.write_image(classfile,i)

	inlst=LSXFile(inptcls, True)

	outlsts=[]
	for lbl in sorted(np.unique(lb)):
		outlst = LSXFile(inptcls.replace(".lst","_{}.lst".format(lbl)))
		for i in range(inlst.n):
			if lb[i]==lbl:
				l=inlst.read(i)
				outlst.write(-1,l[0], l[1],l[2])
		outlst.close()
	inlst.close()

	E2end(logid)


if __name__ == "__main__":
	main()

