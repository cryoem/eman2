#!/usr/bin/env python
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

from future import standard_library
standard_library.install_aliases()
from builtins import range

from EMAN2 import *

from EMAN2_utils import make_path
import numpy as np
import matplotlib.pyplot as plt
import sklearn.decomposition as skdc
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy.signal import find_peaks_cwt
import time
import os
from sys import argv,exit

#def mkpath(base):
	#num = 0
	#while os.path.isdir("{}_{:02d}".format(base,num)):
		#num += 1
	#os.mkdir("{}_{:02d}".format(base,num))
	#return "{}_{:02d}".format(base,num)

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: e2spt_pcasplit.py [options] <spt_XX> <reference>"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--path",type=str,required=True,default=None,help="Path to a folder where results should be stored, following standard naming conventions (default = spt_XX)",guitype='filebox',row=0, col=0, rowspan=1, colspan=2)
	parser.add_argument("--iter",type=int,required=True,help="Iteration number within path. Default is the second to last iteration (-2).",default=-2,guitype='intbox',row=1, col=0, rowspan=1, colspan=1)
	parser.add_argument("--nclass",type=int,required=True,help="Number of classes. Default is 2.",default=2,guitype="intbox",row=1, col=1, rowspan=1, colspan=1)
	parser.add_header(name="orblock1", help='', title="Optional", row=2, col=0, rowspan=1, colspan=2)
	parser.add_argument("--maxres",type=float,help="Filter particles to this resolution (in Angstroms) before classification",default=30.0,guitype="floatbox",row=3, col=0, rowspan=1, colspan=1)
	parser.add_argument("--sym",type=str,help="Apply this symmetry.",default="c1",guitype="strbox",row=3, col=1, rowspan=1, colspan=1)
	parser.add_argument("--mask",type=str,help="Apply this mask. Default is 'mask_tight.hdf' from <--path>_<--iter>. Specify 'none' for no masking",default="",guitype="filebox",row=4, col=0, rowspan=1, colspan=2)
	parser.add_argument("--nbasis",type=int,required=True,help="Number of PCA basis vectors. Default is 3.",default=3,guitype="intbox",row=5, col=0, rowspan=1, colspan=1)
	# parser.add_argument("--keepthresh",type=float,help="Center PCA outliers beyond this value before performing K-means clustering. Default is no threshold (-1).",default=,guitype="floatbox",row=5, col=0, rowspan=1, colspan=1)
	parser.add_argument("--nowedgefill",action='store_true',help="Do not fill the missing wedge before classification.",default=False,guitype="boolbox",row=5, col=1, rowspan=1, colspan=1)
	parser.add_argument("--clean",action='store_true',help="remove outliers before PCA.",default=False,guitype="boolbox",row=6, col=1, rowspan=1, colspan=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--shrink",type=int,help="Shrink particles before classification",default=1)
	parser.add_argument("--dotest",type=int,help="test using N random particles",default=-1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	(options, args) = parser.parse_args()

	if options.path == None:
		print("You must specify the path to an existing spt_XX refinement directory.")
		sys.exit(1)
	if not os.path.isdir(options.path):
		print("Could not locate --path {}".format(options.path))
		sys.exit(1)

	if options.iter==-2:
		fls=[int(i[15:17]) for i in os.listdir(options.path) if i[:15]=="particle_parms_" and str.isdigit(i[15:17])]
		if len(fls)==0 : 
			print("No subtomogram alignment was complered for this refinement. Please try another")
			sys.exit(1)
		options.iter=max(fls)-1

	
	logid=E2init(sys.argv, options.ppid)

	threed = EMData("{}/threed_{:02d}.hdf".format(options.path,options.iter))
	 
	if options.mask=="none":
		msk=threed.copy()
		msk.to_one()
	elif options.mask=="":
		msk = EMData("{}/mask_tight.hdf".format(options.path))
	else:
		msk = EMData(options.mask)

	#refparms=js_open_dict("{}/0_spt_params.json".format(options.path))
	#try:
		#inptcls = refparms["input_ptcls"]
	#except:
		#print("Could not locate input particles. Path and iteration must reference a standard 3D refinement (e2spt_refine.py) and not a subtilt refinement (e2spt_subtilt.py).")
		#sys.exit(1)

	parmsfile = "{}/particle_parms_{:02d}.json".format(options.path,options.iter)
	js=js_open_dict(parmsfile)

	nptcl=len(js.keys())
	irange=np.arange(nptcl, dtype=int)
	if options.dotest>0:
		np.random.shuffle(irange)
		nptcl=min(nptcl, options.dotest)
		irange=irange[:nptcl]
		print("Test with {} particles".format(nptcl))
	inptcls=pname=eval(js.keys()[0])[0]
	

	
	print("Preprocessing {} particles...".format(nptcl))
	
	#### doing everything in fourier space with numpy
	
	data=[]
	wgs=[]
	keys=[]
	
	### to clip fourier space based on resolution
	sz=threed["nx"]//options.shrink
	apix=threed["apix_x"]*options.shrink
	freq=np.fft.fftfreq(sz, apix)[:sz//2]
	clip=sz//2-np.argmax(freq>1./options.maxres)
	if clip==sz//2: clip=0
	#n=EMUtil.get_image_count(pname)
	for i in irange.tolist():
		
		k="('{}', {})".format(pname, i)
		if js.has_key(k)==False: continue
		xf=js[k]['xform.align3d']
		e=EMData(pname, i)
		e.transform(xf)
		e.mult(msk)
		if options.shrink>1:
			e.process_inplace("math.meanshrink",{"n":options.shrink})
	
		en=e.numpy().copy()
		#print(en.shape, clip)
		#### numpy fft is actually significantly slower than EMAN fft. so probably should change to EMAN if I can get the coordinates right..
		
		ef=np.fft.fftshift(np.fft.fftn(en))
		if clip>0:
			ef=ef[clip:-clip,clip:-clip,clip:-clip]
		if len(data)==0:
			sz=len(ef)
			idx=np.indices((sz,sz,sz))-sz//2
			r=np.sqrt(np.sum(idx**2, 0))
			r = r.astype(np.int)
			nr = np.bincount(r.ravel())
			
		efa=abs(ef)
		tbin = np.bincount(r.ravel(), efa.ravel())
		sf = tbin / nr
		sf[sf==0]=1e-5
		div=np.divide(efa,sf[r])
		wdg=np.logical_and(div<1., r>1)
		ef[wdg]=0
		data.append(ef.flatten())
		wgs.append(wdg.flatten())
		keys.append(k)
		
		sys.stdout.write("\r   {}/{} particles".format(len(data),nptcl))
		sys.stdout.flush()
	print()

	js.close()
	data=np.array(data)
	wgs=np.array(wgs)
	print(data.shape)

	#avg=np.mean(data, axis=0)
	avg=np.sum(data, axis=0)
	w=np.sum(1-np.array(wgs), axis=0)+1
	avg=avg/w
	dv=data-avg
	std=np.std(abs(dv))
	dv/=std
	#for i,w in enumerate(wgs):
		#data[i][w]=avg[w]
	#dv=data
	imgsnp=np.hstack([dv.real, dv.imag])
	print(dv.real.shape, imgsnp.shape)
	
	options.outpath = make_path("sptcls")
	print("Output will be written to {}".format(options.outpath))

	#### real space version:
	
	#refft = threed.do_fft()
	#imgs00 = []
	#for i in range(nptcl):

		#sys.stdout.write("\r{}/{} particles".format(i+1,nptcl))
		#sys.stdout.flush()

		#k="('{}', {})".format(pname, i)
		#xf=js[k]['xform.align3d']

		#e=EMData(pname, i)		
		#e.transform(xf)
		#e["score"]=js[k]['score']
		
		#if options.nowedgefill: 
			#enew = e
		#else:
			#eft=e.do_fft()
			#eft.process_inplace("mask.wedgefill",{"fillsource":refft, "thresh_sigma":0.0})
			#enew=eft.do_ift()
		
		#enew.mult(msk)
		#enew.process_inplace("math.meanshrink",{"n":options.shrink})
		#enew.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1/options.maxres})

		#if options.sym != None:
			#enew.process_inplace("xform.applysym",{"averager":"mean.tomo", "sym":options.sym})

		#enew["xform.align3d"]=xf

		#enew.write_image("{}/aliptcls.hdf".format(options.outpath), i)
		#imgs00.append(enew)

	#js.close()
	

	#nptcl = len(imgs00)
	#imgsnp=np.array([m.numpy().copy() for m in imgs00])
	
	print("Performing PCA...")
	#print(imgsnp)
	ptclids=np.arange(nptcl,dtype=int)[:, None]
	nptcl=len(imgsnp)
	
	if options.clean:
		#### do pca twice to remove outliers	
		pca=PCA(options.nbasis)
		pout=pca.fit_transform(imgsnp)
		dst=np.linalg.norm(pout-np.mean(pout, 0), axis=1)
		outlr=dst>np.mean(dst)+np.std(dst)*2
		
		np.savetxt("{}/pca_rmoutlier.txt".format(options.outpath), 
			np.hstack([ptclids, pout]))
		print("Removing {} outliers...".format(np.sum(outlr)))
		
		imgsnp=imgsnp[outlr==False]
		ptclids=ptclids[outlr==False]
	
	
	pca=PCA(options.nbasis)
	pout=pca.fit_transform(imgsnp)
	np.savetxt("{}/pca_ptcls.txt".format(options.outpath), 
		np.hstack([ptclids, pout]))

	basisfile = "{}/pca_basis.hdf".format(options.outpath)
	#threed.process("math.meanshrink",{"n":options.shrink}).write_image(basisfile, 0)
	l=len(data[0])
	for i,c in enumerate(pca.components_):
		eg=c[:l]+c[l:]*1j
		egmap=eg.reshape((sz,sz,sz))
		o=np.real(np.fft.ifftn(np.fft.ifftshift(egmap)))
		m=from_numpy(o.copy())
		m.write_image(basisfile,i)

	print("Classifying particles...")
	kmeans = KMeans(n_clusters=options.nclass).fit(pout)
	lb=kmeans.labels_
	
	if options.clean:
		lbfull=np.zeros(nptcl, dtype=int)
		lbfull[outlr==False]=lb
		lbfull[outlr]=-1
		lb=lbfull
		
	print("Class: Particle count")
	for i in range(kmeans.n_clusters):
		print("{}: {}".format(lb,np.sum(lb==i)))

	
	## subtomogram average particles from each class
	#avgr=[Averagers.get("mean.tomo") for i in range(kmeans.n_clusters)]
	#for i in range(nptcl):
		#avgr[lb[i]].add_image(imgs00[i])

	#classfile = "{}/classes3d.hdf".format(options.outpath)
	#for i in range(kmeans.n_clusters):
		#e=avgr[i].finish()
		#e.process_inplace("normalize")
		#e.write_image(classfile,i)

	inlst=LSXFile(inptcls, True)
	outlsts=[]
	for lbl in sorted(np.unique(lb)):
		#outlst = LSXFile(inptcls.replace(".lst","_{}.lst".format(lbl)))
		outlst = LSXFile("{}/ptcls_cls{:02d}.lst".format(options.outpath, lbl+1))
		for i in range(nptcl):
			if lb[i]==lbl:
				l=inlst.read(i)
				outlst.write(-1,l[0], l[1],l[2])
		outlst.close()
	inlst.close()
	
	js0=js_open_dict(parmsfile)

	#pname=eval(js.keys()[0])[0]
	dics=[{} for i in range(kmeans.n_clusters)]
	for i in range(nptcl):
		if lb[i]>=0:
			k=keys[i]
			dics[lb[i]][k]=js0[k]
			
	js0.close()

	for i,d in enumerate(dics):
		js=js_open_dict("{}/particle_parms_{:02d}.json".format(options.outpath, i+1))
		js.update(d)
		js.close()
		os.system("e2spt_average.py --path {} --iter {} --threads 10 --sym {} --skippostp".format(options.outpath, i+1, options.sym))

	E2end(logid)


if __name__ == "__main__":
	main()

