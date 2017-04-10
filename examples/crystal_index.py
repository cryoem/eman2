#!/usr/bin/env python

#
# Author: Michael Bell, 02/12/2017 (jmbell@bcm.edu).
# Copyright (c) 2017-2020 Baylor College of Medicine
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


from EMAN2 import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage import data, img_as_float
from scipy.spatial import cKDTree as KDTree
from scipy import interpolate
import scipy
import scipy.optimize as optimize

import os
from multiprocessing.pool import ThreadPool, Pool
import itertools

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <crystal image>

	Orient crystals imaged via electron microscopy.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--apix", type=float,default=False, help="Specify the Apix of your input images")
	parser.add_argument("--params", type=str,help="Lattice parameters separated by commas, i.e. '70.3,70.3,32.0,90.0,90.0,90.0'",default="", guitype='intbox', row=11, col=0, rowspan=1, colspan=1, mode="align")
	parser.add_argument("--maxshift", type=float, help="Specify the maximum pixel shift when optimizing peak locations. Default is 5.0 pixels.",default=5.0)
	parser.add_argument("--boxsize", type=int, help="Specify the size of the box within which peaks will be refined. Default is 64 pixels.",default=64)
	parser.add_argument("--exper_weight", type=float, help="Weight of penalty for spots in experimental data not found in the reference lattice. Default is 10.0.",default=10.0)
	parser.add_argument("--data_weight", type=float, help="Weight of penalty for points in reference lattice not found in the experimental data. Default is 1.0.",default=1.0)
	parser.add_argument("--threads",type=int,help="Number of cores over which parallel optimization will be performed. Default is to use all cores.",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	apix = float(options.apix)
	
	try:
		params = [float(p) for p in options.params.split(",")]
		a,b,c,alpha,beta,gamma = params
	except:
		print("Could not read lattice parameters. Please specify as a comma separated list containing 'a,b,c,alpha,beta,gamma'.")
	
	for fn in args:

		try:
			print("READING {}".format(fn))
			orig = EMData(fn)
		except:
			print("Could not find {}".format(fn))
			sys.exit(1)

		# PREPROCESSING

		nx = orig["nx"]
		orig.process_inplace("normalize.edgemean")
		#orig.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
		#orig.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.34})
		nx = min(orig["nx"],orig["ny"]) # clip to min x,y to obtain square image
		reg = Region(0,0,nx,nx)
		orig = orig.get_clip(reg)
		#orig.process_inplace("filter.highpass.gauss",{"cutoff_freq":0.01}) # remove strong signal at origin for improved peak finding
		orig.process_inplace("filter.xyaxes0",{"neighbornorm":2})
		#orig.process_inplace("mask.gaussian",{"outer_radius":orig["nx"]/8}) # overly stringent apodization
		#orig.process_inplace("mask.decayedge2d",{"width":nx/4}) # simple apodization
		orig.process_inplace("mask.gaussian",{"outer_radius":orig["nx"]/8,"exponent":3.0})

		norig = orig.numpy().copy()
		fnorig = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(norig)))

		img = orig.process("math.realtofft") # obtain an amplitude image
		img.process_inplace("normalize")
		img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1}) # strong lowpass
		ra = img.process("math.rotationalaverage")
		img -= ra
		img.process_inplace("threshold.belowtozero",{"minval":img["mean"]+3.0*img["sigma"]})
		#img.process_inplace("threshold.notzero")
		nimg = img.numpy().copy()

		plt.imshow(nimg)
		plt.show()

		print("\nDETECTING SPOTS")

		# dilation of nimg with a size*size structuring element
		image_max = ndi.maximum_filter(nimg, size=np.min([a,b,c]).astype(int), mode='constant') # measured size of spot (rough) # 30
		# compare image_max and nimg to find the coordinates of local maxima
		tmp = peak_local_max(nimg, min_distance=np.max([a,b,c]).astype(int)) # measured minimum distance from origin to first point (rough) # 50
		coords = []
		for coord in tmp:
			if np.abs(np.linalg.norm(c-nimg.shape[0]/2.)) <= (nx/2.):
				coords.append(coord)
		coords = np.asarray(coords)
		refined_coords = refine_peaks(coords,img,options.maxshift,options.boxsize)

		ds = 1/(apix*nx) # angstroms per fourier pixel
		max_radius = np.linalg.norm(refined_coords-nx/2,axis=1).max()
		print("HIGHEST RESOLUTION: {:.2f}A ({} pix)\n".format(1/(ds*max_radius),int(round(max_radius,0))))

		print("DETERMINING ORIENTATION")
		#print("A: {}\nB: {}\nC: {}\nAlpha: {}\nBeta: {}\nGamma: {}".format(a,b,c,alpha,beta,gamma))

		hkl_exper = np.vstack([refined_coords[:,1],refined_coords[:,0]]).T
		hkl_exper = (hkl_exper - np.mean(hkl_exper,axis=0))
		exper_radii = np.linalg.norm(hkl_exper,axis=1)
		hkl_exper = np.c_[hkl_exper,exper_radii]

		hkl_ref = generate_lattice(nx,apix,max_radius,a,b,c,alpha,beta,gamma)

		#close = 10.0
		#rot = rotate(hkl_ref,[0.,0.,0.])
		#plane = rot[np.where(np.logical_and(rot[:,2] >= -1*close/2., rot[:,2] <= close/2.))]
		#ref_radii = plane[:,3] # radii of points in reference volume
		#within = plane[plane[:,3] > 0.0,3].min()
		#print(within)
		#within = 250.0 # furthest experimental datapoint from the origin
		#obs = np.logical_and(ref_radii<=within,ref_radii>0.0) # compare experimental spots to these quasi-coplanar points
		dist_exper,idx_exper = scipy.spatial.cKDTree(hkl_exper[:,:2]).query([0,0],k=9)
		#test_tree_b = scipy.spatial.cKDTree(plane[:,:2])
		#dist_ref,idx_ref = test_tree_b.query([0,0],k=9)
		min_distance = np.min([dist_exper[1],dist_exper[3],dist_exper[5],dist_exper[7]])
		#nmin = 0
		#nmax = 10
		#rngs2 = (slice(nmin,nmax+1,1),slice(nmin,nmax+1,1),slice(1,nmax+1,1))
		#lc = optimize.brute(find_first_reflection,rngs2,finish=None,disp=False,full_output=True,args=(nx,apix,a,b,c,1.0,))

		close = 10.0 # initial thickness of central slab within which we expect to observe spots

		a_mindeg = -180.
		a_maxdeg = 180.
		a_rng = a_maxdeg - a_mindeg
		astep = 0.5

		b_mindeg = -180.
		b_maxdeg = 180.
		b_rng = b_maxdeg - b_mindeg
		bstep = 0.5

		g_mindeg = -180.
		g_maxdeg = 180.
		g_rng = g_maxdeg - g_mindeg
		gstep = 1.0

		rngs = list(itertools.product(
			np.arange(a_mindeg,a_maxdeg+astep,astep),
			np.arange(b_mindeg,b_maxdeg+bstep,bstep),
			np.arange(g_mindeg,g_maxdeg+gstep,gstep)))

		t = time.time()
		mypool = Pool(options.threads)
		res = [mypool.apply_async(cost,(r,hkl_exper,hkl_ref,close,min_distance,options.exper_weight,options.data_weight,)) for r in rngs]
		mypool.close()
		mypool.join()
		dt = time.time()-t

		res1 = rngs[np.argmin(res)]
		refine1 = optimize.fmin(cost,res1,args=(hkl_exper,hkl_ref,close,min_distance,options.data_weight,options.exper_weight,),disp=False)
		az,alt,phi = refine1

		print("Az, Alt, Phi -> {:.2f},{:.2f},{:.2f}".format(az,alt,phi))
		#print("Az: {}\tAlt: {}\tPhi: {}\n".format(az,alt,phi))

		# pln = get_plane(refine1,hkl_ref,close)
		# pln = pln[np.argsort(pln[:,3])] # sort by radius
		# plt.imshow(nimg,cmap=plt.cm.Greys_r)
		# plt.scatter(pln[:,1]+nx/2,pln[:,0]+nx/2)
		# plt.show()

		#sys.exit(1)

		print("REFINING PARAMETERS") #xtol=0.01, ftol=0.01,maxfun=1e5
		
		refine_apix = optimize.fmin(apix_cost,[apix],args=(az,alt,phi,hkl_exper,hkl_ref,close,16.,nx,max_radius,options.data_weight,options.exper_weight,a,b,c,alpha,beta,gamma,))#,disp=False)
		refine_apix = refine_apix[0] #float(refine_apix.x)
		hkl_ref = generate_lattice(nx,refine_apix,max_radius,a,b,c,alpha,beta,gamma)
		ds = 1/(refine_apix*nx) # angstroms per fourier pixel
		print("APIX: {:.2f} -> {:.2f}".format(apix,refine_apix))

		# pln = get_plane(refine1,hkl_ref,close)
		# pln = pln[np.argsort(pln[:,3])] # sort by radius
		# plt.imshow(nimg,cmap=plt.cm.Greys_r)
		# plt.scatter(pln[:,0]+nx/2,pln[:,1]+nx/2)
		# plt.show()

		#sys.exit(1)

		refine_close = optimize.fmin(close_cost,[close],args=(az,alt,phi,hkl_exper,hkl_ref,8.,50.,options.data_weight,options.exper_weight,))#,disp=False)
		if refine_close[0] < 0.5: refine_close[0] = close
		print("SLAB THICKNESS: {:.2f} -> {:.2f}".format(close,refine_close[0]))

		# pln = get_plane(refine1,hkl_ref,refine_close[0])
		# pln = pln[np.argsort(pln[:,3])] # sort by radius
		# plt.imshow(nimg,cmap=plt.cm.Greys_r)
		# plt.scatter(pln[:,0]+nx/2,pln[:,1]+nx/2)
		# plt.show()		

		refine2 = optimize.fmin_cg(cost,refine1,args=(hkl_exper,hkl_ref,refine_close[0],16.,options.data_weight,options.exper_weight,))#,disp=False) # re-refine orientation
		print("ORIENTATION: {}\n".format(refine2))
		print("Az, Alt, Phi -> {:.2f},{:.2f},{:.2f}".format(az,alt,phi))
		
		# pln = get_plane(refine2,hkl_ref,refine_close[0])
		# pln = pln[np.argsort(pln[:,3])] # sort by radius
		# plt.imshow(nimg,cmap=plt.cm.Greys_r)
		# plt.scatter(pln[:,1]+nx/2,pln[:,0]+nx/2)
		# plt.show()

		print("     xc        yc       zc          r           resol    h      k     l      raw_F     raw_p")
		with open("{}.sf".format(fn.split(".")[0]),"w") as sf: # create a quasi.sf file for this image
			for nrow,row in enumerate(pln):
				xc,yc,zc,r,h,k,l = row[:7]
				if r > 0: resol = 1/(r*ds)
				else: resol = "inf"
				if nrow in range(9):
					raw_F,raw_p = get_sf(fnorig,int(xc)+nx/2,int(yc)+nx/2,2)#,show=True)
				else:
					raw_F,raw_p = get_sf(fnorig,int(xc)+nx/2,int(yc)+nx/2,2)#,show=False)
				try:
					print("{:8.1f},{:8.1f},{:8.1f}    {:6.1f}    {:6.1f}    {:4d}    {:4d}    {:4d}    {:8.2f}    {:8.2f}".format(xc,yc,zc,r,resol,int(h),int(k),int(l),raw_F,raw_p))
				except:
					print("{:8.1f},{:8.1f},{:8.1f}    {:6.1f}        inf   {:4d}     {:4d}    {:4d}    {:.2f}    {:.2f}".format(xc,yc,zc,r,int(h),int(k),int(l),raw_F,raw_p))
				sf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(h,k,l,raw_F,raw_p,r,resol))

def wrap_cost(args):
	return cost(*args)

def get_sf(fft,xc,yc,nn=2):
	reg = fft[xc-nn/2:xc+nn/2+1,yc-nn/2:yc+nn/2+1]
	amp = np.absolute(reg).max()
	phase = np.angle(np.sum(np.real(reg).ravel()) + 1j*np.sum(np.imag(reg).ravel()),deg=True)
	return amp,phase

def twoD_Gaussian((x, y), xo, yo, amplitude, sigma_x, sigma_y, theta, offset=0.0):
	xo = float(xo)
	yo = float(yo)
	a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
	b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
	c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
	g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
	return g.ravel()

def refine_peaks(coords,img,ms,bs):
	x = np.linspace(0,bs,bs)
	y = np.linspace(0,bs,bs)
	x,y = np.meshgrid(x,y)
	refined_coords = []
	for (yc,xc) in coords:
		r = Region(xc-bs/2,yc-bs/2,bs,bs)
		ereg = img.get_clip(r)
		nreg = ereg.numpy().copy()
		bds = [(bs/2-ms,bs/2+ms),(bs/2-ms,bs/2+ms),(0,np.max(nreg)*2),(0.1,8.0),(0.1,8.0),(0.,90.)] # lower, upper bounds
		initial_guess = [bs/2.,bs/2.,np.max(nreg),1.,1.,0.]
		try:
			popt,pcov=optimize.curve_fit(twoD_Gaussian,(x,y),nreg.ravel(),p0=initial_guess,bounds=bds)
			popt = [p for p in popt]
		except:
			popt = initial_guess
		popt[0] = popt[0] + xc - bs/2
		popt[1] = popt[1] + yc - bs/2
		refined_coords.append([popt[0],popt[1]])

	return np.asarray(refined_coords)

def generate_lattice(nx,apix,max_radius,a,b,c,alpha,beta,gamma):
	# define pixel spacing relative to apix and unit cell parameters
	const = nx*apix

	astar = const/a
	bstar = const/b
	cstar = const/c

	#create lattice coordinates
	astars = np.arange(0,const+astar,astar,dtype=np.float)
	bstars = np.arange(0,const+bstar,bstar,dtype=np.float)
	cstars = np.arange(0,const+cstar,cstar,dtype=np.float)

	astars -= astars[len(astars)/2] # center about origin
	bstars -= bstars[len(bstars)/2]
	cstars -= cstars[len(cstars)/2]

	h_inds = astars / astar
	k_inds = bstars / bstar
	l_inds = cstars / cstar

	hkl_ref = np.asarray(np.meshgrid(astars,bstars,cstars)).T
	hkl_inds = np.round(np.asarray(np.meshgrid(h_inds,k_inds,l_inds))).astype(int).T

	hkl_inds = np.swapaxes(hkl_inds, 0, 1)
	hkl_inds = np.swapaxes(hkl_inds, 1, 2)

	hkl_ref = np.swapaxes(hkl_ref, 0, 1)
	hkl_ref = np.swapaxes(hkl_ref, 1, 2)

	radii_ref = np.linalg.norm(hkl_ref,axis=-1) # radii of points in reference volume
	hkl_ref = np.concatenate([hkl_ref,np.expand_dims(radii_ref,3),hkl_inds],axis=3)

	# keep only points within a ball the size of measured points in our experimental data
	return hkl_ref[np.logical_and(radii_ref < max_radius, radii_ref > -max_radius)]

# def find_first_reflection(hkl,obj,nx,apix,a,b,c,thresh=1.0):
# 	ds = 1/(apix*nx)
# 	astar = 1/(a*ds) #50.551518297430448
# 	bstar = 1/(b*ds) #50.551518297430448
# 	cstar = 1/(c*ds) #104.66147232670608
# 	mag = np.linalg.norm([hkl[0]*astar,hkl[1]*bstar,hkl[2]*cstar])
# 	return np.abs(obj - mag)

# compare projected reference against experimental data
def compare(exper,data,min_dist,w1,w2):
	diff = len(data)-len(exper)
	if len(data) < len(exper)*0.9:#0.9: 	# better to cover all black dots than to avoid whitespace
		return np.inf #,np.inf,1.0,diff,1.0,0.0,None,None
	if len(data) > 3*len(exper):
		return np.inf #,np.inf,1.0,diff,1.0,0.0,None,None
	dist = []
	a = exper[:,:2]
	b = data[:,:2]
	cdist = scipy.spatial.distance.cdist(a,b)
	nopair_count = 0
	for i in range(np.max(cdist.shape)):
		m = np.min(cdist)
		i,j = np.where(cdist==m)
		if m == np.inf:
			nopair_count+=1
			continue
		else: dist.append(m)
		cdist[i,:] = np.inf
		cdist[:,j] = np.inf
	proximity = np.sum(dist)/len(dist)
	tree_exper = scipy.spatial.cKDTree(exper[:,:2])
	tree_data = scipy.spatial.cKDTree(data[:,:2])
	pairs = tree_exper.query_ball_tree(tree_data,r=min_dist)
	pairs2 = tree_data.query_ball_tree(tree_exper,r=min_dist)
	notpaired = len([p for p in pairs if len(p) == 0]) # number of points in ref without correspondance in the experimental data
	notpaired2 = len([p for p in pairs2 if len(p) == 0])
	energy = proximity + w1 * notpaired + w2 * notpaired2
	return energy#,proximity,w1,notpaired,w2,notpaired2,pairs,pairs2

def cost(params,exper,ref,close,min_distance,w1,w2):
	plane = get_plane(params,ref) # rotate reference
	return compare(exper,plane,min_distance,w1,w2)

def get_plane(params,ref,close=10.0): # rotate reference and return "slab" associated with orientation
	rot = rotate(ref,params) #[az,alt,phi]
	return rot[np.where(np.logical_and(rot[:,2] >= -1*close/2., rot[:,2] <= close/2.))]

def rotate(cloud,params):
	az,alt,phi=params
	t = Transform({"type":"eman","az":az,"alt":alt,"phi":phi})
	mat = np.asarray(t.get_matrix()).reshape(3,4)[:,:3]
	rot = np.dot(cloud[:,:3], mat.T) # only rotate coordinates, not their amp/phase
	return np.hstack([rot,cloud[:,3:]])

def apix_cost(params,az,alt,phi,exper,ref,close,minimum_dist,nx,max_radius,w1,w2,a,b,c,alpha,beta,gamma):
	new = generate_lattice(nx,params[0],max_radius,a,b,c,alpha,beta,gamma)
	return cost([az,alt,phi],exper,new,close,minimum_dist,w1,w2)

def close_cost(params,az,alt,phi,exper,ref,minimum_dist=8.,weight=50.,w1=10.0,w2=1.0):
	c = cost([az,alt,phi],exper,ref,params[0],minimum_dist,w1,w2)
	return params[0]*weight+c

if __name__ == "__main__":
	main()
