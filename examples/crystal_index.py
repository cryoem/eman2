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
#from matplotlib.patches import Rectangle

from scipy import ndimage as ndi
from skimage.feature import peak_local_max
#from skimage import data, img_as_float
from scipy.spatial import cKDTree as KDTree
from scipy import interpolate
import scipy
import scipy.optimize as optimize

import os
import time
import itertools
import threading
import Queue

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <crystal image>

	Orient crystals imaged via electron microscopy.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--apix", required=True,type=float,default=False, help="Specify the Apix of your input images")
	parser.add_argument("--params", required=True,type=str,help="Lattice parameters separated by commas, i.e. '70.3,70.3,32.0,90.0,90.0,90.0'",default="", guitype='intbox', row=11, col=0, rowspan=1, colspan=1, mode="align")
	parser.add_argument("--slab", type=float, help="Specify the thickness of the central slab. Default is 10.0",default=10.0)
	parser.add_argument("--mindeg", type=float, help="Specify the minimum angle for initial orientation search. Default is -180.0",default=-180.0)
	parser.add_argument("--maxdeg", type=float, help="Specify the maximum angle for initial orientation search. Default is 180.0",default=180.0)
	parser.add_argument("--diameter", type=float, help="Specify the minimum spot diameter. Default is 5.0",default=5.0)
	parser.add_argument("--maxshift", type=float, help="Specify the maximum pixel shift when optimizing peak locations. Default is 32.0 pixels.",default=32.0)
	parser.add_argument("--exper_weight", type=float, help="Weight of penalty for spots in experimental data not found in the reference lattice. Default is 10.0.",default=10.0)
	parser.add_argument("--data_weight", type=float, help="Weight of penalty for points in reference lattice not found in the experimental data. Default is 1.0.",default=1.0)
	parser.add_argument("--plot", default=False, help="Show plot of reciprocal reference lattice points overlayed on input image and detected reflections.",action="store_true")
	parser.add_argument("--threads",type=int,help="Number of cores over which parallel optimization will be performed. Default is to use 1 core.",default=1)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	apix = float(options.apix)
	close = options.slab
	
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
		refined_coords = refine_peaks(coords,img,options.maxshift,options.maxshift*2.)

		ds = 1/(apix*nx) # angstroms per fourier pixel
		exper_max_radius = np.linalg.norm(refined_coords-nx/2,axis=1).max()
		print("Highest resolution reflection: {:.2f}A ({} pix)\n".format(1/(ds*exper_max_radius),int(round(exper_max_radius,0))))

		print("DETERMINING ORIENTATION")
		
		resolutions = np.asarray([1000.,100.,50.,25.,20.,18.,16.,10.,8.,5.,4.,3.,2.9,2.8,2.7,2.6,2.5,2.4,2.3,2.2,2.1,2.0,1.9,1.8])
		radii = [float(r) for r in 1/(resolutions*ds) if r <= nx/2.]

		# initial search range
		a_mindeg = options.mindeg
		a_maxdeg = options.maxdeg
		b_mindeg = options.mindeg
		b_maxdeg = options.maxdeg
		g_mindeg = options.mindeg
		g_maxdeg = options.maxdeg

		old_nrefs = 0.0

		count = 0
		for r in radii:
			ang = np.arccos(np.sqrt((r**2-options.diameter)/r**2))*180./np.pi
			#max_radius = r#1/np.sin(ang)**2
			#if r > nx/2: continue
			max_resol = 1/(r*ds)

			hkl_exper = np.vstack([refined_coords[:,1],refined_coords[:,0]]).T
			hkl_exper = (hkl_exper - np.mean(hkl_exper,axis=0))
			all_exper_radii = np.linalg.norm(hkl_exper,axis=1)
			hkl_exper = hkl_exper[all_exper_radii <= r]
			exper_radii = np.linalg.norm(hkl_exper,axis=1)
			hkl_exper = np.c_[hkl_exper,exper_radii]

			if len(hkl_exper) < 2: continue

			if len(hkl_exper) == old_nrefs:
				continue
			else: 
				old_nrefs = len(hkl_exper)

			print("\nAngular step: {:.2f} degrees\tRadius: {:.2f} pixels ({:.2f} Angstroms)\t{} Reflections".format(ang,r,max_resol,len(hkl_exper)))

			hkl_ref = generate_lattice(nx,apix,r,a,b,c,alpha,beta,gamma)

			dist_exper,idx_exper = scipy.spatial.cKDTree(hkl_exper[:,:2]).query([0,0],k=9)
			min_distance = np.min([dist_exper[1],dist_exper[3],dist_exper[5],dist_exper[7]])
			
			if count == 0:
				rngs = list(itertools.product(
					np.arange(a_mindeg,a_maxdeg+ang,ang),
					np.arange(b_mindeg,b_maxdeg+ang,ang),
					np.arange(g_mindeg,g_maxdeg+ang,ang)))
			else:
				rngs = []
				for s in solns:
					print(s)
					a_mindeg = s[0]-ang
					a_maxdeg = s[0]+ang
					b_mindeg = s[1]-ang
					b_maxdeg = s[1]+ang
					c_mindeg = s[2]-ang
					c_maxdeg = s[2]+ang
					rngs.extend(list(itertools.product(np.arange(a_mindeg,a_maxdeg+ang,ang),np.arange(b_mindeg,b_maxdeg+ang,ang),np.arange(g_mindeg,g_maxdeg+ang,ang))))
			count += 1

			start = time.time()
			resq=Queue.Queue(0)

			res=[0]*len(rngs)
			thds = []

			if options.verbose: sys.stdout.write("\rCreating {} threads".format(len(rngs)))
			
			for i in range(len(rngs)):
				if options.verbose and i % 100 == 0:
					sys.stdout.write("\rCreating {}/{} threads".format(i+1,len(rngs)))

				thd = threading.Thread(target=cost_async,args=(rngs[i],hkl_exper,hkl_ref,close,min_distance/10.,options.exper_weight,options.data_weight,i,resq))
				thds.append(thd)
			t0=time.time()

			if options.verbose: sys.stdout.flush()

			minval = np.inf

			thrtolaunch=0
			while thrtolaunch<len(thds) or threading.active_count()>1:
				if thrtolaunch<len(thds):
					while (threading.active_count()==options.threads ) : time.sleep(.1)
					if options.verbose and (thrtolaunch % 100 == 0 or len(thds)-thrtolaunch < 5):
							sys.stdout.write("\rSearched {}/{} orientations".format(thrtolaunch+1,len(thds)))
					thds[thrtolaunch].start()
					thrtolaunch+=1
				else: time.sleep(0.5)

				while not resq.empty():
					i,cx=resq.get()
					if cx < minval: minval = cx
					res[i]=cx

			for th in thds: th.join()

			solns = [rngs[i] for i,v in enumerate(res) if v == minval]

			sys.stdout.write("\t\t\tFound {} solutions:".format(len(solns)))
			print(solns)

			sys.stdout.flush()

		print("\n\nREFINING PARAMETERS")

		best_cost = np.inf
		scost = best_cost
		best_orient = None

		hkl_exper = np.vstack([refined_coords[:,1],refined_coords[:,0]]).T
		hkl_exper = (hkl_exper - np.mean(hkl_exper,axis=0))
		all_exper_radii = np.linalg.norm(hkl_exper,axis=1)
		hkl_exper = hkl_exper[all_exper_radii <= exper_max_radius]
		exper_radii = np.linalg.norm(hkl_exper,axis=1)
		hkl_exper = np.c_[hkl_exper,exper_radii]

		hkl_ref = generate_lattice(nx,apix,exper_max_radius,a,b,c,alpha,beta,gamma)

		for ii,soln in enumerate(solns):
			hkl_ref = generate_lattice(nx,apix,exper_max_radius,a,b,c,alpha,beta,gamma)

			refine1 = optimize.fmin(cost,soln,args=(hkl_exper,hkl_ref,close,min_distance/10.0,options.data_weight,options.exper_weight,),disp=False)
			az,alt,phi = refine1[0],refine1[1],refine1[2]

			refine_apix = optimize.fmin(apix_cost,[apix],args=(az,alt,phi,hkl_exper,hkl_ref,close,16.,nx,exper_max_radius,options.data_weight,options.exper_weight,a,b,c,alpha,beta,gamma,),disp=False)
			refine_apix = float(refine_apix[0])  #float(refine_apix.x)
			hkl_ref = generate_lattice(nx,refine_apix,exper_max_radius,a,b,c,alpha,beta,gamma)
			ds = 1/(refine_apix*nx) # angstroms per fourier pixel

			refine_close = optimize.fmin(close_cost,[close],args=(az,alt,phi,hkl_exper,hkl_ref,8.,5.0,options.data_weight,options.exper_weight,),disp=False)
			refine_close = refine_close[0]
			if refine_close <= 1.0: refine_close = close
			
			refine2 = optimize.fmin(cost,refine1,args=(hkl_exper,hkl_ref,refine_close,min_distance/10.0,options.data_weight,options.exper_weight,),disp=False) # re-refine orientation

			scost =  cost(refine2,hkl_exper,hkl_ref,refine_close,min_distance,options.data_weight,options.exper_weight) 
			if scost < best_cost:
				best_cost = scost
				best_orient = refine2
				best_refine_apix = refine_apix
				best_refine_close = refine_close

		sys.stdout.flush()
		ds = 1/(best_refine_apix*nx) # angstroms per fourier pixel

		print("Refined Apix: {:.2f} -> {:.2f}".format(apix,best_refine_apix))
		print("Refined thickness: {:.2f} -> {:.2f}".format(close,best_refine_close))
		print("Refined orientation: ({:.2f},{:.2f},{:.2f})\n".format(*best_orient))

		hkl_ref = generate_lattice(nx,refine_apix,exper_max_radius,a,b,c,alpha,beta,gamma)
		pln = get_plane(best_orient,hkl_ref,close=best_refine_close)
		pln = pln[np.argsort(pln[:,3])] # sort by radius

		if options.plot:
			plt.imshow(nimg,origin="lower",cmap=plt.cm.Greys_r)
			plt.scatter(hkl_exper[:,1]+nx/2,hkl_exper[:,0]+nx/2,c='b',marker='x')
			plt.scatter(pln[:,1]+nx/2,pln[:,0]+nx/2,c='r',marker='x')
			plt.axis("off")
			plt.title("(Az, Alt, Phi) -> ({:.2f},{:.2f},{:.2f})".format(*best_orient))
			plt.show()

		print("     xc        yc       zc          r     resol      h       k       l       raw_F       raw_p")
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
					print("{:8.1f},{:8.1f},{:8.1f}    {:6.1f}       inf    {:4d}    {:4d}    {:4d}     {:.2f}       {:.2f}".format(xc,yc,zc,r,int(h),int(k),int(l),raw_F,raw_p))
				sf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(h,k,l,raw_F,raw_p,r,resol))

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

# compare projected reference against experimental data
def compare(exper,data,min_dist,w1,w2):
	diff = len(data)-len(exper)
	#if len(data) < len(exper)*0.9:#0.9: 	# better to cover all black dots than to avoid whitespace
	#	return np.inf #,np.inf,1.0,diff,1.0,0.0,None,None
	#if len(data) > 3.*len(exper):
	#	return np.inf #,np.inf,1.0,diff,1.0,0.0,None,None
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

def wrap_cost(args):
	return cost(*args)

def cost(params,exper,ref,close,min_distance,w1,w2):
	#sys.stdout.write("\r{:2f},{:2f},{:2f}".format(*params))
	plane = get_plane(params,ref,close=close) # rotate reference
	return compare(exper,plane,min_distance,w1,w2)

def cost_async(params,exper,ref,close,min_distance,w1,w2,i,out):
	plane = get_plane(params,ref,close=close) # rotate reference
	c = compare(exper,plane,min_distance,w1,w2)
	out.put((i,c))

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
	#sys.stdout.write("\rApix: {}".format(params[0]))
	if params[0] > 1.5:
		return np.inf
	elif params[0] < 0.5:
		return np.inf
	new = generate_lattice(nx,params[0],max_radius,a,b,c,alpha,beta,gamma)
	return cost([az,alt,phi],exper,new,close,minimum_dist,w1,w2)

def close_cost(params,az,alt,phi,exper,ref,minimum_dist,weight=10.,w1=10.0,w2=1.0):
	if params[0]<1.0:
		return np.inf
	#sys.stdout.write("\rSlab: {}".format(params[0]))
	c = cost([az,alt,phi],exper,ref,params[0],minimum_dist,w1,w2)
	return params[0]*weight+c


if __name__ == "__main__":
	main()

# t = time.time()
# mypool = Pool(options.threads)
# res = [mypool.apply_async(cost,(r,hkl_exper,hkl_ref,close,min_distance,options.exper_weight,options.data_weight,)) for r in rngs]
# mypool.close()
# mypool.join()
# dt = time.time()-t
