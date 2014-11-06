#!/usr/bin/env python

# Author: James Michael Bell, 09/2014 (jmbell@bcm.edu), modified by Jesus Galaz-Montoya (jgmontoy@bcm.edu)
# Last modified 03/Nov/2014
# Copyright (c) 2014 Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 2111-1307 USA

from EMAN2 import *
import os
import sys
import time
import copy
import numpy as np

from EMAN2jsondb import JSTask,jsonclasses


def get_usage():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options]
	e2tvrecon3d.py reconstructs a 3D tomogram from a tilt series and tlt tilt 
	angles file. In order to reconstruct the original image, we minimize a 
	function that is the sum of a L2 data fit term and the total variation of 
	the image. Proximal iterations using the Fast Iterative Shrinking & 
	Thresholding Algorithm (FISTA) are used.
	"""
	return usage


def print_usage():
	usage = get_usage()
	print "usage " + usage;
	print "Please run '" + progname + " -h' for detailed options"
	
	return

def main():
	
	parser = EMArgumentParser(usage=get_usage())
	
	parser.add_argument("--tiltseries", default=None, help="""The input projections. 
		Project should usually have the xform.projection header attribute, which is 
		used for slice insertion""")
	parser.add_argument("--tltfile", type=str, default=None, help="""An IMOD-like .tlt file containing 
		alignment angles. If specified slices will be inserted using these angles in the 
		IMOD convention""")	
	parser.add_argument("--output", default="threed.hdf", help="""Output reconstructed 
		tomogram file name.""")
	parser.add_argument("--path",type=str,default='tvrecon_3d',help="""Directory in which 
		results will be stored.""")
	parser.add_argument("--iter", default=10, type=int, help="""Specify the number of 
		iterative reconstructions to complete before returning the final reconstructed volume. 
			The default number is 50.""")
	parser.add_argument("--beta", default=1.0, type=float, help="""Specify the total-variation 
		penalization/regularization weight parameter 'beta'. The default is 5.0.""")
	parser.add_argument("--subpix", default=1, type=int, help="""Specify the number of linear 
		subdivisions used to compute the projection of one image pixel onto a detector pixel.""")
	parser.add_argument("--savesinograms", action="store_true", default=False,help="""If provided,
		this option will save the sinogram for each 2-D slice (along Y) in the reconstruction 
		to disk.""")
	
	parser.add_argument("--inmemory",action='store_true',default=False,help="""If provided,
		this option will keep certain files open in memory instead of writing them and
		reading from disk every time. While this can be faster, it is very memory-intensive.""")
		
	parser.add_argument("--saveslices", action="store_true", default=False,help="""If provided,
		this option will save each reconstructed 2-D slice (along Y) to disk.""")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="""
		verbose level [0-9], higner number means higher level of verboseness.""")
	parser.add_argument("--parallel",type=str,default='thread:1',help="""Default=thread:1. 
		See http://blake.bcm.edu/emanwiki/EMAN2/Parallel""")
	parser.add_argument("--ppid", type=int, help="""Set the PID of the parent process, 
		used for cross platform PPID.""", default=-1)
	
	(options, args) = parser.parse_args()
	
	#Check that the minimum data required are available and sane, otherwise exit
	if not options.tiltseries:
		print "\nERROR: You must specficy --tiltseries"
		sys.exit(1)
	if not options.tltfile:
		print "\nERROR: You must specficy --tlt"
		sys.exit(1)
	if options.beta < 0.0:
		print "\nERROR: Parameter beta must be a positive, real number."
		sys.exit(1)
		
	#Parse and count tilt angles	
	tiltangles = np.asarray([ float( i ) for i in file( options.tltfile , "r" ) ])
	tiltangles = tiltangles.tolist()	
	
	nimgs = EMUtil.get_image_count( options.tiltseries )
	nangles = len( tiltangles )
	if nimgs != nangles:
		print """\nERROR: The number of images in the tiltseries, %d, does not match
			the number of angles in the tlt file, %d""" % ( nimgs, nangles )
		sys.exit(1)
	
	#Read essential info from image header
	hdr = EMData( options.tiltseries, 0 , True)
	apix = hdr['apix_x']
	xsize = hdr['nx']
	ysize = hdr['ny']
	
	#Once all parameters and data have passed wholesomeness checks, initialize logging
	logger = E2init(sys.argv,options.ppid)
	
	#Create new output directory for this run of the program
	options = makepath( options, options.path)
	
	if options.verbose > 2: 
		print "\nGenerating this new directory to save results to:", options.path
	
	options.path = os.getcwd() + "/" + options.path
	
	#Generate one projection operator for all 2D slice reconstructions
	if options.verbose:
		print "\nBuilding projection operator..."
	projection_operator = build_projection_operator( options, tiltangles, xsize, nimgs, None, 0, None )
	
	#Initialize parallelism
	if options.verbose:
		print "\n\n(e2tvrecon.py) INITIALIZING PARALLELISM\n\n"
		
	from EMAN2PAR import EMTaskCustomer
	etc=EMTaskCustomer(options.parallel)

	tasks=[]
	nimgs = len( tiltangles )
	for y in range( ysize ):	
		task=TVReconTask( options, xsize, ysize, y, projection_operator, tiltangles, nimgs )
		tasks.append( task )

	tids=etc.send_tasks(tasks)

	results = get_results(etc,tids,options)
	if options.verbose:
		print "\nThese many results %d were computed because there were these many tasks %d" % ( len(results), len(tasks) )
	
	results.sort()
	np_recons = []
	for i in range(len(results)):
		recon = results[i][-1]
		
		# Store 2D reconstructions in options.path if requested
		if options.saveslices:
			twodpath = options.path + "/slices.hdf"
			from_numpy( recon ).write_image( twodpath, i )
		
		np_recons.append( recon )
		
	reconstack = np.dstack( np_recons )
	threed_recon = from_numpy( reconstack )
	threed_recon['apix_x'] = apix
	threed_recon['apix_y'] = apix
	threed_recon['apix_z'] = apix
	
	threed_recon.rotate( 0, -90, -90 )
	threed_recon.write_image( options.path + '/' + options.output, 0  )
	
	
	E2end(logger)
	return


class TVReconTask(JSTask):
	"""This is a task object for the parallelism system. 
	It is responsible for reconstructing each 2D slice."""	
	
	def __init__(self,options,xsize,ysize,y,projection_operator,tiltangles,nimgs):
		JSTask.__init__(self,"TVRecon",'',{},"")
		self.classoptions={"options":options,"xsize":xsize,"ysize":ysize,"y":y,"projection_operator":projection_operator, "tiltangles":tiltangles,"nimgs":nimgs}
	
	def execute(self,callback=None):
		classoptions = self.classoptions
		#options = self.classoptions['options']
		#nimgs = self.classoptions['nimgs']		
		
		#Generate sinogram	
		sinogram = make_sinogram( classoptions['options'], classoptions['y'], classoptions['xsize'], classoptions['ysize'], classoptions['nimgs'] )			
	
		#Reconstruct sinogram into 2D image
		recon = twod_recon( classoptions['options'], sinogram, classoptions['y'], classoptions['projection_operator'], classoptions['tiltangles'], classoptions['ysize'] )

		return [ classoptions['y'], recon ]
	
	
jsonclasses["TVReconTask"]=TVReconTask.from_jsondict


def get_results(etc,tids,options):
	"""This will get results for a list of submitted tasks. Won't return until it has all requested results.
	aside from the use of options["ptcl"] this is fairly generalizable code. """
	
	# wait for them to finish and get the results
	# results for each will just be a list of (qual,Transform) pairs
	results=[0]*len(tids)		# storage for results
	ncomplete=0
	tidsleft=tids[:]
	while 1:
		time.sleep(5)
		proglist=etc.check_task(tidsleft)
		nwait=0
		for i,prog in enumerate(proglist):
			if prog==-1 : nwait+=1
			if prog==100 :
				r=etc.get_results(tidsleft[i])		#Results for a completed task
				
				if r:
					#print "r is", r
					y=r[0].classoptions["y"]		#Get the slice number from the task rather than trying to work back to it
					results[y] = r[1]				
				ncomplete+=1
		
		tidsleft=[j for i,j in enumerate(tidsleft) if proglist[i]!=100]		# remove any completed tasks from the list we ask about
		if options.verbose:
			print "  %d tasks, %d complete, %d waiting to start        \r"%(len(tids),ncomplete,nwait)
			sys.stdout.flush()
	
		if len(tidsleft)==0: break
		
	return results


"""
MAKE_SINOGRAM
Generates one 2D sinogram for a single specified pixel along the y axis of a tiltseries
and writes it to disk.
"""
def make_sinogram( options, y, xlen, ylen, num_imgs ):
	
	r = Region( 0, y, xlen, 1 )
	sinogram = []
	sinogramname = options.path + "/sinogram_" + str(y).zfill( len( str( ylen ))) + ".hdf"
	for imgnum in range( num_imgs ):
		prj = EMData() 
		prj.read_image( options.tiltseries, imgnum, False, r )
		
		if options.savesinograms:
			prj.write_image( sinogramname , imgnum)
		if options.inmemory:
			sinogram.append( prj )
		
	if options.verbose > 2:
		print "\n(e2tvrecon3d.py)(make_sinogram) Generated sinogram %i of %i" %( y+1, ylen )
	if options.inmemory:
		return sinogram
	else:
		return sinogramname
	

def twod_recon( options, sinogram, y, projection_operator, tiltangles, ylen ):

	npstack = []
	
	if options.inmemory:
		for img in sinogram:
			np_img = img.numpy().copy()
			npstack.append( np_img )
	else:
		n = EMUtil.get_image_count( sinogram )
		for i in range(n):
			img = EMData( sinogram, i )
			np_img = img.numpy().copy()
			npstack.append( np_img )
			
	data = np.asarray( npstack ).astype( np.float32 )
	xlen = img["nx"]
	
	projections = data.ravel()[:, np.newaxis]
	
	if options.savesinograms:
		for i in range(len(tiltangles)):
			from_numpy(projections[i*xlen:(i+1)*xlen]).write_image(options.path + '/projections' + str(y).zfill( len( str( ylen ))) + '.hdf',i)
	
	if options.verbose > 2: 
		print "\nStarting reconstruction for slice", y
		
	t1 = time.time()
	recon, energies = fista_tv( options, tiltangles, projections, projection_operator, None )
	t2 = time.time()
	
	if options.verbose > 3: 
		print "Reconstruction completed in %s s"%(str(t2-t1))
	
	return recon[-1]


def build_projection_operator( options, angles, l_x, n_dir=None, l_det=None, offset=0, pixels_mask=None ):
	try:
		from scipy import sparse
	except:
		print "\nERROR: SciPy not found. Must be installed to run e2tvrecon.py"
	
	"""
	Compute the tomography design matrix.
		
	Parameters
	----------
	angles : array of floats
		angles at which projections will be taken
	
	l_x : int
		linear size of image array
	
	n_dir : int, default l_x
		number of angles at which projections are acquired. n_dir projection 
		angles are regularly spaced between 0 and 180.
		
	l_det : int, default is l_x
		number of pixels in the detector. If l_det is not specified, we suppose 
		that l_det = l_x.
	
	subpix : int, default 1
		number of linear subdivisions used to compute the projection of one 
		image pixel onto a detector pixel.
	
	offset : int, default 0
		width of the strip of image pixels not covered by the detector.
	
	pixels_mask : 1-d ndarray of size l_x**2
		mask of pixels to keep in the matrix (useful for to removeing pixels 
		inside or outside of a circle, for example)
		
	Returns
	-------
	p : sparse matrix of shape (n_dir l_x, l_x**2), in csr format
		Tomography design matrix. The csr (compressed sparse row) allows for 
		efficient subsequent matrix multiplication. The dtype of the elements is 
		float32, in order to save memory.	
	"""
	if l_det is None:
		l_det = l_x
	X, Y = _generate_center_coordinates( options.subpix * l_x)
	X *= 1.0/options.subpix
	Y *= 1.0/options.subpix
	Xbig, Ybig = _generate_center_coordinates(l_det)
	Xbig *= (l_x - 2*offset) / float(l_det)
	orig = Xbig.min()
	labels = None
	if options.subpix > 1:
		# Block-group subpixels
		Xlab, Ylab = np.mgrid[:options.subpix * l_x, :options.subpix * l_x]
		labels = (l_x * (Xlab / options.subpix) + Ylab / options.subpix).ravel()
	if n_dir is None:
		n_dir = l_x
	weights, data_inds, detector_inds = [], [], []
	# Indices for data pixels. For each data, one data pixel
	# will contribute to the value of two detector pixels.
	for i, angle in enumerate(angles):
		# rotate data pixels centers
		Xrot = np.cos(angle*np.pi/180.) * X - np.sin(angle*np.pi/180.) * Y
		# compute linear interpolation weights
		inds, dat_inds, w = _weights_fast(Xrot, dx=(l_x - 2*offset)/float(l_det), orig=orig, labels=labels)
		# crop projections outside the detector
		mask = np.logical_and(inds >= 0, inds < l_det)
		weights.append(w[mask])
		detector_inds.append((inds[mask] + i * l_det).astype(np.int32))
		data_inds.append(dat_inds[mask])
	weights = np.concatenate(weights)
	weights /= options.subpix**2
	detector_inds = np.concatenate(detector_inds)
	data_inds = np.concatenate(data_inds)
	if pixels_mask is not None:
		if pixels_mask.ndim > 1:
			pixels_mask = pixels_mask.ravel()
		mask = pixels_mask[data_inds]
		data_inds = data_inds[mask]
		data_inds = rank_order(data_inds)[0]
		detector_inds = detector_inds[mask]
		weights = weights[mask]
	proj_operator = sparse.coo_matrix((weights, (detector_inds, data_inds)))
	return sparse.csr_matrix(proj_operator)


def _generate_center_coordinates(l_x):
	"""
	Compute the coordinates of pixels centers for an image of
	linear size l_x
	"""
	l_x = float(l_x)
	X, Y = np.mgrid[:l_x, :l_x]
	center = l_x / 2.
	X += 0.5 - center
	Y += 0.5 - center
	return X, Y


def _weights_fast(x, dx=1, orig=0, labels=None):
	"""
	Compute linear interpolation weights for projection array `x`
	and regularly spaced detector pixels separated by `dx` and
	starting at `orig`.
	"""
	x = np.ravel(x)
	floor_x = np.floor((x - orig) / dx).astype(np.int32)
	alpha = ((x - orig - floor_x * dx) / dx).astype(np.float32)
	inds = np.hstack((floor_x, floor_x + 1))
	weights = np.hstack((1 - alpha, alpha))
	data_inds = np.arange(x.size, dtype=np.int32)
	data_inds = np.hstack((data_inds, data_inds))
	if labels is not None:
		data_inds = np.hstack((labels, labels))
		order = np.argsort(inds)
		inds, data_inds, weights = inds[order], data_inds[order], weights[order]
		steps = np.nonzero(np.diff(inds) > 0)[0] + 1
		steps = np.concatenate(([0], steps))
		inds_s, data_inds_s, weights_s = [], [], []
		for i in range(len(steps) - 1):
			d, w = data_inds[steps[i]:steps[i+1]], weights[steps[i]:steps[i+1]]
			count = np.bincount(d, weights=w)
			mask = count>0
			w = count[mask]
			weights_s.append(w)
			datind = np.arange(len(mask))[mask] 
			data_inds_s.append(datind)
			detind = inds[steps[i]]*np.ones(mask.sum()) 
			inds_s.append(detind)
		inds = np.concatenate(inds_s)
		data_inds = np.concatenate(data_inds_s)
		weights = np.concatenate(weights_s)
	return inds, data_inds, weights


# FISTA ALGORITHM
def fista_tv( options, angles, y, H, mask=None):
	try:
		from scipy import sparse
	except:
		print "\nERROR: SciPy not found. Must be installed to run e2tvrecon.py"
	"""
	TV regression using FISTA algorithm
	(Fast Iterative Shrinkage/Thresholding Algorithm)
	
	Parameters
	----------
	
	y : ndarray of floats
		Measures (tomography projection). If H is given, y is a column
		vector. If H is not given, y is a 2-D array where each line
		is a projection along a different angle
	
	beta : float
		weight of TV norm
	
	niter : number of forward-backward iterations to perform
	
	H : sparse matrix
		tomography design matrix. Should be in csr format.
	
	mask : array of bools
	
	Returns
	-------
	
	res : list
		list of iterates of the reconstructed images
	
	energies : list
		values of the function to be minimized at the different
		iterations. Its values should be decreasing.
	
	Notes
	-----
	This algorithm minimizes iteratively the energy
	
	E(x) = 1/2 || H x - y ||^2 + beta TV bet(x) = f(x) + beta TV(x)
	
	by forward - backward iterations:
	
	u_n = prox_{gamma beta TV}(x_n - gamma nabla f(x_n)))
	t_{n+1} = 1/2 * (1 + sqrt(1 + 4 t_n^2))
	x_{n+1} = u_n + (t_n - 1)/t_{n+1} * (u_n - u_{n-1})
	
	References
	----------
	
	A. Beck and M. Teboulle (2009). A fast iterative
	shrinkage-thresholding algorithm for linear inverse problems.
	SIAM J. Imaging Sci., 2(1):183-202.
	
	Nelly Pustelnik's thesis (in French),
	http://tel.archives-ouvertes.fr/tel-00559126_v4/
	Paragraph 3.3.1-c p. 69 , FISTA
	"""
	n_meas, n_pix = H.shape
	if mask is not None:
		l = len(mask)
	else:
		l = int(np.sqrt(n_pix))
	n_angles = len(angles)
	Ht = sparse.csr_matrix(H.transpose())
	x0 = np.zeros(n_pix)[:, np.newaxis]
	res, energies = [], []
	gamma = .9/ (l * n_angles)
	x = x0
	u_old = np.zeros((l, l))
	t_old = 1
	for i in range( options.iter ):
		if options.verbose > 9:
			print "\nFistaTV iteration", i
		eps = 1.e-4
		err = H * x - y
		back_proj = Ht * err
		tmp = x - gamma * back_proj
		if mask is not None:
			tmp2d = np.zeros((l, l))
			tmp2d[mask] = tmp.ravel()
		else:
			tmp2d = tmp.reshape((l, l))
		u_n = tv_denoise_fista(tmp2d, weight=options.beta*gamma, eps=eps)
		t_new = (1 + np.sqrt(1 + 4 * t_old**2))/2.
		t_old = t_new
		x = u_n + (t_old - 1)/t_new * (u_n - u_old)
		u_old = u_n
		res.append(x)
		data_fidelity_err = 1./2 * (err**2).sum()
		tv_value = options.beta * tv_norm(x)
		# tv_norm_anisotropic(x)
		# tv_l0_norm(x)
		energy = data_fidelity_err + tv_value
		energies.append(energy)
		if mask is not None:
			x = x[mask][:, np.newaxis]
		else:
			x = x.ravel()[:, np.newaxis]
	return res, energies


def tv_denoise_fista(im, weight=50, eps=5.e-5, n_iter_max=200, check_gap_frequency=3):
	"""
	Perform total-variation denoising on 2-d and 3-d images
	
	Find the argmin `res` of
		1/2 * ||im - res||^2 + weight * TV(res),
	
	where TV is the isotropic l1 norm of the gradient.
	
	Parameters
	----------
	im: ndarray of floats (2-d or 3-d)
		input data to be denoised. `im` can be of any numeric type,
		but it is cast into an ndarray of floats for the computation
		of the denoised image.
	
	weight: float, optional
		denoising weight. The greater ``weight``, the more denoising (at
		the expense of fidelity to ``input``)
	
	eps: float, optional
		precision required. The distance to the exact solution is computed
		by the dual gap of the optimization problem and rescaled by the l2
		norm of the image (for contrast invariance).
	
	n_iter_max: int, optional
		maximal number of iterations used for the optimization.
	
	Returns
	-------
	out: ndarray
		denoised array
	
	Notes
	-----
	The principle of total variation denoising is explained in
	http://en.wikipedia.org/wiki/Total_variation_denoising
	
	The principle of total variation denoising is to minimize the
	total variation of the image, which can be roughly described as
	the integral of the norm of the image gradient. Total variation
	denoising tends to produce "cartoon-like" images, that is,
	piecewise-constant images.
	
	This function implements the FISTA (Fast Iterative Shrinkage
	Thresholding Algorithm) algorithm of Beck et Teboulle, adapted to
	total variation denoising in "Fast gradient-based algorithms for
	constrained total variation image denoising and deblurring problems"
	(2009).
	"""
	if not im.dtype.kind == 'f':
		im = im.astype(np.float)
	shape = [im.ndim, ] + list(im.shape)
	grad_im = np.zeros(shape)
	grad_aux = np.zeros(shape)
	t = 1.
	i = 0
	while i < n_iter_max:
		error = weight * div(grad_aux) - im
		grad_tmp = gradient(error)
		grad_tmp *= 1./ (8 * weight)
		grad_aux += grad_tmp
		grad_tmp = _projector_on_dual(grad_aux)
		t_new = 1. / 2 * (1 + np.sqrt(1 + 4 * t**2))
		t_factor = (t - 1) / t_new
		grad_aux = (1 + t_factor) * grad_tmp - t_factor * grad_im
		grad_im = grad_tmp
		t = t_new
		if (i % check_gap_frequency) == 0:
			gap = weight * div(grad_im)
			new = im - gap
			dgap = dual_gap(im, new, gap, weight)
			if dgap < eps:
				break
		i += 1
	return new


# TV NORMS
def tv_norm(img):
	"""Compute the (isotropic) TV norm of an image"""
	grad_x1 = np.diff(img, axis=0)
	grad_x2 = np.diff(img, axis=1)
	return np.sqrt(grad_x1[:, :-1]**2 + grad_x2[:-1, :]**2).sum()

def tv_l0_norm( img ):
	"""Compute the (isotropic) TV norm of a 2D image"""
	grad_x1 = np.diff(img, axis=0)
	grad_x2 = np.diff(img, axis=1)
	return (grad_x1[:, :-1]**2 + grad_x2[:-1, :]**2 > 0).mean()

def tv_norm_anisotropic( img ):
	"""Compute the anisotropic TV norm of an image"""
	grad_x1 = np.diff(img, axis=0)
	grad_x2 = np.diff(img, axis=1)
	return np.abs(grad_x1[:, :-1]).sum() + np.abs(grad_x2[:-1, :]).sum()


def rank_order(image):
	"""Return an image of the same shape where each pixel is the index of the 
	pixel value in the ascending order of the unique values of `image`, aka the 
	rank-order value.
	
	Parameters
	----------
	image: ndarray
	
	Returns
	-------
	labels: ndarray of type np.uint32, of shape image.shape
		New array where each pixel has the rank-order value of the corresponding 
		pixel in `image`. Pixel values are between 0 and n - 1, where n is the 
		number of distinct unique values in ''image`.
	
	original_values: 1-d ndarray
		Unique original values of `image`
	
	Examples
	--------
	>>> a = np.array([[1, 4, 5], [4, 4, 1], [5, 1, 1]])
	>>> a
	array([[1, 4, 5],
		[4, 4, 1],
		[5, 1, 1]])
	>>> rank_order(a)
	(array([[0, 1, 2],
		[1, 1, 0],
		[2, 0, 0]], dtype=uint32), array([1, 4, 5]))
	>>> b = np.array([-1., 2.5, 3.1, 2.5])
	>>> rank_order(b)
	(array([0, 1, 2, 1], dtype=uint32), array([-1. ,  2.5,  3.1]))
	"""
	flat_image = image.ravel()
	sort_order = flat_image.argsort().astype(numpy.uint32)
	flat_image = flat_image[sort_order]
	sort_rank = numpy.zeros_like(sort_order)
	is_different = flat_image[:-1] != flat_image[1:]
	numpy.cumsum(is_different, out=sort_rank[1:])
	original_values = numpy.zeros((sort_rank[-1] + 1,), image.dtype)
	original_values[0] = flat_image[0]
	original_values[1:] = flat_image[1:][is_different]
	int_image = numpy.zeros_like(sort_order)
	int_image[sort_order] = sort_rank
	return (int_image.reshape(image.shape), original_values)


def div(grad):
	""" Compute divergence of image gradient """
	res = np.zeros(grad.shape[1:])
	for d in range(grad.shape[0]):
		this_grad = np.rollaxis(grad[d], d)
		this_res = np.rollaxis(res, d)
		this_res[:-1] += this_grad[:-1]
		this_res[1:-1] -= this_grad[:-2]
		this_res[-1] -= this_grad[-2]
	return res


def gradient(img):
	""" 
	Compute gradient of an image
	
	Parameters
	===========
	img: ndarray
		N-dimensional image
	
	Returns
	=======
	gradient: ndarray
		Gradient of the image: the i-th component along the first axis is the 
		gradient along the i-th axis of the original array img
	"""
	shape = [img.ndim, ] + list(img.shape)
	gradient = np.zeros(shape, dtype=img.dtype)
	# 'Clever' code to have a view of the gradient with dimension i stop at -1
	slice_all = [0, slice(None, -1),]
	for d in range(img.ndim):
		gradient[slice_all] = np.diff(img, axis=d)
		slice_all[0] = d + 1
		slice_all.insert(1, slice(None))
	return gradient


def _projector_on_dual(grad):
	"""
	Modifies in place the gradient to project iton the L2 unit ball
	"""
	norm = np.maximum(np.sqrt(np.sum(grad**2, 0)), 1.)
	for grad_comp in grad:
		grad_comp /= norm
	return grad


def dual_gap(im, new, gap, weight):
	"""
	dual gap of total variation denoising
	see "Total variation regularization for fMRI-based prediction of behavior", 
	by Michel et al. (2011) for a derivation of the dual gap
	"""
	im_norm = (im**2).sum()
	gx, gy = np.zeros_like(new), np.zeros_like(new)
	gx[:-1] = np.diff(new, axis=0)
	gy[:, :-1] = np.diff(new, axis=1)
	if im.ndim == 3:
		gz = np.zeros_like(new)
		gz[..., :-1] = np.diff(new, axis=2)
		tv_new = 2 * weight * np.sqrt(gx**2 + gy**2 + gz**2).sum()
	else:
		tv_new = 2 * weight * np.sqrt(gx**2 + gy**2).sum()
	dual_gap = (gap**2).sum() + tv_new - im_norm + (new**2).sum()
	
	try:
		
		ret = 0
		divisor = float( im_norm ) * float( dual_gap )
		if divisor:
			ret = 0.5 / divisor
			
		return ret
	except:
		return 0


def makepath(options, stem=''):
	if options.verbose > 5: print "makepath function called"
	if options.path and ("/" in options.path or "#" in options.path):
		print "Path specifier should be the name of a subdirectory to use in the current directory."
		print "Neither '/' or '#' can be included. Please edit your --path argument accordingly."
		sys.exit(1)
	if not options.path:
		options.path = stem + '_01'
		if options.verbose > 5:
			print "--path was not specified, therefore it will have the default value"
	files=os.listdir(os.getcwd())
	while options.path in files:
		if '_' not in options.path:
			options.path = options.path + '_00'
		else:
			#jobtag=''
			components=options.path.split('_')
			if components[-1].isdigit():
				components[-1] = str(int(components[-1])+1).zfill(2)
			else:
				components.append('00')
			options.path = '_'.join(components)
	if options.verbose > 5: print "The new options.path is", options.path
	if options.path not in files:
		if options.verbose > 5:
			print "Creating the following path: ", options.path
		os.system('mkdir ' + options.path)
	return options


if __name__=="__main__":
	main()

