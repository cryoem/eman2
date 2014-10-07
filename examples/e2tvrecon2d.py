#!/usr/bin/env python

#
# Author: James Michael Bell, 2014 (jmbell@bcm.edu)
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
#

from EMAN2 import *
import os
import sys
import scipy
import time
import ntpath
import numpy as np
from scipy import sparse
from scipy import ndimage


def get_usage():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options]
	e2tvrecon.py reconstructs a 2D image from a incomplete set of projections. 
	In addition, gaussian noise may be added to the image before projections are
	made to more accurately simulate experimentally collected data. In order to 
	reconstruct the original image, we minimize a function that is the sum of a 
	L2 data fit and the total variation of the image. Proximal iterations using 
	the FISTA scheme are used. The original source of FISTA implemented can be 
	found at https://github.com/emmanuelle/tomo-tv.
	"""
	return usage


def print_usage():
	usage = get_usage()
	print "usage " + usage;
	print "Please run '" + progname + " -h' for detailed options"


def main():
	parser = EMArgumentParser(usage=get_usage())
	parser.add_argument("--tiltseries", default=None, help="The input projections. Project should usually have the xform.projection header attribute, which is used for slice insertion")
	parser.add_argument("--imgnum", default=None, type=int, help="The image number which will be read from the stack when reconstructing an image from a user specified tiltseries.")
	parser.add_argument("--testdata", default=None, help="A 2D image to project a number of times (specified by --nslices and --tiltrange) and then reconstructed via compressed sensing.")
	parser.add_argument("--tlt", default=None, type=str, help="An imod tlt file containing alignment angles. If specified slices will be inserted using these angles in the IMOD convention")
	parser.add_argument("--nslices", default=120, type=int, help="Specify the number slices into which an image will be projected. Only applicable when using the --testdata option.")
	parser.add_argument("--tiltrange", default='60.0', type=str, help="Specify the range of degrees over which data was collected. This defaults to 60 degrees, resulting in the generation of projections from -60.0 to 60.0 degrees.")
	parser.add_argument("--output", default="recon.hdf", help="Output reconstructed tomogram file name.")
	parser.add_argument("--noise",action="store_true",default=False, help="If true, noise will be added to the image before reconstruction.")
	parser.add_argument("--noisiness",default=0.1, type=float, help="Multiply noise by a specified factor. The default value is 0.1")
	parser.add_argument("--path",type=str,default='recon',help="Directory in which results will be stored.")
	parser.add_argument("--niters", default=100, type=int, help="Specify the number of iterative reconstructions to complete before returning the final reconstructed volume.")
	parser.add_argument("--beta", default=0.2, type=float, help="Specify the total-variation regularization/penalization weight parameter 'beta'. The default value is 0.2. Note that this parameter must be greater than 0, but values much greater than 1 will produce cartoon-like reconstructions as a result of the total variation denoising procedure.")
	parser.add_argument("--subpix", default=1, type=int, help="Specify the number of linear subdivisions used to compute the projection of one image pixel onto a detector pixel. Note that this parameter must be a positive integer.")
	parser.add_argument("--fsc",action="store_true",default=False, help="Generate a fourier shell correlation plot comparing the input and output data.")
	parser.add_argument("--norm",default=None, type=str, help="Choose between 'regular', 'anisotropic', and 'l0' TV norms. The default is 'regular'.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID", default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	(options, args) = parser.parse_args()
	
	if options.output:
		outfile = options.output

	if options.tiltseries:
		infile = options.tiltseries
	elif options.testdata:
		infile = options.testdata
		nslices = options.nslices
		tiltrange = options.tiltrange
	else:
		print "ERROR: You must speficy either --testdata OR --tiltseries."
		exit(1)
	
	if options.imgnum:
		imgnum = int(options.imgnum)
	else: 
		if options.testdata:
			if EMUtil.get_image_count(options.testdata) != 1:
				print "Using the 0th image by default."
		imgnum = 0
	
	if options.norm:
		if options.norm == 'regular':
			print "This instance will utilize the regular (default) TV norm"
		elif options.norm == 'anisotropic':
			print "This instance will utilize an anisotropic TV norm"
		elif options.norm == 'l0':
			print "This instance will utilize an l0 TV norm"
		else:
			print "The option --norm was specified improperly. Please specify either anisotropic, l0, or regular. Regular is specified by default."
			exit(1)
	
	if options.tlt != None:
		fangles = np.asarray([ float( i ) for i in file( options.tlt , "r" ) ])
		tiltangles = fangles.tolist()
		nslices = len( tiltangles )
		pass
	elif ( options.nslices and options.tiltrange ):
		tiltrange = float(options.tiltrange)
		nslices = int(options.nslices)
		print "Using tiltrange from -%s, %s degrees consisting of %i slices."%(options.tiltrange, options.tiltrange, options.nslices)
	elif options.testdata:
		print "You must specify --nslices AND --tiltrange when using --testdata."
		exit(1)
		tiltangles = np.linspace(tiltrange,-1.*tiltrange,nslices).tolist()
	else:
		print "You must specify --tlt when using --tiltseries"
		exit(1)
	
	
	if options.niters:
		niters = int(options.niters)
	
	if options.beta:
		beta = float(options.beta)
		if beta < 0:
			print "Parameter beta (--beta) must be greater than 0."
			exit(1)
	
	if options.subpix:
		subpix = int(options.subpix)
	else: 
		subpix = 1
	
	noisiness = 0.0
	if options.noise:
		if options.noisiness != 0.0:
			noisiness = options.noisiness
	
	if options.output:
		outfile = options.output
	
	if options.verbose > 1: print "e2tvrecon.py"
	logger=E2init(sys.argv,options.ppid)
	
	# Create new output directory for this instance
	options = makepath( options, options.path)
	rootpath = os.getcwd()
	options.path = rootpath + "/" + options.path
	
	# Link original data file to output directory
	if options.verbose > 7: print "Linking input data to instance directory..."
	
	if options.testdata != None:
		pathname = os.path.dirname(os.path.abspath( options.testdata ))
		filename = ntpath.basename( options.testdata )
		linkfrom = pathname + "/" + filename
		linkto = options.path + "/" + filename
		os.symlink( linkfrom, linkto )
	
	if options.tiltseries != None:
		pathname = os.path.dirname(os.path.abspath( options.tiltseries ))
		filename = ntpath.basename( options.tiltseries )
		linkfrom = pathname + "/" + filename
		linkto = options.path + "/" + filename
		os.symlink( linkfrom, linkto )
	
	if options.tlt != None:
		pathname = os.path.dirname(os.path.abspath( options.tlt ))
		filename = ntpath.basename( options.tlt )
		linkfrom = pathname + "/" + filename
		linkto = options.path + "/" + filename
		os.symlink( linkfrom, linkto )
	
	# Get image/projection data
	data, dim = get_data( options, nslices, noisiness, imgnum )
	xlen = dim[0]
	
	# Projection operator and projections data
	if options.verbose > 2: print "Building Projection Operator..."
	projection_operator = build_projection_operator( tiltangles, xlen, nslices, None, subpix, 0, None )
	
	if options.tiltseries:
		projections = data.ravel()[:, np.newaxis]
	else:
		projections = projection_operator * data.ravel()[:, np.newaxis]
	
	if options.verbose > 9: print "Writing Projections to Disk... "
	outpath = options.path + "/" + "projections.hdf"
	for i in range( nslices ):
		from_numpy(projections[i*xlen:(i+1)*xlen]).write_image( outpath, i )
	
	# Reconstruction
	if options.verbose > 2: print "Starting Reconstruction..."
	t1 = time.time()
	recon, energies = fista_tv( options, tiltangles, projections, beta, niters, projection_operator )
	t2 = time.time()
	if options.verbose > 3: print "Reconstruction completed in %s s"%(str(t2-t1))
	
	# Store reconstruction in instance outfile directory
	outpath = options.path + "/" + outfile
	from_numpy( recon[-1] ).write_image( outpath )
	
	if options.fsc != False:
		if options.verbose > 3: print "Generating an FSC plot..."
		fscpath = options.path + "/" + "fsc.txt"
		datapath = options.testdata
		os.popen("e2proc3d.py %s %s --calcfsc %s"%( outpath, fscpath, datapath ))
	
	E2end(logger)
	if options.verbose > 1: print "Exiting"
	return


def get_data( options, nslices, noisiness, imgnum=0 ):
	"""Read an input image as a numpy array return its length in x"""
	if options.testdata:
		if options.verbose > 3:
			print "Generating Projections of %s"%(options.testdata)
		testdata = EMData( options.testdata, imgnum )
		dim = [testdata.get_xsize(), testdata.get_ysize(), testdata.get_zsize()]
		data = testdata.numpy()
	elif options.tiltseries:
		if options.verbose > 3:
			print "Reading Projections from %s"%(options.tiltseries)
		npstack = []
		for i in range( nslices ):
			img = EMData( options.tiltseries, i )
			np_img = img.numpy().copy()
			if (options.noise != False) and (noisiness != 0.0):	# Add Noise to Projections
				if options.verbose > 2: print "Adding Noise to Input Data..."
				np_img += noisiness * np.random.randn(*np_img.shape)
				from_numpy(np_img).write_image(options.path + "noisy_img.hdf")
			npstack.append( np_img )
		data = np.asarray( npstack )
		dim = [img["nx"],img["ny"],img["nz"]]
	return data.astype( np.float32 ), dim


def build_projection_operator( angles, l_x, n_dir=None, l_det=None, subpix=1, offset=0, pixels_mask=None ):
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
	X, Y = _generate_center_coordinates(subpix*l_x)
	X *= 1./subpix
	Y *= 1./subpix
	Xbig, Ybig = _generate_center_coordinates(l_det)
	Xbig *= (l_x - 2*offset) / float(l_det)
	orig = Xbig.min()
	labels = None
	if subpix > 1:
		# Block-group subpixels
		Xlab, Ylab = np.mgrid[:subpix * l_x, :subpix * l_x]
		labels = (l_x * (Xlab / subpix) + Ylab / subpix).ravel()
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
	weights /= subpix**2
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
def fista_tv(options, angles, y, beta, niter, H, verbose=0, mask=None):
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
	for i in range(niter):
		if verbose:
			print i
		eps = 1.e-4
		err = H * x - y
		back_proj = Ht * err
		tmp = x - gamma * back_proj
		if mask is not None:
			tmp2d = np.zeros((l, l))
			tmp2d[mask] = tmp.ravel()
		else:
			tmp2d = tmp.reshape((l, l))
		u_n = tv_denoise_fista(tmp2d, weight=beta*gamma, eps=eps)
		t_new = (1 + np.sqrt(1 + 4 * t_old**2))/2.
		t_old = t_new
		x = u_n + (t_old - 1)/t_new * (u_n - u_old)
		u_old = u_n
		res.append(x)
		data_fidelity_err = 1./2 * (err**2).sum()
		if options.norm == 'anisotropic':
			tv_value = beta * tv_norm_anisotropic(x)
		elif options.norm == 'l0':
			tv_value = beta * tv_l0_norm(x)
		else:
			tv_value = beta * tv_norm(x)
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
	return 0.5 / im_norm * dual_gap


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
