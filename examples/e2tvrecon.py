#!/usr/bin/env python

#
# Author: James Michael Bell, 2014 (jmbell@bcm.edu)
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

from EMAN2 import *
import os
import sys
import scipy
import time
import numpy as np
import matplotlib.pyplot as plt
from numpy import loadtxt
from scipy import sparse
from scipy import ndimage
from scipy import fftpack


def get_usage():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options]
	e2tvrecon.py reconstructs an image from its tomographic projections with an incomplete set of projections. In addition, noise may be added to the projections. In order to reconstruct the original image, we minimize a function that is the sum of a L2 data fit term and the total variation of the image. Proximal iterations using the FISTA scheme are used.
	
	For the original 2D implementation of this algorithm, visit https://github.com/emmanuelle/tomo-tv
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
	parser.add_argument("--testdata", default=None, help="A 2D image to project a number of times (specified by --nslices) and then reconstructed via compressed sensing.")
	parser.add_argument("--tlt", default=None, type=str, help="An imod tlt file containing alignment angles. If specified slices will be inserted using these angles in the IMOD convention")
	parser.add_argument("--nslices", default=120, type=int, help="Specify the number slices into which an image will be projected. Only applicable when using the --testdata option.")
	parser.add_argument("--tiltrange", default='60.0', type=str, help="Specify the range of degrees over which data was collected. If two comma-separated numbers are specified, they will act as a lower and upper bound respectively. For example --tiltrange=50.5 OR --tiltrange=-30.3,65.0.")
	parser.add_argument("--output", default="recon.hdf", help="Output reconstructed tomogram file name.")
	parser.add_argument("--noise",action="store_true",default=False, help="If true, noise will be added to projections.")
	parser.add_argument("--noisiness",default=2, type=int, help="Multiply noise by a specified factor.")
	parser.add_argument("--path",type=str,default='recon',help="Directory in which results will be stored.")
	parser.add_argument("--niters", default=100, type=int, help="Specify the number of iterative reconstructions to complete before returning the final reconstructed volume.")
	#parser.add_argument("--crossval", default=None, help="Use cross validaton to specify the parameter 'beta'. Input 0 for the best beta as determined by cross-validation or 1 for the best beta for segmentation as compared to a complete data set (i.e. for use with test data).")
	parser.add_argument("--beta", default=20.0, type=float, help="Specify the total-variation regularization weight parameter 'beta' without performing cross validation.")
	parser.add_argument("--subpix", default=1, type=int, help="Specify the number of linear subdivisions used to compute the projection of one image pixel onto a detector pixel.")
	parser.add_argument("--fsc",action="store_true",default=False, help="If true, an fourier shell correlation plot will be generated comparing the input and output data.")
	parser.add_argument("--plots",action="store_true",default=False, help="If true, python plots will be shown.from_numpy")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID", default=-1)
	(options, args) = parser.parse_args()

	# parse options
	if options.output : outfile = options.output
	if options.tiltseries and options.testdata:
		print "A tiltseries and testdata may not be specified simultaneously."
		exit(1)
	
	if options.tiltseries : 
		infile = options.tiltseries
	elif options.testdata : infile = options.testdata
	else: 
		print "You must speficy either --testdata OR --tiltseries"
		exit(1)
	
	if options.imgnum != None: imgnum = options.imgnum
	else:
		imgnum = 0
	
	if options.tiltrange == None: options.tiltrange = np.pi / 3.0
	else:
		options.tiltrange = map(float, options.tiltrange.split(','))
			
	if options.niters : niters = int(options.niters)
	if options.beta : 
		beta = float(options.beta)
		if beta == 0:
			print "Parameter beta cannot equal 0."
			exit(1)
			
	if options.subpix : subpix = int(options.subpix)
	else: subpix = 1
		
	if options.tlt: tltfile = options.tlt
	
	if options.noisiness : noisiness = options.noisiness
	else: noisiness = 2.
		
	if options.nslices: nslices = options.nslices
	elif options.tlt:
		angles = get_angles( tiltfile )
		nslices = len( angles )
	elif options.tiltseries and not options.tlt:
		nslices = EMUtil.get_image_count( infile )
	else:
		print "You must specify --nslices explicitly or implicitly by supplying a tiltseries or tlt file"
		exit(1)
		
	if options.output : outfile = options.output
	
	if options.verbose > 1: print "e2tvrecon.py"
	logger=E2init(sys.argv,options.ppid)
	
	# Create new output directory for this instance
	options = makepath( options, options.path)
	rootpath = os.getcwd()
	options.path = rootpath + "/" + options.path
	
	# Link original data file to output directory
	if options.testdata != None:
		os.symlink(infile, options.path + "/input.hdf")
	
	#if options.crossval != None:
	#	betas = get_best_betas( infile, nslices, niters)
	#	if options.crossval == 1:
	#		beta = betas[1]
	#	else:
	#		beta = betas[0]
	
	# Get Tomogram data
	img, dim, options = get_tomo_data( options, imgnum )
	xlen = dim[0]
	#ylen = dim[1]
	#zlen = dim[2]
	
	# Projection operator and projections data, with noise
	projection_operator = build_projection_operator( options, xlen, nslices, None, subpix, 0, None)
	projections = projection_operator * img.ravel()[:, np.newaxis]
	
	# Generate stack of projections
	outpath = options.path + "/prjs.hdf"
	for i in range( nslices ):
		from_numpy(projections[i*dim[0]:(i+1)*dim[0]]).write_image( outpath, i )

	if options.noise != False:	# add Noise to Projections
		projections += 2*np.random.randn(*projections.shape)
		# Generate stack of noisy projections	
		outpath = options.path + "/noisy-prjs.hdf"
		for i in range( nslices ):
			from_numpy(projections[i*dim[0]:(i+1)*dim[0]]).write_image( outpath, i )
			
	# Reconstruction
	t1 = time.time()
	recon, energies = fista_tv(projections, beta, niters, projection_operator) 
	t2 = time.time()
	if options.verbose > 2: print "reconstruction completed in %f s" %(t2 - t1)

	# Fraction of errors of segmented image wrt 'complete' data set
	err = [np.abs(img - (reconstruction > 0.5)).mean() for reconstruction in recon]
	
	# Store reconstruction in instance outfile directory
	outpath = options.path + "/" + outfile
	from_numpy(recon[-1]).write_image( outpath )
	
	if options.fsc != False:
		fscpath = options.path + "/fsc.txt"
		datapath = options.testdata
		os.popen("e2proc3d.py %s %s --calcfsc %s"%(outpath, fscpath, datapath))

	if options.plots != False: # Display comparison of images
		plt.figure('Qualitative Comparison of Data and Reconstruction')
		plt.subplot(221)
		plt.imshow(img, interpolation='nearest', vmin=0, vmax=1)
		plt.title('original data (%i x %i x %i)'%(dim[0],dim[1],dim[2]))
		plt.axis('off')
		plt.subplot(222)
		plt.imshow(recon[-1], interpolation='nearest', vmin=0, vmax=1)
		plt.title('Reconstruction \n %i iterations, B = %f'%(niters,beta))
		plt.axis('off')
		plt.subplot(223)
		plt.loglog(energies, 'o')
		plt.xlabel('iteration number')
		plt.title('energy')
		plt.subplot(224)
		plt.loglog(err, 'o')
		plt.xlabel('iteration number')
		plt.title('error fraction')
		plt.show()
	
	E2end(logger)
	if options.verbose > 1: print "Exiting"
	return
	
def tv_norm(im):
	"""Compute the (isotropic) TV norm of an image"""
	grad_x1 = np.diff(im, axis=0)
	grad_x2 = np.diff(im, axis=1)
	#grad_x3 = np.diff(im, axis=2)
	#return np.sqrt(grad_x1[:, :-1, :-1]**2 + grad_x2[:-1, :, :-1]**2 + grad_x3[:-1,:-1,:]**2).sum()
	return np.sqrt(grad_x1[:, :-1]**2 + grad_x2[:-1, :]**2).sum()


def tv_norm_anisotropic(im):
	"""Compute the anisotropic TV norm of an image"""
	grad_x1 = np.diff(im, axis=0)
	grad_x2 = np.diff(im, axis=1)
	#grad_x3 = np.diff(im, axis=2)
	#return np.abs(grad_x1[:, :-1, :-1]).sum() + np.abs(grad_x2[:-1, :, :-1]).sum() + np.abs(grad_x3[:-1, :-1, :]).sum()
	return np.abs(grad_x1[:, :-1]).sum() + np.abs(grad_x2[:-1, :]).sum()


def fista_tv(y, beta, niter, H, verbose=0, mask=None):
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

	E(x) = 1/2 || H x - y ||^2 + beta TV(x) = f(x) + beta TV(x)

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
	n_angles = n_meas / l
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
		tv_value = beta * tv_norm(x)
		energy = data_fidelity_err + tv_value
		energies.append(energy)
		if mask is not None:
			x = x[mask][:, np.newaxis]
		else:
			x = x.ravel()[:, np.newaxis]
	return res, energies


def get_tomo_data( options, imgnum=0 ):
	"""Read a tomogram as a numpy array return its dimensions"""
	if options.tiltseries != None:
		tomo = EMData( options.tiltseries , imgnum )
		outpath = options.path + "/input.hdf"
		tomo.write_image( outpath )
		options.testdata = outpath
	else:
		tomo = EMData( options.testdata , imgnum )
	dim = [tomo.get_xsize(), tomo.get_ysize(), tomo.get_zsize()]
	tomo_array = tomo.numpy()
	return tomo_array.astype(np.float32), dim, options


def make_tlt( options ):	
	# form a 'complete' list of tilt angles
	options.range
	lower_angles=[]
	while angle < lower_bound:		# generate angles below data
		lower_angles.append( angle )
		angle = angle + tltstep
		
	upper_angles=[]
	angle = upper_bound + 1.0	
	while angle <= 90.0:			# generate angles above data
		upper_angles.append( angle )
		angle = angle + tltstep
	
	new_angles = lower_angles + angles + upper_angles
	# generate a new tilt angles file in wedgecomp_0# directory
	
	options.tlt = options.path + "/" + os.path.splitext(options.tlt)[0] + "_new.tlt"
	
	with open( options.tlt ,"w") as newtlt:
		for a in angles:
			newtlt.write( str( a ) + "\n")
			
	return options, lower_bound, upper_bound


def tv_l0_norm( img ):
	"""Compute the (isotropic) TV norm of a 2D image"""
	grad_x1 = np.diff(im, axis=0)
	grad_x2 = np.diff(im, axis=1)
	return (grad_x1[:, :-1]**2 + grad_x2[:-1, :]**2 > 0).mean()


def compute_sparsity( img ):
	l_x = len(img)
	X, Y = np.ogrid[:l_x, :l_x]
	mask = ((X - l_x/2)**2 + (Y - l_x/2)**2 <= (l_x/2)**2)
	grad1 = ndimage.morphological_gradient(img, footprint=np.ones((3, 3)))
	grad2 = ndimage.morphological_gradient(img, footprint=ndimage.generate_binary_structure(2, 1))
	return (grad1[mask] > 0).mean(), (grad2[mask] > 0).mean()


def generate_synthetic_data(l_x=128, seed=None, crop=True, n_pts=25):
	if seed is None:
		seed = 0
	# Fix the seed for reproducible results
	rs = np.random.RandomState(seed)
	x, y = np.ogrid[:l_x, :l_x]
	mask = np.zeros((l_x, l_x))
	points = l_x * rs.rand(2, n_pts)
	mask[(points[0]).astype(np.int), (points[1]).astype(np.int)] = 1
	mask = ndimage.gaussian_filter(mask, sigma=l_x / (4. * np.sqrt(n_pts)))
	# Limit the non-zero data to a central circle
	if crop:
		mask_outer = (x - l_x / 2) ** 2 + (y - l_x / 2) ** 2 < (l_x / 2) ** 2
		mask = np.logical_and(mask > mask.mean(), mask_outer)
	else:
		mask = mask > mask.mean()
	return mask.astype(np.float32)


def makepath(options, stem=''):
	if options.verbose > 5: print "makepath function called"
	if options.path and ("/" in options.path or "#" in options.path):
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. Please edit your --path argument accordingly."
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
			jobtag=''
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


def get_angles( options ):
	angles = []
	if options.tlt:
		angles=[ float( i ) for i in file( options.tlt , "r" ) ]
	elif options.tiltseries:
		n = EMUtil.get_image_count( options.tiltseries )
		for i in range( n ):
			headderangle = EMData( options.tiltseries, i, True )['tiltangle']
			angles.append( headderangle )
	else:
		print "No tlt file or tiltseries. Returning empty array."
	return angles


# --------------- Data projection operator  --------------------
def build_projection_operator( options, l_x, n_dir=None, l_det=None, subpix=1, offset=0, pixels_mask=None ):
	"""
	Compute the tomography design matrix.
		
	Parameters
	----------
	
	l_x : int
		linear size of image array

	n_dir : int, default l_x
		number of angles at which projections are acquired. n_dir projection angles are regularly spaced between 0 and 180.
		
	l_det : int, default is l_x
		number of pixels in the detector. If l_det is not specified, we suppose that l_det = l_x.
	
	subpix : int, default 1
		number of linear subdivisions used to compute the projection of one image pixel onto a detector pixel.
    
	offset : int, default 0
		width of the strip of image pixels not covered by the detector.

	pixels_mask : 1-d ndarray of size l_x**2
		mask of pixels to keep in the matrix (useful for to removeing pixels inside or outside of a circle, for example)
    
    Returns
	-------
	p : sparse matrix of shape (n_dir l_x, l_x**2), in csr format Tomography design matrix. The csr (compressed sparse row) allows for efficient subsequent matrix multiplication. The dtype of the elements is float32, in order to save memory.	
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
	tiltrange = options.tiltrange
	if options.tlt != None:
		angles = get_angles( options )
	elif len(tiltrange) > 1:
		angles = np.linspace(tiltrange[0], tiltrange[1], n_dir, endpoint=False)
	else:
		angles = np.linspace(-1*tiltrange[0], tiltrange[0], n_dir, endpoint=False)
	weights, data_inds, detector_inds = [], [], []
	# Indices for data pixels. For each data, one data pixel
	# will contribute to the value of two detector pixels.
	for i, angle in enumerate(angles):
		# rotate data pixels centers
		Xrot = np.cos(angle) * X - np.sin(angle) * Y
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


def _weights(x, dx=1, orig=0, ravel=True, labels=None):
	"""
	Compute linear interpolation weights for projection array `x`
	and regularly spaced detector pixels separated by `dx` and
	starting at `orig`.
	"""
	if ravel:
		x = np.ravel(x)
	floor_x = np.floor((x - orig) / dx).astype(np.int32)
	alpha = ((x - orig - floor_x * dx) / dx).astype(np.float32)
	inds = np.hstack((floor_x, floor_x + 1))
	weights = np.hstack((1 - alpha, alpha))
	data_inds = np.arange(x.size, dtype=np.int32)
	data_inds = np.hstack((data_inds, data_inds))
	if labels is not None:
		data_inds = np.hstack((labels, labels))
		w = np.histogram2d(data_inds, inds,
			bins=(np.arange(data_inds.max()+1.5), np.arange(inds.max()+1.5)),
			weights=weights)[0]
		data_inds, inds = np.argwhere(w>0).T
		weights = w[w>0]
	return inds, data_inds, weights


def _weights_nn(x, dx=1, orig=0, ravel=True):
	"""
	Nearest-neighbour interpolation
	"""
	if ravel:
		x = np.ravel(x)
	floor_x = np.floor(x - orig)
	return floor_x.astype(np.float32)


def _generate_center_coordinates_3d(dim):
	"""
	Compute the coordinates of pixels centers for an image of
	linear size l_x
	"""
	lx = float(dim[0])
	lx = float(dim[1])
	lx = float(dim[2])
	X, Y, Z= np.mgrid[:lx, :ly, :lz]
	center_x = lx / 2.
	center_y = ly / 2.
	center_z = lz / 2.
	X += 0.5 - center_x
	Y += 0.5 - center_y
	Z += 0.5 - center_z
	return X, Y, Z


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


# ----------------- Direct projection method -------------------------
# (without computing explicitely the design matrix)
def back_projection(projections):
	"""
	Back-projection (without filtering)

	Parameters
	----------
	projections: ndarray of floats, of shape n_dir x l_x
		Each line of projections is the projection of a data image
		acquired at a different angle. The projections angles are
		supposed to be regularly spaced between 0 and 180.

	Returns
	-------
	recons: ndarray of shape l_x x l_x
		Reconstructed array

	Notes
	-------
	A linear interpolation is used when rotating the back-projection.
	This function uses ``scipy.ndimage.rotate`` for the rotation.
	"""
	n_dir, l_x = projections.shape
	recons = np.zeros((l_x, l_x), dtype=np.float)
	angles = np.linspace(0, 180, n_dir, endpoint=False)
	for angle, line in zip(angles, projections):
		# BP: repeat the detector line along the direction of the beam
		tmp = np.tile(line[:, np.newaxis], (1, l_x))
		# Rotate the back-projection of the detector line, and add
		# it to the reconstructed image
		recons += ndimage.rotate(tmp, -angle, order=1, reshape=False)
	return recons


def projection(im, n_dir=None, interpolation='nearest'):
	"""
	Tomography projection of an image along n_dir directions.

	Parameters
	----------
	im : ndarray of square shape l_x x l_x
		Image to be projected

	n_dir : int
		Number of projection angles. Projection angles are regularly spaced
		between 0 and 180.

	interpolation : str, {'interpolation', 'nearest'}
		Interpolation method used during the projection. Default is
		'nearest'.

	Returns
	-------
	projections: ndarray of shape n_dir x l_x
		Array of projections.

	Notes
	-----
	The centers of the data pixels are projected onto the detector, then
	the contribution of a data pixel to detector pixels is computed
	by nearest neighbor or linear interpolation. The function 
	np.bincount`` is used to compute the projection, with weights
	corresponding to the values of data pixels, multiplied by interpolation
	weights in the case of linear interpolation.
	"""
	l_x = len(im)
	if n_dir is None:
		n_dir = l_x
	im = im.ravel()
	projections = np.empty((n_dir, l_x))
	X, Y = _generate_center_coordinates(l_x)
	angles = np.linspace(0, np.pi, n_dir, endpoint=False)
	for i, angle in enumerate(angles):
		Xrot = np.cos(angle) * X - np.sin(angle) * Y 
		if interpolation == 'nearest':
			inds = _weights_nn(Xrot, dx=1, orig=X.min())
			mask = inds>= 0
			w = im[mask]
		elif interpolation == 'linear':
			inds, _, w = _weights(Xrot, dx=1, orig=X.min())
			w[:l_x**2] *= im
			w[l_x**2:] *= im
			mask = inds >= 0
			w = w[mask]
		projections[i] = np.bincount(inds[mask].astype(np.int), weights=w)[:l_x]
	return projections


# -----------------Filtered back-projection----------------------
def filter_projections(proj_set, reg=False):
	"""
	Ramp filter used in the filtered back projection.
	We use zero padding.

	Parameters
	----------
	proj_set: 2-d ndarray
		each line is one projection (1 line of the detector) to be filtered

	Returns
	-------

	res: 2-d ndarray
		filtered projections

	Notes
	-----

	We use zero padding. However, we do not use any filtering (hanning, etc.)
	in the FFT yet.
	"""
	nb_angles, l_x = proj_set.shape
	#Assume l is even for now
	ramp = 1./l_x * np.hstack((np.arange(l_x), np.arange(l_x, 0, -1)))
	return fftpack.ifft(ramp * fftpack.fft(proj_set, 2*l_x, axis=1), axis=1)[:,:l_x]


def rank_order(image):
	"""Return an image of the same shape where each pixel is the index of the pixel value in the ascending order of the unique values of `image`, aka the rank-order value.
	
	Parameters
	----------
	image: ndarray
	
	Returns
	-------
	labels: ndarray of type np.uint32, of shape image.shape
		New array where each pixel has the rank-order value of the corresponding pixel in `image`. Pixel values are between 0 and n - 1, where n is the number of distinct unique values in ''image`.
	
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
		Gradient of the image: the i-th component along the first axis is the gradient along the i-th axis of the original array img
	"""
	shape = [img.ndim, ] + list(img.shape)
	gradient = np.zeros(shape, dtype=img.dtype)
	# 'Clever' code to have a view of the gradient with dimension i stop
	# at -1
	slice_all = [0, slice(None, -1),]
	for d in range(img.ndim):
		gradient[slice_all] = np.diff(img, axis=d)
		slice_all[0] = d + 1
		slice_all.insert(1, slice(None))
	return gradient


def _projector_on_dual(grad):
	"""
	modifies in place the gradient to project iton the L2 unit ball
	"""
	norm = np.maximum(np.sqrt(np.sum(grad**2, 0)), 1.)
	for grad_comp in grad:
		grad_comp /= norm
	return grad


def dual_gap(im, new, gap, weight):
	"""
	dual gap of total variation denoising
	see "Total variation regularization for fMRI-based prediction of behavior", by Michel et al. (2011) for a derivation of the dual gap
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
	
	
def get_tomo_data( options, imgnum=0 ):
	"""Read a tomogram as a numpy array return its dimensions"""
	if options.tiltseries != None:
		tomo = EMData( options.tiltseries , imgnum )
		outpath = options.path + "/input.hdf"
		tomo.write_image( outpath )
		options.testdata = outpath
	else:
		tomo = EMData( options.testdata , imgnum )
	dim = [tomo.get_xsize(), tomo.get_ysize(), tomo.get_zsize()]
	tomo_array = tomo.numpy()
	return tomo_array.astype(np.float32), dim, options


def make_tlt( options ):	
	# form a 'complete' list of tilt angles
	options.range
	lower_angles=[]
	while angle < lower_bound:		# generate angles below data
		lower_angles.append( angle )
		angle = angle + tltstep
		
	upper_angles=[]
	angle = upper_bound + 1.0	
	while angle <= 90.0:			# generate angles above data
		upper_angles.append( angle )
		angle = angle + tltstep
	
	new_angles = lower_angles + angles + upper_angles
	# generate a new tilt angles file in wedgecomp_0# directory
	
	options.tlt = options.path + "/" + os.path.splitext(options.tlt)[0] + "_new.tlt"
	
	with open( options.tlt ,"w") as newtlt:
		for a in angles:
			newtlt.write( str( a ) + "\n")
			
	return options, lower_bound, upper_bound


def tv_l0_norm( img ):
	"""Compute the (isotropic) TV norm of a 2D image"""
	grad_x1 = np.diff(im, axis=0)
	grad_x2 = np.diff(im, axis=1)
	return (grad_x1[:, :-1]**2 + grad_x2[:-1, :]**2 > 0).mean()


def compute_sparsity( img ):
	l_x = len(img)
	X, Y = np.ogrid[:l_x, :l_x]
	mask = ((X - l_x/2)**2 + (Y - l_x/2)**2 <= (l_x/2)**2)
	grad1 = ndimage.morphological_gradient(img, footprint=np.ones((3, 3)))
	grad2 = ndimage.morphological_gradient(img, footprint=ndimage.generate_binary_structure(2, 1))
	return (grad1[mask] > 0).mean(), (grad2[mask] > 0).mean()


def generate_synthetic_data(l_x=128, seed=None, crop=True, n_pts=25):
	if seed is None:
		seed = 0
	# Fix the seed for reproducible results
	rs = np.random.RandomState(seed)
	x, y = np.ogrid[:l_x, :l_x]
	mask = np.zeros((l_x, l_x))
	points = l_x * rs.rand(2, n_pts)
	mask[(points[0]).astype(np.int), (points[1]).astype(np.int)] = 1
	mask = ndimage.gaussian_filter(mask, sigma=l_x / (4. * np.sqrt(n_pts)))
	# Limit the non-zero data to a central circle
	if crop:
		mask_outer = (x - l_x / 2) ** 2 + (y - l_x / 2) ** 2 < (l_x / 2) ** 2
		mask = np.logical_and(mask > mask.mean(), mask_outer)
	else:
		mask = mask > mask.mean()
	return mask.astype(np.float32)


def makepath(options, stem=''):
	if options.verbose > 5: print "makepath function called"
	if options.path and ("/" in options.path or "#" in options.path):
		print "Path specifier should be the name of a subdirectory to use in the current directory. Neither '/' or '#' can be included. Please edit your --path argument accordingly."
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
			jobtag=''
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
	

# BROKEN! Sorry... :(
def get_best_betas( imgpath, nslices=None, niters=400 ):
	x, dim = get_tomo_data( imgpath )
	l = dim[0]
	if nslices == None:
		n_dir = l / 3.
	else: 
		n_dir = nslices
		
	def rec_error(beta, niters):
		"""cross-validation"""
		res, energies = fista_tv(y1, beta, niters, H)
		yres = H * res[-1].ravel()[:, np.newaxis]
		return (((yres - y2)**2).mean()), res[-1], energies

	# Projection operator and projections data, with 2 realizations of the noise
	H = build_projection_operator(l, n_dir)
	y = H * x.ravel()[:, np.newaxis]
	np.random.seed(0)
	y1 = y + 2*np.random.randn(*y.shape)  # 1st realization
	y2 = y + 2*np.random.randn(*y.shape)  # 2nd realization
	
	# Range of beta parameter
	betas = 2**np.arange(2, 6, 0.25)

	results = Parallel(n_jobs=-1)(delayed(rec_error)(beta, niters) for beta in betas)
	errors = [res[0] for res in results]
	images = [res[1] for res in results]
	energies = [res[2] for res in results]
	
	# Segmentation compared to ground truth
	segmentation_error = [np.abs((image > 0.5) - x).mean() for image in images]
	
	print "best beta from cross-validation %f" %(betas[np.argmin(errors)])
	print "best beta for segmentation compared to ground truth %f"%(betas[np.argmin(segmentation_error)])
	best_betas = [betas[np.argmin(errors)], betas[np.argmin(segmentation_error)]]
	return best_betas
	

if __name__=="__main__":
	main()
