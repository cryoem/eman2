#!/usr/bin/env python

#
# Author: James Michael Bell, 8/9/15 (jmbell@bcm.edu)
# Copyright (c) 2015- Baylor College of Medicine
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

from __future__ import division, print_function
from EMAN2 import *
from math import *
from numpy.fft import fft2, fftshift, ifft2
from scipy import ndimage
import math
import matplotlib.pyplot as plt
import numpy
import numpy as np
import scipy as sp
try: import scipy.ndimage.interpolation as ndii
except ImportError: import ndimage.interpolation as ndii

def main():
	maxshift = 20
	noise = 2.5
	fname = "/Users/jmbell/workspace/data/centering/ribosome.hdf"
	
	n = np.random.random_integers(0,EMUtil.get_image_count(fname))
	
	im0 = EMData(fname,n) #test_image(0)
	im0 = im0.process('normalize.edgemean')
	im0 = im0.numpy().astype(np.float32)
	
	im1 = EMData(fname,n) #test_image(0)
	im1 = im1.process('normalize.edgemean')
	t = [np.random.random_integers(-maxshift,maxshift), np.random.random_integers(-maxshift,maxshift)]
	im1.translate(t[0],t[1],0)
	im1=im1.process('math.addnoise',{'noise':noise})
	im1 = im1.numpy().astype(np.float32)
	
	im2, scale, angle, (t0, t1), ir = similarity(im0, im1)
	
	lbls = ['Image ({},{})'.format(round(t[0],0),round(t[1],0)),'Target (0,0)', 'Phase Correlation', 'Aligned ({},{})\nScale: {}, Angle: {}'.format(t[0]-t1,t[1]-t0,scale,angle)]
	fig, axar = plt.subplots(2,2)
	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])	
	ar_plts = [im1,im0,ir,im2]
	for i,ax in enumerate(axar.flat):
		im = ax.imshow(ar_plts[i], interpolation='nearest', origin='lower',cmap=plt.cm.Greys_r)
		ax.grid(False)
		ax.set_yticks([])
		ax.set_xticks([])
		ax.set_xlabel(lbls[i])
	fig.colorbar(im, cax=cbar_ax,cmap=plt.cm.Greys_r)
	
	plt.show()

def translation(im0, im1):
	"""Return translation vector to register images."""
	shape = im0.shape
	f0 = fft2(im0)
	f1 = fft2(im1)
	ir = abs(ifft2((f0 * f1.conjugate()) / (abs(f0) * abs(f1))))
	t0, t1 = numpy.unravel_index(numpy.argmax(ir), shape)
	if t0 > shape[0] // 2: t0 -= shape[0]
	if t1 > shape[1] // 2: t1 -= shape[1]
	return [t0, t1]

def similarity(im0, im1):
	"""Return similarity transformed image im1 and transformation parameters.

	Transformation parameters are: isotropic scale factor, rotation angle (in
	degrees), and translation vector.

	A similarity transformation is an affine transformation with isotropic
	scale and without shear.

	Limitations:
	Image shapes must be equal and square.
	All image areas must have same scale, rotation, and shift.
	Scale change must be less than 1.8.
	No subpixel precision.

	"""
	if im0.shape != im1.shape: raise ValueError("Images must have same shapes.")
	elif len(im0.shape) != 2: raise ValueError("Images must be 2 dimensional.")

	f0 = fftshift(abs(fft2(im0)))
	f1 = fftshift(abs(fft2(im1)))

	h = highpass(f0.shape)
	f0 *= h
	f1 *= h
	del h

	f0, log_base = logpolar(f0)
	f1, log_base = logpolar(f1)

	f0 = fft2(f0)
	f1 = fft2(f1)
	r0 = abs(f0) * abs(f1)
	ir = abs(ifft2((f0 * f1.conjugate()) / r0))
	i0, i1 = numpy.unravel_index(numpy.argmax(ir), ir.shape)
	angle = 180.0 * i0 / ir.shape[0]
	scale = log_base ** i1

	if scale > 1.8:
		ir = abs(ifft2((f1 * f0.conjugate()) / r0))
		i0, i1 = numpy.unravel_index(numpy.argmax(ir), ir.shape)
		angle = -180.0 * i0 / ir.shape[0]
		scale = 1.0 / (log_base ** i1)
		if scale > 1.8:
			raise ValueError("Images are not compatible. Scale change > 1.8")

	if angle < -90.0: angle += 180.0
	elif angle > 90.0: angle -= 180.0

	im2 = ndii.zoom(im1, 1.0/scale)
	im2 = ndii.rotate(im2, angle)

	if im2.shape < im0.shape:
		t = numpy.zeros_like(im0)
		t[:im2.shape[0], :im2.shape[1]] = im2
		im2 = t
	elif im2.shape > im0.shape:
		im2 = im2[:im0.shape[0], :im0.shape[1]]

	f0 = fft2(im0)
	f1 = fft2(im2)
	ir = abs(ifft2((f0 * f1.conjugate()) / (abs(f0) * abs(f1))))
	t0, t1 = numpy.unravel_index(numpy.argmax(ir), ir.shape)

	if t0 > f0.shape[0] // 2: t0 -= f0.shape[0]
	if t1 > f0.shape[1] // 2: t1 -= f0.shape[1]

	im2 = ndii.shift(im2, [t0, t1])

	# correct parameters for ndimage's internal processing
	if angle > 0.0:
		d = int((int(im1.shape[1] / scale) * math.sin(math.radians(angle))))
		t0, t1 = t1, d+t0
	elif angle < 0.0:
		d = int((int(im1.shape[0] / scale) * math.sin(math.radians(angle))))
		t0, t1 = d+t1, d+t0
	scale = (im1.shape[1] - 1) / (int(im1.shape[1] / scale) - 1)

	return im2, scale, angle, [-t0, -t1], ir

def similarity_matrix(scale, angle, vector):
	"""Return homogeneous transformation matrix from similarity parameters.

	Transformation parameters are: isotropic scale factor, rotation angle (in
	degrees), and translation vector (of size 2).

	The order of transformations is: scale, rotate, translate.

	"""
	S = numpy.diag([scale, scale, 1.0])
	R = numpy.identity(3)
	angle = math.radians(angle)
	R[0, 0] = math.cos(angle)
	R[1, 1] = math.cos(angle)
	R[0, 1] = -math.sin(angle)
	R[1, 0] = math.sin(angle)
	T = numpy.identity(3)
	T[:2, 2] = vector
	return numpy.dot(T, numpy.dot(R, S))

def logpolar(image, angles=None, radii=None):
	"""Return log-polar transformed image and log base."""
	shape = image.shape
	center = shape[0] / 2, shape[1] / 2
	if angles is None: angles = shape[0]
	if radii is None: radii = shape[1]
	theta = numpy.empty((angles, radii), dtype=numpy.float64)
	theta.T[:] = -numpy.linspace(0, numpy.pi, angles, endpoint=False)
	#d = radii
	d = numpy.hypot(shape[0]-center[0], shape[1]-center[1])
	log_base = 10.0 ** (math.log10(d) / (radii))
	radius = numpy.empty_like(theta)
	radius[:] = numpy.power(log_base, numpy.arange(radii,
												   dtype=numpy.float64)) - 1.0
	x = radius * numpy.sin(theta) + center[0]
	y = radius * numpy.cos(theta) + center[1]
	output = numpy.empty_like(x)
	ndii.map_coordinates(image, [x, y], output=output)
	return output, log_base

def highpass(shape):
	"""Return highpass filter to be multiplied with fourier transform."""
	x = numpy.outer(
		numpy.cos(numpy.linspace(-math.pi/2., math.pi/2., shape[0])),
		numpy.cos(numpy.linspace(-math.pi/2., math.pi/2., shape[1])))
	return (1.0 - x) * (2.0 - x)

if __name__ == "__main__":
	main()
