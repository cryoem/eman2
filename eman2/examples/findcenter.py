#!/usr/bin/env python

#
# Authors: James Michael Bell, 06/03/2015
# Copyright (c) 2015 Baylor College of Medicine
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
import os
import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import scipy
import scipy.signal


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """findcenter.py finds the center of images in a stack.
	"""

	parser = EMArgumentParser(usage=usage, version=EMANVERSION)

	parser.add_argument("--path",type=str,default=None,help="Specify the path to the particle stack you wish to center.",required=True)
	parser.add_argument("--output",type=str,default='centered.hdf',help="Specify the path to the particle stack you wish to center.")
	parser.add_argument("--predisplay",action="store_true",default=False,help="Display particles before filtering.")
	parser.add_argument("--postdisplay",action="store_true",default=False,help="Choose not to display particles after filtering.")
	parser.add_argument("--nowrite", action="store_true", default=False, help="Choose not to display particles after filtering.")
	parser.add_argument("--nobinarize", action="store_true", default=False, help="Choose not to binarize particles after filtering.")
	parser.add_argument("--deviations", type=float, help="Set the number of standard deviations to exclude when binarizing.",default=3.0)
	parser.add_argument("--nopreprocess",action="store_true",default=False,help="Choose not to apply filtering to particles.")
	parser.add_argument("--noedgenorm",action="store_true",default=False,help="Display particles before filtering.")
	parser.add_argument("--filter",choices=['wavelet','bilateral','rbf','bandpass'],default='bilateral',help="Select the type of wavelet you wish to utilize.")
	parser.add_argument("--wavelet",choices=['harr','daub','bspl'],default='daub',help="Select the type of wavelet you wish to utilize.")
	parser.add_argument("--kernel_mu", type=int, help="Set the mean (mu) of the gaussian kernel for rbf filtering.",default=64)
	parser.add_argument("--kernel_sigma", type=int, help="Set the standard deviation (sigma) of the gaussian kernel for rbf filtering.",default=8)
	parser.add_argument("--distance_sigma", type=float, help="Set distance_sigma value for bilateral filtering.",default=1.0)
	parser.add_argument("--value_sigma", type=float, help="Set the value_sigma parameter for bilateral filtering.",default=1.0)
	parser.add_argument("--half_width", type=int, help="Set the half_width parameter for bilateral filtering.",default=1)
	parser.add_argument("--order",type=int,choices=[2,4,6,8,10,12,14,16,18,20,103,105,202,204,206,208,301,303,305,307,309],default=4,help="Select the type of wavelet you wish to utilize.")
	parser.add_argument("--meanshrink", type=int, help="Specify by how many times you wish to mean shrink the particles.",default=1)
	parser.add_argument("--ntimes", type=int, help="Specify by how many times you wish to apply the wavelet filter.",default=1)
	parser.add_argument("--step", type=int, help="Specify the step size for computing projections (1 every N degrees). The default 30.",default=30)
	parser.add_argument("--threshold", type=float, help="Set the minval to threshold in wavelet space or after filtering via a radial basis function.",default=10.0)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if not os.path.isfile(options.path):
		print("The image stack you specified does not exist. Please check the --path parameter and try again.")
		sys.exit(1)

	if options.wavelet == 'harr':
		if options.order != 2:
			print("Harr wavelets must be of order 2.")
			sys.exit(1)
	if options.wavelet == 'bspl':
		if len(str(options.order)) != 3:
			print("B-Spline wavelets must have a 3 digit order. Options include: 103,105,202,204,206,208,301,303,305,307, and 309.")
			sys.exit(1)
	if options.wavelet == 'daub':
		if options.order < 4 or options.order > 20:
			print("Daubche wavelets must be of even order between 4 and 20.")
			sys.exit(1)

	if options.output.split('.')[-1] != 'hdf':
		print("Replacing your output file extension ('.{}') with '.hdf'".format(options.output.split('.')[-1]))
		options.output = options.output.split('.')[-2] + '.hdf'

	if os.path.isfile(options.output):
		os.remove(options.output)
	if os.path.isfile(options.path[:-4]+'_proc.hdf'):
		os.remove(options.path[:-4]+'_proc.hdf')

	hdr = EMData(options.path,0,True)
	nx = hdr.get_xsize() # we assume square boxes so we won't bother with ny
	for i in xrange(options.meanshrink):
		nx = nx / 2
	if nx < 8:
		print("You're trying to shrink the image a too much. Try using a lower value for the --meanshrink parameter.")
		sys.exit(1)
	if options.meanshrink != 1 and options.meanshrink % 2 != 0:
		print("The --meanshrink parameter must be an even number.")
		sys.exit(1)

	if options.filter == 'rbf':
		gauss1d = scipy.signal.gaussian(options.kernel_mu, options.kernel_sigma)
		kernel = np.outer(gauss1d,gauss1d)

	nimgs = EMUtil.get_image_count(options.path)
	if options.predisplay: display([EMData(options.path,i) for i in xrange(nimgs)])
	if options.verbose > 2: print("Applying a {} filter {} times to {}.".format(options.filter, options.ntimes, os.path.basename(options.path)))

	if not options.nopreprocess:
		for i in xrange(nimgs):

			if options.verbose > 5: print("Pre-processing image {}/{}".format(i+1,nimgs))

			ptcl = EMData(options.path,i)

			if options.meanshrink > 1: ptcl.process_inplace('math.meanshrink',{'n':options.meanshrink})
			if not options.noedgenorm: ptcl.process_inplace('normalize.edgemean')

			for i in range(options.ntimes):
				if options.filter == 'bilateral':
					params = {'distance_sigma':options.distance_sigma, 'half_width':options.half_width, 'niter':1, 'value_sigma':options.value_sigma}
					ptcl.process_inplace('filter.bilateral',params)
				if options.filter == 'wavelet':
					ptcl.process_inplace('basis.wavelet',{'dir':1,'ord':options.order,'type':options.wavelet})
					ptcl.process_inplace('threshold.belowtozero',{'minval':options.threshold})
					ptcl.process_inplace('basis.wavelet',{'dir':-1,'ord':options.order,'type':options.wavelet})
				if options.filter == 'rbf':
					ptcl_array = ptcl.numpy()
					conv2d = scipy.signal.convolve2d(ptcl_array, kernel, boundary='symm', mode='same')
					ptcl = from_numpy(np.float32(conv2d))
					ptcl.process_inplace('threshold.belowtozero',{'minval':options.threshold})
				if options.filter == 'bandpass':
					apix = ptcl.get_attr('ctf').apix
					sigma = ptcl.get_attr('sigma') / 20.0
					ptcl.process_inplace('filter.bandpass.gauss',{'apix':apix,'center':0.0,'cutoff_abs':0.005,'sigma':sigma})
					ptcl.process_inplace('mask.soft',{'outer_radius': -0.20 * ptcl['nx']})

			if not options.noedgenorm: ptcl.process_inplace('normalize.edgemean')
			if not options.nobinarize: ptcl.process_inplace('threshold.binary',{'value':ptcl.get_attr('sigma') * options.deviations})

			if not options.nowrite: ptcl.write_image(options.path[:-4]+'_proc.hdf',-1)

		if options.postdisplay: display([EMData(options.output,i) for i in xrange(nimgs)])

	theta = 2*np.pi/360 * np.array(range(0, 360, options.step) + [0])

	with PdfPages(options.path[:-4] + '.pdf') as pdf:
		for i in xrange(nimgs):
			if options.verbose > 5: print("Centering image {}/{}".format(i+1,nimgs))

			orig = EMData(options.path,i)

			if not options.nopreprocess: ptcl = EMData(options.path[:-4]+'_proc.hdf',i)
			else: ptcl = orig
			parr = ptcl.numpy()

			fig = plt.figure(figsize=(10,10))

			ax0 = fig.add_subplot(221)
			img = ax0.imshow(orig.numpy(),cmap=plt.cm.Greys_r)
			#fig.colorbar(img, orientation='horizontal')
			ax0.invert_yaxis()
			ax0.xaxis.set_ticklabels([])
			ax0.yaxis.set_ticklabels([])
			ax0.set_title('Original')

			eman_centered = orig.process('xform.center')
			t = Transform({'type':'2d','alpha':270.0})
			eman_centered.transform(t)
			ax2 = fig.add_subplot(222)
			ax2.imshow(eman_centered.numpy().transpose(),cmap=plt.cm.Greys_r) # mu
			ax2.set_title('Centered (EMAN)')

			ax1 = fig.add_subplot(223)
			proc = ax1.imshow(parr,cmap=plt.cm.Greys_r)
			#fig.colorbar(proc, orientation='horizontal')
			ax1.invert_yaxis()
			ax1.xaxis.set_ticklabels([])
			ax1.yaxis.set_ticklabels([])
			ax1.set_title('Binarized')

			ax3 = fig.add_subplot(224)
			centered = orig.copy()

			#mu = np.mean(parr.nonzero(),axis=1) # rough geometric center
			mu = (np.max(parr.nonzero(),axis=1) + np.min(parr.nonzero(),axis=1))/2.

			tx = -(mu[1] - nx)
			ty = -(mu[0] - nx)
			if options.verbose > 8: print("Translating image {} by ({},{})".format(i,tx,ty))
			t = Transform({'type':'eman','tx':tx,'ty':ty})
			centered.transform(t)
			final = ax3.imshow(centered.numpy(),cmap=plt.cm.Greys_r)
			ax3.invert_yaxis()
			#fig.colorbar(final, orientation='horizontal')
			ax3.xaxis.set_ticklabels([])
			ax3.yaxis.set_ticklabels([])
			ax3.set_title('Centered (New)')

			# r = []
			# s = []
			# for i, angle in enumerate(theta):
			# 	proj = ptcl.project("standard",Transform({'type':'2d','alpha':angle}))
			# 	projarr = proj.numpy()
			# 	projarr = projarr/np.max(projarr)
			# 	nz = np.nonzero(projarr)
			# 	mu = np.mean(nz)
			# 	sd = np.std(nz)
			# 	r.append(mu)
			# 	s.append(sd)
			# r = np.array(r/np.max(r))
			# s = np.array(s)
			#
			# ax2 = fig.add_subplot(223,polar=True)
			# ax2.plot(theta, r, "ro") # mu
			# ax2.errorbar(theta, r, yerr=s, xerr=0, capsize=0) #sigma
			# ax2.set_ylim([0,1])
			# ax2.set_title('Projections')

			pdf.savefig(fig)

			plt.close()

			if not options.nowrite: ptcl.write_image(options.output,-1)

if __name__ == "__main__":
	main()
