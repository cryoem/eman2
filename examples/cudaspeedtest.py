#!/usr/bin/env python

#
# Author: David Woolford March 2nd 2009
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
#


import EMAN2
from EMAN2 import *
import sys
import math
import os
from time import time

def test_main():
	
	test_dims = [64*i for i in [2,3,4,6,8]]
	
	test_range = range(20)
	
	gpu_times = []
	cpu_times = []
	
	print "Testing pixel multiplication by a constant (2D)"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims))
		a._copy_cpu_to_gpu_rw()
		t = time()
		for i in test_range:
			a.mult_cuda(2.0)
			a.mult_cuda(0.5)
		gpu_times.append(time()-t)
		
		a = test_image(0,size=(dims,dims))
		t = time()
		for i in test_range:
			a.mult(2.0)
			a.mult(0.5)
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1]
	
	
	print "Testing FFT/IFT (2D)"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims))
		tmp = a.do_fft_cuda()
		t = time()
		for i in test_range:
			#a = aa.copy()
			b = a.do_fft_cuda()
			c = b.do_ift_cuda()
		gpu_times.append(time()-t)
		a = test_image(0,size=(dims,dims))
		t = time()
		for i in test_range:
			#a = aa.copy()
			b = a.do_fft()
			c = b.do_ift()
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1]

	print "Testing Fourier correlation (2D)"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims)).do_fft_cuda()
		b = test_image(0,size=(dims,dims)).do_fft_cuda()
		
		t = time()
		for i in test_range:
			c = a.calc_ccf_cuda(b)
		gpu_times.append(time()-t)
		a = test_image(0,size=(dims,dims)).do_fft()
		b = test_image(0,size=(dims,dims)).do_fft()
		t = time()
		for i in test_range:
			c = a.calc_ccf(b)
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1]


	test_range = range(10)
	print "Testing 3D real space projection"
	print "Dims","\t", "GPU speedup"
	#for dims in test_dims:
		#a = test_image_3d(5,size=(dims,dims,dims))
		#trans = Transform()
		#t = time()
		#for i in test_range:
			#p = a.project("cuda_standard",trans)
		#gpu_times.append(time()-t)
		
		#a = test_image_3d(5,size=(dims,dims,dims))
		#t = time()
		#for i in test_range:
			#p = a.project("standard",trans)
		#cpu_times.append(time()-t)
		#print dims,"\t", cpu_times[-1]/gpu_times[-1]

if __name__ == '__main__':
	if EMUtil.cuda_available(): test_main()