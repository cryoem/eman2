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
	
	test_dims = [64*i for i in [1,2,3,4]]
	test_dims_3d = [64*i for i in [1,2,3,4]]
	
	test_range = range(30)
	
	gpu_times = []
	cpu_times = []
	
	
	print "Testing 180 rotation"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		for d in a: d.set_gpu_rw_current()
		t = time()
		for i in test_range:
			a[i].process_inplace("math.rotate.180")
			
		gpu_times.append(time()-t)
		
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		t = time()
		for i in test_range:
			a[i].process_inplace("math.rotate.180")
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
	
#	print "Testing phase comparison"
#	for dims in test_dims:
#		a = [test_image(0,size=(dims,dims)) for i in test_range]
#		for d in a: d.set_gpu_rw_current()
#		tr = Transform()
#		
#		t = time()
#		for i in test_range:
#			b = a[i].cmp("phase",a[i])
#			#c.print_this()
#
#		gpu_times.append(time()-t)
#		#print dims, "B"
#		a = [test_image(0,size=(dims,dims)) for i in test_range]
#		t = time()
#		for i in test_range:
#			b = a[i].cmp("phase",a[i])
#		cpu_times.append(time()-t)
#		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
#	
#	alis = ["translational","rotational","rotate_translate","rotate_translate_flip"]
#	for ali in alis:
#		print "Testing",ali,"alignment"
#		print "Dims","\t", "GPU speedup"
#		for dims in test_dims:
#			a = [test_image(0,size=(dims,dims)) for i in test_range]
#			b = [test_image(0,size=(dims,dims)) for i in test_range]
#			for d in a: d.set_gpu_rw_current()
#			for d in b: d.set_gpu_rw_current()
#			
#			t = time()
#			for i in test_range:
#				a[i].set_gpu_rw_current()
#				b[i].set_gpu_rw_current()
#				c = a[i].align(ali,b[i],{},"phase",{})
#				#c.print_this()
#	
#			gpu_times.append(time()-t)
#			#print dims, "B"
#			a = [test_image(0,size=(dims,dims)) for i in test_range]
#			b = [test_image(0,size=(dims,dims)) for i in test_range]
#			t = time()
#			for i in test_range:
#				c = a[i].align(ali,b[i],{},"phase",{})
#			cpu_times.append(time()-t)
#			print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
	
	print "Testing transform (2D)"
	for dims in test_dims:
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		for d in a: d.set_gpu_rw_current()
		tr = Transform()
		
		t = time()
		for i in test_range:
			b = a[i].process("math.transform",{"transform":tr})
			#c.print_this()

		gpu_times.append(time()-t)
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		t = time()
		for i in test_range:
			b = a[i].process("math.transform",{"transform":tr})
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
		
	print "Testing transform (3D)"
	for dims in [64,96,128,256]:
		a = [test_image_3d(0,size=(dims,dims,dims)) for i in range(1)]
		for d in a: d.set_gpu_rw_current()
		tr = Transform()
		
		t = time()
		for i in range(10):
			b = a[0].process("math.transform",{"transform":tr})
			#c.print_this()

		gpu_times.append(time()-t)
		#print dims, "B"
		a = [test_image_3d(0,size=(dims,dims,dims)) for i in range(1)]
		t = time()
		for i in range(10):
			b = a[0].process("math.transform",{"transform":tr})
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
	
	print "Testing calc_ccfx"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		for i in a: i.set_gpu_rw_current()
		t = time()
		for i in test_range:
			c = a[i].calc_ccfx(a[i],0,-1,True)
			#c.print_this()

		gpu_times.append(time()-t)
		#print dims, "B"
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		#a.print_this()
		t = time()
		for i in test_range:
			c = a[i].calc_ccfx(a[i],0,-1,True)
			#print "the other stuff is ", a.get_cuda_handle()
			#a.print_this()
			#b.print_this()
			#c.print_this()
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
	
	print "Testing calc_ccfx column sum"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		for i in a: i.set_gpu_rw_current()
		t = time()
		for i in test_range:
			c = a[i].calc_ccfx(a[i])

		gpu_times.append(time()-t)
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		t = time()
		for i in test_range:
			c = a[i].calc_ccfx(a[i])

		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
	
	print "Testing get_clip, lose edge pixels"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		for i in a: i.set_gpu_rw_current()
		r = Region(1,1,dims-2,dims-2)
		t = time()
		for i in test_range:
			b = a[i].get_clip(r)
		gpu_times.append(time()-t)
		
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		t = time()
		for i in test_range:
			b = a[i].get_clip(r)
			
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
		
	
	
	print "Testing phase origin"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims))
		a.set_gpu_rw_current()
		t = time()
		for i in test_range:
			a.process_inplace("xform.phaseorigin.tocenter")
			
		gpu_times.append(time()-t)
		
		a = test_image(0,size=(dims,dims))
		t = time()
		for i in test_range:
			a.process_inplace("xform.phaseorigin.tocenter")
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
	
	
	print "Testing make rotational footprint"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		for d in a:
			d.set_gpu_rw_current()
		t = time()
		for i in test_range:
			a[i].set_gpu_rw_current()
			c = a[i].make_rotational_footprint_cuda()
				
		gpu_times.append(time()-t)
		a = [test_image(0,size=(dims,dims)) for i in test_range]
		t = time()
		for i in test_range:
			c = a[i].make_rotational_footprint_e1()
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
	
	print "Testing unwrap"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims))
		a.set_gpu_rw_current()
		
		t = time()
		for i in test_range:
			c = a.unwrap()
			
		gpu_times.append(time()-t)
		
		a = test_image(0,size=(dims,dims))
		
		for i in test_range:
			c = a.unwrap()
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1],'\t',cpu_times[-1],'\t',gpu_times[-1]
	
	print "Testing pixel multiplication by a constant (3D)"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims_3d:
		a = test_image_3d(0,size=(dims,dims,dims))
		a.set_gpu_rw_current()
		t = time()
		for i in test_range:
			a.mult(2.0)
			a.mult(0.5)
		gpu_time = time()-t
		
		a = test_image_3d(0,size=(dims,dims,dims))
		t = time()
		for i in test_range:
			a.mult(2.0)
			a.mult(0.5)
		cpu_time = time()-t
		print dims,"\t", cpu_time/gpu_time,'\t',cpu_time,'\t',gpu_time
	
	print "Testing pixel multiplication by a constant (2D)"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims))
		a.set_gpu_rw_current()
		t = time()
		for i in test_range:
			a.mult(2.0)
			a.mult(0.5)
		gpu_time = time()-t
		
		a = test_image(0,size=(dims,dims))
		t = time()
		for i in test_range:
			a.mult(2.0)
			a.mult(0.5)
		cpu_time = time()-t
		print dims,"\t", cpu_time/gpu_time,'\t',cpu_time,'\t',gpu_time

	print "Testing FFT (2D)"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims))
		#tmp = a.do_fft_cuda()
		t = time()
		for i in test_range:
			#a = aa.copy()
			b = a.do_fft_cuda()
		gpu_times.append(time()-t)
		a = test_image(0,size=(dims,dims))
		t = time()
		for i in test_range:
			#a = aa.copy()
			b = a.do_fft()
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1]
		
	print "Testing IFT (2D)"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims))
		b = a.do_fft_cuda()
		#tmp = a.do_fft_cuda()
		t = time()
		for i in test_range:
			c = b.do_ift_cuda()
		gpu_times.append(time()-t)
		a = test_image(0,size=(dims,dims))
		b =a.do_fft()
		t = time()
		for i in test_range:
			c = b.do_ift()
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

	print "Testing Fourier correlation (2D) NO TEXTURE"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims)).do_fft_cuda()
		b = test_image(0,size=(dims,dims)).do_fft_cuda()
		
		t = time()
		for i in test_range:
			c = a.calc_ccf_cuda(b,False,1)
			
		gpu_times.append(time()-t)
		a = test_image(0,size=(dims,dims)).do_fft()
		b = test_image(0,size=(dims,dims)).do_fft()
		t = time()
		for i in test_range:
			c = a.calc_ccf(b,fp_flag.CIRCULANT,1)
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1]
		
	#print "Testing Fourier correlation (2D) NO TEXTURE - 853x1272"
	#print "Dims","\t", "GPU speedup"
	#for dims in [[853,1272]]:
		#a = test_image(0,size=(dims[0],dims[1])).do_fft_cuda()
		#b = test_image(0,size=(dims[0],dims[1])).do_fft_cuda()
		
		#t = time()
		#for i in test_range:
			#c = a.calc_ccf_cuda(b,True)
			
		#gpu_times.append(time()-t)
		#a = test_image(0,size=(dims[0],dims[1])).do_fft()
		#b = test_image(0,size=(dims[0],dims[1])).do_fft()
		#t = time()
		#for i in test_range:
			#c = a.calc_ccf(b)
		#cpu_times.append(time()-t)
		#print dims[0],'x',dims[1],"\t", cpu_times[-1]/gpu_times[-1]
		
	print "Testing Fourier correlation (2D) TEXTURE"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims:
		a = test_image(0,size=(dims,dims)).do_fft_cuda()
		b = test_image(0,size=(dims,dims)).do_fft_cuda()
		
		t = time()
		for i in test_range:
			c = a.calc_ccf_cuda(b,True,1)
			
		gpu_times.append(time()-t)
		a = test_image(0,size=(dims,dims)).do_fft()
		b = test_image(0,size=(dims,dims)).do_fft()
		t = time()
		for i in test_range:
			c = a.calc_ccf(b,fp_flag.CIRCULANT,1)
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1]

	test_range = range(10)
	print "Testing Fourier correlation (3D) NO TEXTURE"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims_3d:
		a = test_image_3d(0,size=(dims,dims,dims)).do_fft_cuda()
		b = test_image_3d(0,size=(dims,dims,dims)).do_fft_cuda()
		
		t = time()
		for i in test_range:
			c = a.calc_ccf_cuda(b,False,1)
			
		gpu_times.append(time()-t)
		a = test_image_3d(0,size=(dims,dims,dims)).do_fft()
		b = test_image_3d(0,size=(dims,dims,dims)).do_fft()
		t = time()
		for i in test_range:
			c = a.calc_ccf(b,fp_flag.CIRCULANT,1)
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1]
		
	#print "Testing Fourier correlation (3D) TEXTURE"
	#print "Dims","\t", "GPU speedup"
	#for dims in test_dims_3d:
		#a = test_image_3d(0,size=(dims,dims,dims)).do_fft_cuda()
		#b = test_image_3d(0,size=(dims,dims,dims)).do_fft_cuda()
		
		#t = time()
		#for i in test_range:
			#c = a.calc_ccf_cuda(b,True)
			
		#gpu_times.append(time()-t)
		#a = test_image_3d(0,size=(dims,dims,dims)).do_fft()
		#b = test_image_3d(0,size=(dims,dims,dims)).do_fft()
		#t = time()
		#for i in test_range:
			#c = a.calc_ccf(b)
		#cpu_times.append(time()-t)
		#print dims,"\t", cpu_times[-1]/gpu_times[-1]

	
	print "Testing 3D real space projection"
	print "Dims","\t", "GPU speedup"
	for dims in test_dims_3d:
		a = test_image_3d(5,size=(dims,dims,dims))
		trans = Transform()
		t = time()
		for i in test_range:
			p = a.project("cuda_standard",trans)
		gpu_times.append(time()-t)
		
		a = test_image_3d(5,size=(dims,dims,dims))
		t = time()
		for i in test_range:
			p = a.project("standard",trans)
		cpu_times.append(time()-t)
		print dims,"\t", cpu_times[-1]/gpu_times[-1]

if __name__ == '__main__':
	if EMUtil.cuda_available(): test_main()